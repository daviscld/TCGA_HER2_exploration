# ============================================================
# Multimodal HER2 Definition: RNA vs DNA vs Combined
# Modalities: ERBB2 mRNA (VST) | GISTIC2 CNV | FGA | Mutation count
# Outputs: ROC/AUC, logistic coefficients, SHAP, concordance plots
# ============================================================
# Expected inputs (from previous QC + DESeq2 scripts):
#   metadata_qc_passed.csv  – must contain: her2_clinical (HER2pos/HER2neg),
#                             fraction_genome_altered, mutation_count,
#                               ideally also patient_id for alignment 
#                               (if not, rownames will be used),
#                               and also ERBB2 copy number column, -2 to 2
#   vst_her2_brca.rds       – VST SummarizedExperiment from DESeq2 script

# ============================================================

library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pROC)           # ROC / AUC
library(caret)          # cross-validation
library(glmnet)         # logistic regression with regularisation
library(randomForest)
library(shapviz)        # SHAP for R
library(kernelshap)     # model-agnostic SHAP backend
library(pheatmap)
library(RColorBrewer)
library(yardstick)      # tidy metrics (PR AUC)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)

out_dir <- "multimodal_plots";   dir.create(out_dir, showWarnings = FALSE)
res_dir <- "multimodal_results"; dir.create(res_dir, showWarnings = FALSE)

cv_folds     <- 5
random_seed  <- 42
erbb2_gene   <- "ERBB2"
palette      <- c("HER2pos" = "#E64B35", "HER2neg" = "#4DBBD5")

# Helper: save plots consistently
save_plot <- function(p, fname, w=8, h=6, dpi=300) {
  ggsave(file.path(out_dir, fname), p, width=w, height=h, dpi=dpi,
         bg = "white")
  invisible(p)
}

# ═══════════════════════════════════════════════════════════
# 1. LOAD & ALIGN DATA
# ═══════════════════════════════════════════════════════════

meta <- read.csv("metadata_qc_passed.csv", row.names = 1)
vst  <- readRDS("vst_her2_brca.rds")            # RangedSummarizedExperiment
meta_ids <- if ("patient_id" %in% colnames(meta)) meta$patient_id else rownames(meta)
cnv <- meta$erbb2_copy_number %>%
  setNames(meta_ids) %>%
  as.data.frame()
colnames(cnv) <- "ERBB2_CN"

# Keep only samples with definitive clinical FISH/IHC label
meta <- meta[meta$her2_clinical %in% c("HER2pos", "HER2neg"), ]
meta$her2_binary <- as.integer(meta$her2_clinical == "HER2pos")
meta_model_ids <- if ("patient_id" %in% colnames(meta)) meta$patient_id else rownames(meta)

# ═══════════════════════════════════════════════════════════
# 2. BUILD FEATURE MATRICES
# ═══════════════════════════════════════════════════════════

vst_mat <- assay(vst)   # genes × samples

# ── 2A. RNA: ERBB2 mRNA ───────────────────────────────────
erbb2_row <- intersect(erbb2_gene, rownames(vst_mat))
if (length(erbb2_row) == 0) stop("ERBB2 gene symbol not found in VST matrix")

X_rna_erbb2 <- data.frame(
  ERBB2_mRNA = as.numeric(vst_mat[erbb2_row[1], ]),
  row.names  = colnames(vst_mat)
)

# ── 2B. RNA: Expanded — ERBB2 amplicon neighbours ─────────
# GRB7, STARD3, MIEN1 are co-amplified on chr17q12; biologically meaningful

#doi 10.1016/j.molonc.2012.10.012
#The smallest common region of amplification found in all of the 
#71 tumors analyzed was 78.61 Kbp, including six genes; 
#STARD3, TCAP, PNMT, PERLD1, HER2, and C17orf37. 
#Ninety percent (64/71) of the tumors shared a 255.74 Kbp 
#amplification region consisting of ten genes; 
#NEUROD2, PPP1R1B, STARD3, TCAP, PNMT, PERLD1, HER2, C17orf37, GRB7 and ZNFN1A3
small_amplicon_genes <- c("ERBB2", "TCAP", "PNMT", "PERLD1", "C17orf37")
larger_amplicon_genes <- c(small_amplicon_genes, c("NEUROD2", "PPP1R1B", 
                                "STARD3", "TCAP", 
                                "PNMT", "PERLD1",
                            "ERBB2", "C17orf37", "GRB7", "ZNFN1A3"))    

amp_rows <- intersect(larger_amplicon_genes, rownames(vst_mat))
if (length(amp_rows) == 0) stop("None of the amplicon gene symbols were found in VST matrix")

X_rna_amplicon <- t(vst_mat[amp_rows, ]) %>%
  as.data.frame() %>%
  setNames(colnames(.))

# ── 2C. DNA: ERBB2 score ──────────────────────────────────
# CNV score: -2=deep deletion, -1=shallow del, 0=neutral, 1=gain, 2=amp
cnv_col <- intersect(c("ERBB2", "ERBB2_CN"), colnames(cnv))[1]
if (is.na(cnv_col)) stop("Cannot find ERBB2 column in CNV file")

X_cnv <- data.frame(
  ERBB2_score = cnv[[cnv_col]],
  row.names    = rownames(cnv)
)

# Binary amplification flag (score >= 1 = gain/amplification)
X_cnv$ERBB2_amp_flag <- as.integer(X_cnv$ERBB2_score >= 1)

# ── 2D. Genomic instability features ─────────────────────
X_genomic <- meta %>%
  dplyr::select(Fraction_Genome_Altered, Mutation_Count) %>%
  mutate(
    fga_scaled  = scale(Fraction_Genome_Altered)[, 1],
    mut_log     = log10(Mutation_Count + 1)
  ) %>%
  dplyr::select(fga_scaled, mut_log)
rownames(X_genomic) <- meta_model_ids

combine_by_sample <- function(...) {
  x_list <- list(...)
  x_tbls <- lapply(x_list, function(df) {
    tibble::rownames_to_column(as.data.frame(df), var = "sample_id")
  })
  out <- Reduce(function(a, b) full_join(a, b, by = "sample_id"), x_tbls)
  out <- as.data.frame(out)
  rownames(out) <- out$sample_id
  out$sample_id <- NULL
  out
}

# ── 2E. Combined feature matrices ─────────────────────────
# M4 (glmnet) and M5 (RF) both use the full amplicon + CNV + genomic set.
# glmnet's regularisation handles collinearity between amplicon neighbours;
# using the same features makes M4 vs M5 a pure model-class comparison.
X_combined <- combine_by_sample(
  X_rna_amplicon,   # ERBB2 + all chr17q12 amplicon neighbours
  X_cnv,
  X_genomic
) %>% as.data.frame()

# X_rf_all is identical to X_combined — kept as a named alias for clarity
X_rf_all <- X_combined

# Outcome vector (aligned)
y        <- meta$her2_binary
y_factor <- factor(ifelse(y == 1, "HER2pos", "HER2neg"),
                   levels = c("HER2neg", "HER2pos"))
y_factor_by_id <- setNames(y_factor, meta_model_ids)

# ── 2F. Quick plot: ERBB2 mRNA by copy number ─────────────
erbb2_by_cn <- data.frame(
  ERBB2_mRNA = X_rna_erbb2$ERBB2_mRNA,
  ERBB2_CN   = X_cnv$ERBB2_score,
  row.names  = rownames(meta)
) %>%
  filter(!is.na(ERBB2_CN))

my_comparisons <- list(
  "0 vs 1" = c("0", "1"),
  "0 vs 2" = c("0", "2"),
  "1 vs 2" = c("1", "2"),
  "-1 vs 0" = c("-1", "0"),
  "-1 vs 1" = c("-1", "1"),
  "-1 vs 2" = c("-1", "2")
)

p_mrna_by_cn <- ggplot(erbb2_by_cn, aes(factor(ERBB2_CN), ERBB2_mRNA)) +
  geom_violin(fill = "#3C5488", alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.size = 1) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  stat_compare_means(method = "kruskal.test",
                     label.y = max(erbb2_by_cn$ERBB2_mRNA) + 0.5) +
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons,
                     label = "p.signif",
                     p.adjust.method = "BH",
                     tip.length = 0.01,
                     step.increase = 0.1) +
  labs(title    = "ERBB2 mRNA Expression by Copy Number",
       x = "ERBB2 Copy Number",
       y = "ERBB2 VST Expression") +
  theme_bw(base_size = 11)

save_plot(p_mrna_by_cn, "00_ERBB2_mRNA_by_CN.png", w = 7, h = 5)

# ═══════════════════════════════════════════════════════════
# 3. CROSS-VALIDATED MODEL EVALUATION
# ═══════════════════════════════════════════════════════════
# M1: ERBB2 mRNA alone
# M2: ERBB2 CNV score alone
# M3: ERBB2 CNV + amplification flag
# M4: Amplicon + CNV + genomic (glmnet) ─┐ same features,
# M5: Amplicon + CNV + genomic (RF)      ─┘ different model class

set.seed(random_seed)

# ── M6: PLS-DA (Partial Least Squares Discriminant Analysis) ──
# Uses same features as M4/M5 for direct comparison

if (!requireNamespace("pls", quietly = TRUE)) install.packages("pls")
library(pls)

fit_plsda <- function(X, label) {
  model_ids <- intersect(rownames(X), names(y_factor_by_id))
  df <- as.data.frame(X[model_ids, , drop = FALSE])
  df$outcome <- y_factor_by_id[model_ids]
  df <- df[complete.cases(df), , drop = FALSE]
  
  if (nrow(df) < cv_folds) stop(label, ": too few complete-case samples")
  if (length(unique(df$outcome)) < 2) stop(label, ": only one outcome class")
  
  fold_idx <- createFolds(df$outcome, k = cv_folds, returnTrain = TRUE)
  ctrl <- trainControl(
    method = "cv", number = cv_folds,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", sampling = "up",
    index = fold_idx,
    indexOut = lapply(fold_idx, function(i) setdiff(seq_len(nrow(df)), i))
  )
  
  fit <- train(outcome ~ ., data = df, method = "pls",
               trControl = ctrl, metric = "ROC",
               tuneGrid = expand.grid(ncomp = 1:5))
  
  attr(fit, "n_samples") <- nrow(df)
  attr(fit, "n_her2_pos") <- sum(df$outcome == "HER2pos")
  attr(fit, "n_her2_neg") <- sum(df$outcome == "HER2neg")
  fit
}

cat("\nFitting PLS-DA model...\n")
m6_plsda <- fit_plsda(X_combined, "M6: Amplicon + CNV + genomic (PLS-DA)")

models[["M6: Amplicon + CNV + genomic (PLS-DA)"]] <- m6_plsda

# ── Helpers ───────────────────────────────────────────────
roc_auc_score_r <- function(obs, prob) {
  as.numeric(roc(as.integer(obs == "HER2pos"), prob, quiet = TRUE)$auc)
}

pr_auc_score_r <- function(obs, prob) {
  yardstick::pr_auc(
    data.frame(truth = factor(obs, levels = c("HER2neg", "HER2pos")),
               prob = prob),
    truth, prob, event_level = "second"
  )$.estimate
}

# Compute all metrics from a fold or aggregated predictions data frame.
# Expects columns: obs (factor), HER2pos (numeric probability)
calc_metrics <- function(df) {
  pred_class <- factor(ifelse(df$HER2pos >= 0.5, "HER2pos", "HER2neg"),
                       levels = c("HER2neg", "HER2pos"))
  obs_class  <- factor(df$obs, levels = c("HER2neg", "HER2pos"))
  cm <- confusionMatrix(pred_class, obs_class, positive = "HER2pos")
  data.frame(
    ROC_AUC           = roc_auc_score_r(df$obs, df$HER2pos),
    PR_AUC            = pr_auc_score_r(df$obs, df$HER2pos),
    Accuracy          = as.numeric(cm$overall["Accuracy"]),
    Balanced_Accuracy = as.numeric(cm$byClass["Balanced Accuracy"]),
    Sensitivity       = as.numeric(cm$byClass["Sensitivity"]),
    Specificity       = as.numeric(cm$byClass["Specificity"]),
    F1                = as.numeric(cm$byClass["F1"])
  )
}

metric_names <- c("ROC_AUC","PR_AUC","Accuracy","Balanced_Accuracy",
                  "Sensitivity","Specificity","F1")

# ── Model fitting ─────────────────────────────────────────
fit_model <- function(X, label, method = "glm") {
  model_ids <- intersect(rownames(X), names(y_factor_by_id))
  df <- as.data.frame(X[model_ids, , drop = FALSE])
  df$outcome <- y_factor_by_id[model_ids]
  df <- df[complete.cases(df), ]

  if (nrow(df) < cv_folds)          stop(label, ": too few complete-case samples")
  if (length(unique(df$outcome)) < 2) stop(label, ": only one outcome class")

  # Shared folds: real and shuffled models see identical splits
  fold_idx <- createFolds(df$outcome, k = cv_folds, returnTrain = TRUE)
  ctrl <- trainControl(
    method = "cv", number = cv_folds,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", sampling = "up",  # upsample — preserves all HER2+ cases
    index = fold_idx,
    indexOut = lapply(fold_idx, function(i) setdiff(seq_len(nrow(df)), i))
  )

  fit          <- train(outcome ~ ., data = df,         method = method,
                        trControl = ctrl, metric = "ROC",
                        preProcess = c("center","scale"))
  fit_shuffled <- train(outcome ~ ., data = {df$outcome <- sample(df$outcome); df},
                        method = method, trControl = ctrl, metric = "ROC",
                        preProcess = c("center","scale"))

  attr(fit, "n_samples")    <- nrow(df)
  attr(fit, "n_her2_pos")   <- sum(df$outcome == "HER2pos")
  attr(fit, "n_her2_neg")   <- sum(df$outcome == "HER2neg")
  attr(fit, "shuffled_fit") <- fit_shuffled
  fit
}

cat("\nFitting models...\n")
m1_rna      <- fit_model(X_rna_erbb2,          "M1: ERBB2 mRNA")
m2_cnv      <- fit_model(X_cnv["ERBB2_score"], "M2: ERBB2 CNV")
m3_cnv_bin  <- fit_model(X_cnv,                "M3: ERBB2 CNV + amp flag")
m4_combined <- fit_model(X_combined,            "M4: Amplicon + CNV + genomic (glmnet)", method = "glmnet")
m5_rf       <- fit_model(X_rf_all,             "M5: Amplicon + CNV + genomic (RF)",     method = "rf")

models <- list(
  "M1: ERBB2 mRNA"                        = m1_rna,
  "M2: ERBB2 CNV"                         = m2_cnv,
  "M3: ERBB2 CNV + amp flag"              = m3_cnv_bin,
  "M4: Amplicon + CNV + genomic (glmnet)" = m4_combined,
  "M5: Amplicon + CNV + genomic (RF)"     = m5_rf
)

# ── Sample size table ─────────────────────────────────────
model_sample_sizes <- data.frame(
  Model     = names(models),
  N         = vapply(models, \(m) attr(m, "n_samples"),  numeric(1)),
  N_HER2pos = vapply(models, \(m) attr(m, "n_her2_pos"), numeric(1)),
  N_HER2neg = vapply(models, \(m) attr(m, "n_her2_neg"), numeric(1))
)
print(model_sample_sizes)
write.csv(model_sample_sizes, file.path(res_dir, "model_sample_sizes.csv"), row.names = FALSE)

# ── Aggregate CV metrics ──────────────────────────────────
cv_metrics <- lapply(names(models), function(nm) {
  preds <- models[[nm]]$pred %>%
    group_by(rowIndex) %>%
    summarise(HER2pos = mean(HER2pos), obs = first(obs), .groups = "drop")
  cbind(Model = nm, calc_metrics(preds))
}) %>% bind_rows()

# ── Per-fold real vs shuffled (paired Wilcoxon) ───────────
cv_real_vs_shuffled <- lapply(names(models), function(nm) {
  real_folds     <- models[[nm]]$pred             %>% group_by(Resample) %>% group_modify(~ calc_metrics(.x)) %>% ungroup()
  shuffled_folds <- attr(models[[nm]], "shuffled_fit")$pred %>% group_by(Resample) %>% group_modify(~ calc_metrics(.x)) %>% ungroup()

  lapply(metric_names, function(metric) {
    rv <- real_folds[[metric]];  sv <- shuffled_folds[[metric]]
    p  <- tryCatch(
      wilcox.test(rv, sv, paired = TRUE, alternative = "greater", exact = FALSE)$p.value,
      error = \(e) NA_real_
    )
    data.frame(Model = nm, Metric = metric,
               Real_Mean = mean(rv, na.rm=TRUE), Shuffled_Mean = mean(sv, na.rm=TRUE),
               Delta = mean(rv - sv, na.rm=TRUE), P_Value = p)
  }) %>% bind_rows()
}) %>% bind_rows()

cv_real_vs_shuffled$Adj_P_Value <- p.adjust(cv_real_vs_shuffled$P_Value, method = "BH")

cv_metrics_significant <- cv_real_vs_shuffled %>%
  filter(!is.na(Adj_P_Value), Adj_P_Value < 0.05, Delta > 0) %>%
  arrange(Model, Adj_P_Value, desc(Delta))

# ── Save & report ─────────────────────────────────────────
print(cv_metrics)
write.csv(cv_metrics,             file.path(res_dir, "model_cv_metrics.csv"),                          row.names = FALSE)
write.csv(cv_real_vs_shuffled,    file.path(res_dir, "model_cv_metrics_real_vs_shuffled_all.csv"),     row.names = FALSE)
write.csv(cv_metrics_significant, file.path(res_dir, "model_cv_metrics_significant_vs_shuffled.csv"), row.names = FALSE)

cat("\nMetrics significantly better than shuffled labels (BH-adjusted p < 0.05):\n")
if (nrow(cv_metrics_significant) == 0) cat("  None — consider stress-testing on IHC 2+ equivocal subset.\n") else print(cv_metrics_significant)

# ═══════════════════════════════════════════════════════════
# 4. ROC CURVES — all models overlaid
# ═══════════════════════════════════════════════════════════
roc_list <- lapply(names(models), function(nm) {
  preds <- models[[nm]]$pred %>%
    group_by(rowIndex) %>%
    summarise(prob = mean(HER2pos), obs = first(obs), .groups = "drop")
  roc_obj <- roc(as.integer(preds$obs == "HER2pos"), preds$prob, quiet = TRUE)
  data.frame(
    FPR   = 1 - roc_obj$specificities,
    TPR   = roc_obj$sensitivities,
    Model = nm,
    AUC   = round(as.numeric(roc_obj$auc), 3)
  )
}) %>% bind_rows() %>%
  mutate(Model_AUC = paste0(Model, "  (AUC=", AUC, ")"))

p_roc <- ggplot(roc_list, aes(FPR, TPR, color = Model_AUC)) +
  geom_line(linewidth = 0.9) +
  geom_abline(linetype = "dashed", color = "grey60") +
  scale_color_brewer(palette = "Set1") +
  labs(title    = "ROC Curves — HER2+ Prediction by Modality",
       subtitle  = paste0(cv_folds, "-fold cross-validated"),
       x = "False Positive Rate", y = "True Positive Rate",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.62, 0.22),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8))
save_plot(p_roc, "01_ROC_curves.png", w = 8, h = 6)

# ── Precision-Recall curves ───────────────────────────────
if (!requireNamespace("PRROC", quietly = TRUE)) install.packages("PRROC")
library(PRROC)

pr_list <- lapply(names(models), function(nm) {
  preds <- models[[nm]]$pred %>%
    group_by(rowIndex) %>%
    summarise(prob = mean(HER2pos), obs = first(obs), .groups = "drop")
  pr <- PRROC::pr.curve(
    scores.class0 = preds$prob[preds$obs == "HER2pos"],
    scores.class1 = preds$prob[preds$obs == "HER2neg"],
    curve = TRUE)
  data.frame(Recall = pr$curve[, 1], Precision = pr$curve[, 2],
             Model = nm, PRAUC = round(pr$auc.integral, 3))
}) %>% bind_rows() %>%
  mutate(Model_AUC = paste0(Model, "  (PRAUC=", PRAUC, ")"))

baseline_prev <- mean(models[["M4: Amplicon + CNV + genomic (glmnet)"]]$trainingData$.outcome == "HER2pos")

p_pr <- ggplot(pr_list, aes(Recall, Precision, color = Model_AUC)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = baseline_prev, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.8, y = baseline_prev + 0.02,
           label = paste0("Baseline (prevalence=", round(baseline_prev, 2), ")"),
           size = 3, color = "grey40") +
  scale_color_brewer(palette = "Set1") +
  labs(title    = "Precision-Recall Curves — HER2+ Prediction",
       subtitle  = paste0(cv_folds, "-fold cross-validated"),
       x = "Recall", y = "Precision", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.38, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8))
save_plot(p_pr, "02_PR_curves.png", w = 8, h = 6)

# ── AUC bar chart ─────────────────────────────────────────
p_auc_bar <- cv_metrics %>%
  pivot_longer(c(ROC_AUC, PR_AUC), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = reorder(Model, Value), y = Value, fill = Metric)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = round(Value, 3)),
            position = position_dodge(0.65), hjust = -0.1, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("ROC_AUC" = "#3C5488", "PR_AUC" = "#E64B35")) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(title = "Cross-validated AUC by Modality",
       x = NULL, y = "AUC", fill = NULL) +
  theme_bw(base_size = 11)
save_plot(p_auc_bar, "03_AUC_comparison.png", w = 9, h = 5)

# ═══════════════════════════════════════════════════════════
# 5. LOGISTIC REGRESSION COEFFICIENTS
# ═══════════════════════════════════════════════════════════
# Refit glmnet on full data for coefficient interpretation.
# Uses the same X_combined as M4 so coefficients are directly
# interpretable in the context of the CV-validated model.

prep_glmnet <- function(X) {
  scale(as.matrix(X))   # z-score so coefficients are comparable
}

fit_lr_full <- function(X, label) {
  model_ids <- intersect(rownames(X), names(y_factor_by_id))
  df <- as.data.frame(X[model_ids, , drop = FALSE])
  df$outcome <- y_factor_by_id[model_ids]
  df <- df[complete.cases(df), , drop = FALSE]

  Xs <- prep_glmnet(dplyr::select(df, -outcome))
  y_local <- as.integer(df$outcome == "HER2pos")

  cv_fit <- cv.glmnet(Xs, y_local, family = "binomial", alpha = 0.5,
                      nfolds = cv_folds, type.measure = "auc",
                      standardize = FALSE)
  coefs <- coef(cv_fit, s = "lambda.1se") %>%
    as.matrix() %>% as.data.frame() %>%
    rownames_to_column("Feature") %>%
    rename(Coefficient = s1) %>%
    filter(Feature != "(Intercept)", Coefficient != 0) %>%
    mutate(Model = label,
           Direction = ifelse(Coefficient > 0, "HER2+ associated", "HER2– associated"))
  list(fit = cv_fit, coefs = coefs)
}

lr_combined <- fit_lr_full(X_combined, "M4: Amplicon + CNV + genomic (glmnet)")

p_coef <- lr_combined$coefs %>%
  ggplot(aes(x = reorder(Feature, Coefficient),
             y = Coefficient, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = c("HER2+ associated" = "#E64B35",
                                "HER2– associated" = "#4DBBD5")) +
  labs(title    = "Logistic Regression Coefficients (glmnet, lambda.1se)",
       subtitle  = "Standardised features — coefficients are directly comparable",
       x = NULL, y = "Coefficient (log-odds per SD)", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
save_plot(p_coef, "04_logistic_coefficients.png", w = 8, h = 5)

# ═══════════════════════════════════════════════════════════
# 6. SHAP VALUES (Random Forest — M5)
# ═══════════════════════════════════════════════════════════
# kernelshap works with any predict() method; shapviz provides plots

rf_df <- as.data.frame(X_rf_all)
rf_model_ids <- intersect(rownames(rf_df), names(y_factor_by_id))
rf_df <- rf_df[rf_model_ids, , drop = FALSE]
rf_df$outcome <- y_factor_by_id[rf_model_ids]
rf_df <- rf_df[complete.cases(rf_df), , drop = FALSE]

X_rf <- dplyr::select(rf_df, -outcome)
y_rf <- rf_df$outcome

# Refit RF on full data for SHAP
set.seed(random_seed)
rf_full <- randomForest(
  x = scale(X_rf),
  y = y_rf,
  ntree = 500,
  importance = TRUE
)

# SHAP via kernelshap (model-agnostic)
bg_idx <- sample(nrow(X_rf), min(100, nrow(X_rf)))
X_rf_scaled <- as.data.frame(scale(X_rf))

ks <- kernelshap(
  rf_full,
  X    = X_rf_scaled,
  bg_X = X_rf_scaled[bg_idx, ],
  pred_fun = function(m, X) predict(m, X, type = "prob")[, "HER2pos"]
)
sv <- shapviz(ks)

# SHAP beeswarm (summary plot)
png(file.path(out_dir, "05_SHAP_beeswarm.png"),
    width = 900, height = 700, res = 130)
sv_importance(sv, kind = "beeswarm", max_display = 15,
              viridis_args = list(option = "C")) +
  ggtitle("SHAP Values — M5 Random Forest (HER2+ probability)")
dev.off()

# SHAP bar (mean |SHAP|)
png(file.path(out_dir, "06_SHAP_importance_bar.png"),
    width = 800, height = 600, res = 130)
sv_importance(sv, kind = "bar", max_display = 15,
              fill = "#3C5488") +
  ggtitle("Mean |SHAP| — Feature Importance (RF)")
dev.off()

# SHAP dependence plot for top 2 features
top_feats <- sv_importance(sv, kind = "bar")$data %>%
  arrange(desc(value)) %>% pull(feature) %>% head(2)

for (feat in top_feats) {
  p_dep <- sv_dependence(sv, v = feat, color_var = "auto") +
    theme_bw(base_size = 11) +
    labs(title = paste0("SHAP Dependence: ", feat))
  save_plot(p_dep,
            paste0("07_SHAP_dependence_", make.names(feat), ".png"),
            w = 6, h = 5)
}

# RF permutation importance as cross-check
imp_df <- importance(rf_full, type = 1) %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  rename(MDA = MeanDecreaseAccuracy) %>%
  arrange(desc(MDA))

p_imp <- ggplot(imp_df, aes(x = reorder(Feature, MDA), y = MDA)) +
  geom_col(fill = "#7E6148") +
  coord_flip() +
  labs(title = "RF Permutation Importance (Mean Decrease Accuracy)",
       x = NULL, y = "Mean Decrease Accuracy") +
  theme_bw(base_size = 11)
save_plot(p_imp, "08_RF_permutation_importance.png", w = 7, h = 5)

# ═══════════════════════════════════════════════════════════
# 7. CONCORDANCE / DISCORDANCE PLOTS
# ═══════════════════════════════════════════════════════════

# ── 7A. RNA vs CNV scatter coloured by IHC label ──────────
conc_df <- data.frame(
  ERBB2_mRNA   = X_rna_erbb2$ERBB2_mRNA,
  ERBB2_GISTIC = X_cnv$ERBB2_score,
  IHC_label    = meta$her2_clinical,
  row.names    = rownames(meta)
) %>%
  mutate(
    concordant = case_when(
      IHC_label == "HER2pos" & ERBB2_GISTIC >= 1 & ERBB2_mRNA > median(ERBB2_mRNA) ~ "True Positive",
      IHC_label == "HER2neg" & ERBB2_GISTIC < 1  & ERBB2_mRNA <= median(ERBB2_mRNA) ~ "True Negative",
      IHC_label == "HER2pos" & ERBB2_GISTIC < 1  ~ "IHC+, CNV neutral (RNA-driven?)",
      IHC_label == "HER2pos" & ERBB2_mRNA <= median(ERBB2_mRNA) ~ "IHC+, low mRNA",
      IHC_label == "HER2neg" & ERBB2_GISTIC >= 1 ~ "IHC−, CNV amplified",
      TRUE ~ "Other"
    )
  )

conc_pal <- c(
  "True Positive"                   = "#E64B35",
  "True Negative"                   = "#4DBBD5",
  "IHC+, CNV neutral (RNA-driven?)" = "#F39B7F",
  "IHC+, low mRNA"                  = "#984EA3",
  "IHC−, CNV amplified"             = "#FF7F00",
  "Other"                           = "grey70"
)

p_conc <- ggplot(conc_df,
                 aes(ERBB2_GISTIC, ERBB2_mRNA,
                     color = concordant, shape = IHC_label)) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.75) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = median(conc_df$ERBB2_mRNA),
             linetype = "dashed", color = "grey40") +
  scale_color_manual(values = conc_pal) +
  scale_shape_manual(values = c("HER2pos" = 17, "HER2neg" = 16)) +
  labs(title    = "ERBB2: RNA vs CNV Concordance",
       subtitle  = "Coloured by concordance with IHC clinical label",
       x = "ERBB2 GISTIC2 Score", y = "ERBB2 VST Expression",
       color = "Concordance class", shape = "IHC label") +
  theme_bw(base_size = 11) +
  theme(legend.text = element_text(size = 8))
save_plot(p_conc, "09_RNA_CNV_concordance_scatter.png", w = 9, h = 6)

# ── 7B. UpSet-style concordance counts ────────────────────
conc_df$RNA_high <- as.integer(conc_df$ERBB2_mRNA > median(conc_df$ERBB2_mRNA))
conc_df$CNV_amp  <- as.integer(conc_df$ERBB2_GISTIC >= 1)
conc_df$IHC_pos  <- as.integer(conc_df$IHC_label == "HER2pos")

upset_df <- conc_df %>%
  count(RNA_high, CNV_amp, IHC_pos) %>%
  mutate(
    combo = paste0(
      ifelse(RNA_high, "RNA↑ ", "RNA↓ "),
      ifelse(CNV_amp,  "CNV+ ", "CNV- "),
      ifelse(IHC_pos,  "IHC+",  "IHC-")
    )
  )

p_upset <- ggplot(upset_df, aes(x = reorder(combo, n), y = n,
                                 fill = factor(IHC_pos))) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("1" = "#E64B35", "0" = "#4DBBD5"),
                    labels = c("1" = "IHC HER2+", "0" = "IHC HER2-")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "RNA / CNV / IHC Concordance Counts",
       x = NULL, y = "Number of samples", fill = "IHC label") +
  theme_bw(base_size = 11)
save_plot(p_upset, "10_concordance_counts.png", w = 8, h = 5)

# ── 7C. Distribution of ERBB2 mRNA by GISTIC tier & IHC ──
conc_df$GISTIC_tier <- factor(
  case_when(
    conc_df$ERBB2_GISTIC == -2 ~ "Deep del (-2)",
    conc_df$ERBB2_GISTIC == -1 ~ "Shallow del (-1)",
    conc_df$ERBB2_GISTIC ==  0 ~ "Neutral (0)",
    conc_df$ERBB2_GISTIC ==  1 ~ "Gain (1)",
    conc_df$ERBB2_GISTIC ==  2 ~ "Amplified (2)"
  ),
  levels = c("Deep del (-2)", "Shallow del (-1)", "Neutral (0)", "Gain (1)", "Amplified (2)")
)

p_violin <- ggplot(conc_df,
                   aes(GISTIC_tier, ERBB2_mRNA, fill = IHC_label)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = palette) +
  labs(title    = "ERBB2 mRNA by GISTIC2 Copy-Number Tier",
       subtitle  = "Filled by IHC clinical label",
       x = "GISTIC2 Score", y = "ERBB2 VST Expression",
       fill = "IHC label") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_violin, "11_mRNA_by_GISTIC_tier.png", w = 9, h = 5)

# ── 7D. Confusion matrices — RNA-only vs CNV-only vs combined
plot_cm <- function(model, title) {
  preds <- model$pred %>%
    group_by(rowIndex) %>%
    summarise(pred_class = names(which.max(table(pred))),
              obs = first(obs), .groups = "drop")
  cm <- confusionMatrix(factor(preds$pred_class, levels = c("HER2neg", "HER2pos")),
                        preds$obs)
  cm_df <- as.data.frame(cm$table) %>%
    rename(Predicted = Prediction, Actual = Reference)

  ggplot(cm_df, aes(Actual, Predicted, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#3C5488") +
    labs(title = title,
         subtitle = paste0("Sens=", round(cm$byClass["Sensitivity"], 3),
                           "  Spec=", round(cm$byClass["Specificity"], 3),
                           "  PPV=",  round(cm$byClass["Pos Pred Value"], 3))) +
    theme_bw(base_size = 12) + theme(legend.position = "none")
}

p_cm_rna  <- plot_cm(m1_rna,      "Confusion Matrix — ERBB2 mRNA alone")
p_cm_cnv  <- plot_cm(m2_cnv,      "Confusion Matrix — ERBB2 CNV alone")
p_cm_comb <- plot_cm(m4_combined, "Confusion Matrix — Combined (M4 glmnet)")

save_plot(p_cm_rna | p_cm_cnv | p_cm_comb,
          "12_confusion_matrices.png", w = 14, h = 5)

# ── 7E. FGA and mutation count by HER2 label ─────────────
p_fga <- ggplot(meta, aes(her2_clinical, Fraction_Genome_Altered,
                           fill = her2_clinical)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = palette) +
  labs(title = "Fraction Genome Altered by HER2 Status",
       x = NULL, y = "FGA") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_mut <- ggplot(meta, aes(her2_clinical, log10(Mutation_Count + 1),
                           fill = her2_clinical)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = palette) +
  labs(title = "Mutation Count (log10) by HER2 Status",
       x = NULL, y = "log10(mutation count + 1)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

save_plot(p_fga + p_mut, "13_FGA_mutation_by_HER2.png", w = 9, h = 5)

# ═══════════════════════════════════════════════════════════
# 8. SUMMARY TABLE
# ═══════════════════════════════════════════════════════════
cat("\n══════════════════════════════════════════\n")
cat("  Modality Comparison — Cross-validated AUC\n")
cat("══════════════════════════════════════════\n")
print(cv_metrics, digits = 3, row.names = FALSE)

write.csv(conc_df,           file.path(res_dir, "sample_concordance_table.csv"))
write.csv(imp_df,            file.path(res_dir, "RF_feature_importance.csv"),    row.names = FALSE)
write.csv(lr_combined$coefs, file.path(res_dir, "logistic_coefficients.csv"),   row.names = FALSE)

cat("\n✓ Plots saved to:", out_dir, "\n")
cat("✓ Tables saved to:", res_dir, "\n")