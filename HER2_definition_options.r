# ============================================================
# Multimodal HER2 Definition: RNA vs DNA vs Combined
# Modalities: ERBB2 mRNA (VST) | GISTIC2 CNV | FGA | Mutation count
# Outputs: ROC/AUC, logistic coefficients, SHAP, concordance plots
# ============================================================
# Expected inputs (from previous QC + DESeq2 scripts):
#   metadata_qc_passed.csv  – must contain: her2_clinical (HER2pos/HER2neg),
#                             fraction_genome_altered, mutation_count
#   vst_her2_brca.rds       – VST SummarizedExperiment from DESeq2 script
#   erbb2_gistic2.csv       – GISTIC2 scores, rows=samples, col="ERBB2"
#                             values: -2 (deep del) / -1 / 0 / 1 / 2 (amp)
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

OUT_DIR <- "multimodal_plots";   dir.create(OUT_DIR, showWarnings = FALSE)
RES_DIR <- "multimodal_results"; dir.create(RES_DIR, showWarnings = FALSE)

CV_FOLDS     <- 5
RANDOM_SEED  <- 42
ERBB2_ENS    <- "ENSG00000141736"   # adjust if IDs carry version suffix e.g. .13
PALETTE      <- c("HER2pos" = "#E64B35", "HER2neg" = "#4DBBD5")

save_plot <- function(p, fname, w = 8, h = 6, dpi = 150)
  ggsave(file.path(OUT_DIR, fname), p, width = w, height = h, dpi = dpi, bg = "white")

# ═══════════════════════════════════════════════════════════
# 1. LOAD & ALIGN DATA
# ═══════════════════════════════════════════════════════════

meta <- read.csv("metadata_qc_passed.csv", row.names = 1)
vst  <- readRDS("vst_her2_brca.rds")            # RangedSummarizedExperiment
cnv  <- read.csv("erbb2_gistic2.csv", row.names = 1)   # samples × genes

# ── Harmonise TCGA barcodes to 15 chars ───────────────────
trim_bc <- function(x, n = 15) substr(x, 1, n)

rownames(meta)   <- trim_bc(rownames(meta))
colnames(vst)    <- trim_bc(colnames(vst))
rownames(cnv)    <- trim_bc(rownames(cnv))

# Keep only samples with definitive IHC label
meta <- meta[meta$her2_clinical %in% c("HER2pos", "HER2neg"), ]
meta$her2_binary <- as.integer(meta$her2_clinical == "HER2pos")

# Intersect across all layers
common_samples <- Reduce(intersect, list(
  rownames(meta),
  colnames(vst),
  rownames(cnv)
))

meta <- meta[common_samples, ]
cnv  <- cnv[common_samples, , drop = FALSE]
vst  <- vst[, common_samples]

cat(sprintf("Common samples: %d  |  HER2+: %d  |  HER2-: %d\n",
            length(common_samples),
            sum(meta$her2_binary), sum(1 - meta$her2_binary)))

# ═══════════════════════════════════════════════════════════
# 2. BUILD FEATURE MATRICES
# ═══════════════════════════════════════════════════════════

vst_mat <- assay(vst)   # genes × samples

# ── 2A. RNA: ERBB2 mRNA ───────────────────────────────────
erbb2_row <- rownames(vst_mat)[str_detect(rownames(vst_mat), ERBB2_ENS)]
if (length(erbb2_row) == 0) stop("ERBB2 Ensembl ID not found in VST matrix")

X_rna_erbb2 <- data.frame(
  ERBB2_mRNA = as.numeric(vst_mat[erbb2_row[1], ]),
  row.names  = colnames(vst_mat)
)

# ── 2B. RNA: Expanded — ERBB2 amplicon neighbours ─────────
# GRB7, STARD3, MIEN1 are co-amplified on chr17q12; biologically meaningful
amplicon_genes <- c("ENSG00000141736",  # ERBB2
                    "ENSG00000177455",  # GRB7
                    "ENSG00000132481",  # STARD3
                    "ENSG00000132470")  # MIEN1

amp_rows <- rownames(vst_mat)[
  str_detect(rownames(vst_mat),
             paste(amplicon_genes, collapse = "|"))
]

X_rna_amplicon <- t(vst_mat[amp_rows, ]) %>%
  as.data.frame() %>%
  setNames(str_remove(colnames(.), "\\.\\d+$"))

# ── 2C. DNA: ERBB2 GISTIC2 score ─────────────────────────
# GISTIC2 score: -2=deep deletion, -1=shallow del, 0=neutral, 1=gain, 2=amp
cnv_col <- intersect(c("ERBB2", "ERBB2_gistic"), colnames(cnv))[1]
if (is.na(cnv_col)) stop("Cannot find ERBB2 column in CNV file")

X_cnv <- data.frame(
  ERBB2_GISTIC = cnv[[cnv_col]],
  row.names    = rownames(cnv)
)

# Binary amplification flag (GISTIC >= 1 = gain/amplification)
X_cnv$ERBB2_amp_flag <- as.integer(X_cnv$ERBB2_GISTIC >= 1)

# ── 2D. Genomic instability features ─────────────────────
X_genomic <- meta %>%
  dplyr::select(fraction_genome_altered, mutation_count) %>%
  mutate(
    fga_scaled  = scale(fraction_genome_altered)[, 1],
    mut_log     = log10(mutation_count + 1)
  ) %>%
  dplyr::select(fga_scaled, mut_log)

# ── 2E. Combined feature matrix ───────────────────────────
X_combined <- bind_cols(
  X_rna_erbb2,
  X_cnv,
  X_genomic
) %>% as.data.frame()

# Outcome vector (aligned)
y        <- meta$her2_binary
y_factor <- factor(ifelse(y == 1, "HER2pos", "HER2neg"),
                   levels = c("HER2neg", "HER2pos"))

# ═══════════════════════════════════════════════════════════
# 3. CROSS-VALIDATED MODEL EVALUATION
# ═══════════════════════════════════════════════════════════
# We compare 5 models:
#   M1: ERBB2 mRNA alone
#   M2: ERBB2 GISTIC2 score alone
#   M3: ERBB2 GISTIC2 + amplification flag
#   M4: ERBB2 mRNA + CNV + genomic instability
#   M5: ERBB2 amplicon genes + CNV + genomic instability (RF)

set.seed(RANDOM_SEED)
cv_ctrl <- trainControl(
  method          = "cv",
  number          = CV_FOLDS,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling        = "up"          # upsample minority class (HER2+ ~20%)
)

fit_model <- function(X, label, method = "glm") {
  df <- bind_cols(X, outcome = y_factor)
  train(outcome ~ .,
        data      = df,
        method    = method,
        trControl = cv_ctrl,
        metric    = "ROC",
        preProcess = c("center", "scale"))
}

cat("\nFitting models...\n")
m1_rna     <- fit_model(X_rna_erbb2,           "M1: ERBB2 mRNA")
m2_cnv     <- fit_model(X_cnv["ERBB2_GISTIC"], "M2: ERBB2 CNV")
m3_cnv_bin <- fit_model(X_cnv,                 "M3: ERBB2 CNV + amp flag")
m4_combined <- fit_model(X_combined,            "M4: mRNA + CNV + genomic",
                         method = "glmnet")
m5_rf       <- fit_model(
  bind_cols(X_rna_amplicon, X_cnv, X_genomic),
  "M5: Amplicon + CNV + genomic (RF)",
  method = "rf"
)

models <- list(
  "M1: ERBB2 mRNA"              = m1_rna,
  "M2: ERBB2 CNV (GISTIC2)"     = m2_cnv,
  "M3: ERBB2 CNV + amp flag"    = m3_cnv_bin,
  "M4: mRNA + CNV + genomic"    = m4_combined,
  "M5: Amplicon + CNV + RF"     = m5_rf
)

# Extract CV metrics
cv_metrics <- lapply(names(models), function(nm) {
  preds <- models[[nm]]$pred
  # one row per sample per fold; average across resamples
  preds %>%
    group_by(rowIndex) %>%
    summarise(prob = mean(HER2pos), obs = first(obs), .groups = "drop") %>%
    summarise(
      Model    = nm,
      ROC_AUC  = roc_auc_score_r(obs, prob),    # see helper below
      PR_AUC   = pr_auc_score_r(obs, prob)
    )
}) %>% bind_rows()

# R helper wrappers (pROC + yardstick)
roc_auc_score_r <- function(obs, prob) {
  roc(as.integer(obs == "HER2pos"), prob, quiet = TRUE)$auc %>% as.numeric()
}
pr_auc_score_r <- function(obs, prob) {
  df <- data.frame(truth = factor(obs, levels = c("HER2neg","HER2pos")),
                   prob  = prob)
  yardstick::pr_auc(df, truth, prob, event_level = "second")$.estimate
}

# Recompute now helpers are defined
cv_metrics <- lapply(names(models), function(nm) {
  preds <- models[[nm]]$pred %>%
    group_by(rowIndex) %>%
    summarise(prob = mean(HER2pos), obs = first(obs), .groups = "drop")
  data.frame(
    Model   = nm,
    ROC_AUC = roc_auc_score_r(preds$obs, preds$prob),
    PR_AUC  = pr_auc_score_r(preds$obs, preds$prob)
  )
}) %>% bind_rows()

print(cv_metrics)
write.csv(cv_metrics, file.path(RES_DIR, "model_cv_metrics.csv"), row.names = FALSE)

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
       subtitle  = paste0(CV_FOLDS, "-fold cross-validated"),
       x = "False Positive Rate", y = "True Positive Rate",
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.62, 0.22),
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.text = element_text(size = 8))
save_plot(p_roc, "01_ROC_curves.png", w = 8, h = 6)

# ── Precision-Recall curves ───────────────────────────────
pr_list <- lapply(names(models), function(nm) {
  preds <- models[[nm]]$pred %>%
    group_by(rowIndex) %>%
    summarise(prob = mean(HER2pos), obs = first(obs), .groups = "drop")
  pr <- pr.curve(scores.class0 = preds$prob[preds$obs == "HER2pos"],
                 scores.class1 = preds$prob[preds$obs == "HER2neg"],
                 curve = TRUE)
  data.frame(Recall = pr$curve[,1], Precision = pr$curve[,2],
             Model  = nm, PRAUC = round(pr$auc.integral, 3))
}) %>% bind_rows() %>%
  mutate(Model_AUC = paste0(Model, "  (PRAUC=", PRAUC, ")"))

# PRROC needed for pr.curve — install if absent
if (!requireNamespace("PRROC", quietly=TRUE)) install.packages("PRROC")
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

baseline_prev <- mean(y)

p_pr <- ggplot(pr_list, aes(Recall, Precision, color = Model_AUC)) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = baseline_prev, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.8, y = baseline_prev + 0.02,
           label = paste0("Baseline (prevalence=", round(baseline_prev,2), ")"),
           size = 3, color = "grey40") +
  scale_color_brewer(palette = "Set1") +
  labs(title    = "Precision-Recall Curves — HER2+ Prediction",
       subtitle  = paste0(CV_FOLDS, "-fold cross-validated"),
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
# Refit a single glmnet on full data for coefficient interpretation
# (CV already confirmed model quality; coefficients on full data
#  are the most stable estimates)

prep_glmnet <- function(X) {
  scale(as.matrix(X))   # z-score so coefficients are comparable
}

fit_lr_full <- function(X, label) {
  Xs <- prep_glmnet(X)
  cv_fit <- cv.glmnet(Xs, y, family = "binomial", alpha = 0.5,
                      nfolds = CV_FOLDS, type.measure = "auc",
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

lr_combined <- fit_lr_full(X_combined, "M4: mRNA + CNV + genomic")

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

X_rf <- bind_cols(X_rna_amplicon, X_cnv, X_genomic) %>% as.data.frame()

# Refit RF on full data for SHAP
set.seed(RANDOM_SEED)
rf_full <- randomForest(
  x = scale(X_rf),
  y = y_factor,
  ntree = 500,
  importance = TRUE
)

# SHAP via kernelshap (model-agnostic)
# Use a background sample of 100 for speed
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
png(file.path(OUT_DIR, "05_SHAP_beeswarm.png"),
    width = 900, height = 700, res = 130)
sv_importance(sv, kind = "beeswarm", max_display = 15,
              viridis_args = list(option = "C")) +
  ggtitle("SHAP Values — M5 Random Forest (HER2+ probability)")
dev.off()

# SHAP bar (mean |SHAP|)
png(file.path(OUT_DIR, "06_SHAP_importance_bar.png"),
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
imp_df <- importance(rf_full, type = 1) %>%   # type=1: mean decrease accuracy
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
  ERBB2_mRNA  = X_rna_erbb2$ERBB2_mRNA,
  ERBB2_GISTIC = X_cnv$ERBB2_GISTIC,
  IHC_label   = meta$her2_clinical,
  row.names   = rownames(meta)
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
  "True Positive"                  = "#E64B35",
  "True Negative"                  = "#4DBBD5",
  "IHC+, CNV neutral (RNA-driven?)"= "#F39B7F",
  "IHC+, low mRNA"                 = "#984EA3",
  "IHC−, CNV amplified"            = "#FF7F00",
  "Other"                          = "grey70"
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
# For each sample: is RNA high? CNV amplified? IHC+?
conc_df$RNA_high   <- as.integer(conc_df$ERBB2_mRNA > median(conc_df$ERBB2_mRNA))
conc_df$CNV_amp    <- as.integer(conc_df$ERBB2_GISTIC >= 1)
conc_df$IHC_pos    <- as.integer(conc_df$IHC_label == "HER2pos")

upset_df <- conc_df %>%
  count(RNA_high, CNV_amp, IHC_pos) %>%
  mutate(
    combo = paste0(
      ifelse(RNA_high, "RNA↑ ", "RNA↓ "),
      ifelse(CNV_amp,  "CNV+ ",  "CNV- "),
      ifelse(IHC_pos,  "IHC+",   "IHC-")
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
  levels = c("Deep del (-2)","Shallow del (-1)","Neutral (0)","Gain (1)","Amplified (2)")
)

p_violin <- ggplot(conc_df,
                   aes(GISTIC_tier, ERBB2_mRNA, fill = IHC_label)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5,
               position = position_dodge(0.9)) +
  scale_fill_manual(values = PALETTE) +
  labs(title    = "ERBB2 mRNA by GISTIC2 Copy-Number Tier",
       subtitle  = "Filled by IHC clinical label",
       x = "GISTIC2 Score", y = "ERBB2 VST Expression",
       fill = "IHC label") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_violin, "11_mRNA_by_GISTIC_tier.png", w = 9, h = 5)

# ── 7D. Confusion matrices — best RNA-only vs best DNA-only
plot_cm <- function(model, title, fname) {
  preds <- model$pred %>%
    group_by(rowIndex) %>%
    summarise(pred_class = names(which.max(table(pred))),
              obs = first(obs), .groups = "drop")
  cm <- confusionMatrix(factor(preds$pred_class, levels = c("HER2neg","HER2pos")),
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

p_cm_rna <- plot_cm(m1_rna,  "Confusion Matrix — ERBB2 mRNA alone",  "12a")
p_cm_cnv <- plot_cm(m2_cnv,  "Confusion Matrix — ERBB2 CNV alone",   "12b")
p_cm_comb <- plot_cm(m4_combined, "Confusion Matrix — Combined (M4)", "12c")

save_plot(p_cm_rna | p_cm_cnv | p_cm_comb,
          "12_confusion_matrices.png", w = 14, h = 5)

# ── 7E. FGA and mutation count by HER2 label ─────────────
p_fga <- ggplot(meta, aes(her2_clinical, fraction_genome_altered,
                           fill = her2_clinical)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = PALETTE) +
  labs(title = "Fraction Genome Altered by HER2 Status",
       x = NULL, y = "FGA") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_mut <- ggplot(meta, aes(her2_clinical, log10(mutation_count + 1),
                           fill = her2_clinical)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = PALETTE) +
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

write.csv(conc_df,    file.path(RES_DIR, "sample_concordance_table.csv"))
write.csv(imp_df,     file.path(RES_DIR, "RF_feature_importance.csv"), row.names = FALSE)
write.csv(lr_combined$coefs, file.path(RES_DIR, "logistic_coefficients.csv"), row.names = FALSE)

cat("\n✓ Plots saved to:", OUT_DIR, "\n")
cat("✓ Tables saved to:", RES_DIR, "\n")