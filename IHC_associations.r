# ============================================================
# RNA vs DNA: which is more predictive of IHC HER2 label?
# Two univariate logistic regressions, head-to-head ROC comparison
# ============================================================
# Requires from previous scripts:
#   meta  – with columns:
#             her2_ihc_status  ("Positive" / "Negative" / other)
#             erbb2_copy_number (-2 to 2)
#   vst   – VST SummarizedExperiment
# ============================================================

library(DESeq2); library(SummarizedExperiment)
library(dplyr); library(tibble); library(ggplot2); 
library(pROC); library(patchwork)

out_dir <- "rna_vs_dna"; dir.create(out_dir, showWarnings = FALSE)

# ── 1. OUTCOME ─────────────────────────────────────────────
# Keep only unambiguous Positive / Negative IHC calls
df_outcome <- meta %>%
  filter(IHC_HER2 %in% c("Positive", "Negative")) %>%
  mutate(y = as.integer(IHC_HER2 == "Positive"))

cat(sprintf("Samples after filtering: %d  |  IHC+: %d  |  IHC-: %d\n",
            nrow(df_outcome), sum(df_outcome$y), sum(df_outcome$y == 0)))

# ── 2. FEATURES ────────────────────────────────────────────
vst_mat <- assay(vst)

erbb2_row <- grep("^ERBB2$", rownames(vst_mat), value = TRUE)[1]
if (is.na(erbb2_row)) stop("ERBB2 not found in VST matrix rownames")

sample_ids <- rownames(df_outcome)
df_outcome$rna <- as.numeric(vst_mat[erbb2_row, sample_ids])
df_outcome$dna <- df_outcome$erbb2_copy_number   

# Drop any samples missing either feature
df_model <- df_outcome %>%
  filter(!is.na(rna), !is.na(dna)) %>%
  mutate(rna_z = scale(rna)[, 1],     # standardise
         dna_z = scale(dna)[, 1])

cat(sprintf("Complete cases for modelling: %d\n", nrow(df_model)))

# ── 3. LOGISTIC REGRESSION ─────────────────────────────────
fit_rna <- glm(y ~ rna_z, data = df_model, family = binomial)
fit_dna <- glm(y ~ dna_z, data = df_model, family = binomial)
fit_both <- glm(y ~ rna_z + dna_z, data = df_model, family = binomial)

# ── 4. ROC CURVES ──────────────────────────────────────────
roc_rna  <- roc(df_model$y, predict(fit_rna,  type = "response"), quiet = TRUE)
roc_dna  <- roc(df_model$y, predict(fit_dna,  type = "response"), quiet = TRUE)
roc_both <- roc(df_model$y, predict(fit_both, type = "response"), quiet = TRUE)

auc_rna  <- round(as.numeric(auc(roc_rna)),  3)
auc_dna  <- round(as.numeric(auc(roc_dna)),  3)
auc_both <- round(as.numeric(auc(roc_both)), 3)

# DeLong test: is the difference between RNA and DNA AUC significant?
delong <- roc.test(roc_rna, roc_dna, method = "delong")

cat(sprintf("\nAUC — RNA: %.3f  |  DNA: %.3f  |  Combined: %.3f\n",
            auc_rna, auc_dna, auc_both))
cat(sprintf("DeLong test (RNA vs DNA): D = %.3f, p = %.4f\n",
            delong$statistic, delong$p.value))

# ── 5. PLOT: ROC CURVES ────────────────────────────────────
roc_df <- bind_rows(
  data.frame(FPR = 1 - roc_rna$specificities,  TPR = roc_rna$sensitivities,
             Model = paste0("RNA (AUC=", auc_rna, ")")),
  data.frame(FPR = 1 - roc_dna$specificities,  TPR = roc_dna$sensitivities,
             Model = paste0("DNA (AUC=", auc_dna, ")")),
  data.frame(FPR = 1 - roc_both$specificities, TPR = roc_both$sensitivities,
             Model = paste0("Combined (AUC=", auc_both, ")"))
)

p_roc <- ggplot(roc_df, aes(FPR, TPR, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed", color = "grey60") +
  annotate("text", x = 0.6, y = 0.15,
           label = sprintf("DeLong p = %.4f", delong$p.value),
           size = 3.5, color = "grey30") +
  scale_color_manual(values = c("#E64B35", "#4DBBD5", "#7E6148")) +
  labs(title   = "ROC Curves: RNA vs DNA predicting IHC HER2+",
       x = "False Positive Rate", y = "True Positive Rate", color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.65, 0.2),
        legend.background = element_rect(fill = "white", color = "grey80"))

# ── 6. PLOT: FEATURE DISTRIBUTIONS BY IHC LABEL ───────────
df_long <- df_model %>%
  dplyr::select(y, IHC_HER2, rna, dna) %>%
  tidyr::pivot_longer(c(rna, dna), names_to = "Modality", values_to = "Value") %>%
  mutate(Modality = ifelse(Modality == "rna", "ERBB2 mRNA (VST)", "ERBB2 CNV"))

p_dist <- ggplot(df_long, aes(IHC_HER2, Value, fill = IHC_HER2)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = c("Positive" = "#E64B35", "Negative" = "#4DBBD5")) +
  facet_wrap(~Modality, scales = "free_y") +
  labs(title = "ERBB2 RNA and DNA by IHC Label",
       x = NULL, y = "Feature Value") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# ── 7. PLOT: COEFFICIENTS ──────────────────────────────────
coef_df <- bind_rows(
  broom::tidy(fit_rna,  conf.int = TRUE) %>% mutate(Model = "RNA only"),
  broom::tidy(fit_dna,  conf.int = TRUE) %>% mutate(Model = "DNA only"),
  broom::tidy(fit_both, conf.int = TRUE) %>% mutate(Model = "Combined")
) %>% filter(term != "(Intercept)")

p_coef <- ggplot(coef_df, aes(x = Model, y = estimate,
                                group = term, fill = term,
                                ymin = conf.low, ymax = conf.high,
                                color = term)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(position = position_dodge(0.5), size = 0.6) +
  scale_color_manual(values = c("rna_z" = "#E64B35",
                                 "dna_z" = "#4DBBD5",
                                 "both_z" = "#7E6148")) +
  labs(title   = "Logistic Regression Coefficients (standardised)",
       x = NULL, y = "Log-odds (95% CI)", color = NULL) +
  theme_bw(base_size = 12)

# ── 8. SAVE ────────────────────────────────────────────────
ggsave(file.path(out_dir, "01_ROC_RNA_vs_DNA.png"),        p_roc,  width=7, h=6, dpi=300, bg="white")
ggsave(file.path(out_dir, "02_distributions_by_IHC.png"),  p_dist, width=8, h=5, dpi=300, bg="white")
ggsave(file.path(out_dir, "03_logistic_coefficients.png"), p_coef, width=7, h=4, dpi=300, bg="white")

results <- data.frame(
  Model    = c("RNA only", "DNA only", "Combined"),
  AUC      = c(auc_rna, auc_dna, auc_both),
  Coef     = c(coef(fit_rna)["rna_z"], coef(fit_dna)["dna_z"], NA),
  DeLong_p = c(delong$p.value, delong$p.value, NA)
)
write.csv(results, file.path(out_dir, "RNA_vs_DNA_results.csv"), row.names = FALSE)

cat("\n✓ Done. Results in:", out_dir, "\n")


#exploring
# ============================================================
# Diagnose the ROC jump: who are the "obvious" HER2+ cases?
# ============================================================

# ── 1. Find the jump threshold ─────────────────────────────
# The jump corresponds to the predicted probability at which
# TPR rises steeply before FPR has moved much.
# Operationally: find the threshold where TPR > 0.3 and FPR < 0.05

roc_df_rna <- data.frame(
  threshold = roc_rna$thresholds,
  TPR       = roc_rna$sensitivities,
  FPR       = 1 - roc_rna$specificities
)

jump_threshold <- roc_df_rna %>%
  filter(FPR < 0.05, TPR > 0.3) %>%
  slice_max(TPR, n = 1, with_ties = FALSE) %>%
  pull(threshold)

cat(sprintf("Jump threshold (RNA model): %.3f predicted probability\n", jump_threshold))

# ── 2. Label each HER2+ sample as "obvious" or "hard" ──────
df_model$prob_rna <- predict(fit_rna, type = "response")
df_model$prob_dna <- predict(fit_dna, type = "response")

df_model <- df_model %>%
  mutate(
    case_type = case_when(
      y == 0                         ~ "HER2-",
      prob_rna >= jump_threshold     ~ "Obvious HER2+",
      TRUE                           ~ "Hard HER2+"
    )
  )

cat("\nCase breakdown:\n")
print(table(df_model$case_type))

# ── 3. What distinguishes obvious from hard cases? ─────────
p_rna <- ggplot(df_model %>% filter(y == 1),
                aes(x = case_type, y = rna, fill = case_type)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = c("Obvious HER2+" = "#E64B35", "Hard HER2+" = "#F39B7F")) +
  labs(title = "ERBB2 mRNA", x = NULL, y = "VST expression") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

p_dna <- ggplot(df_model %>% filter(y == 1),
                aes(x = case_type, y = dna, fill = case_type)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.12, outlier.size = 0.5) +
  scale_fill_manual(values = c("Obvious HER2+" = "#4DBBD5", "Hard HER2+" = "#91D1C2")) +
  labs(title = "ERBB2 CNV (GISTIC2)", x = NULL, y = "GISTIC2 score") +
  theme_bw(base_size = 11) + theme(legend.position = "none")

# Scatter: RNA vs DNA, coloured by case type (all samples)
p_scatter <- ggplot(df_model, aes(dna, rna, color = case_type)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("HER2-"         = "grey70",
                                 "Obvious HER2+" = "#E64B35",
                                 "Hard HER2+"    = "#F39B7F")) +
  labs(title   = "RNA vs DNA — obvious vs hard HER2+ cases",
       x = "ERBB2 CNV (GISTIC2)", y = "ERBB2 mRNA (VST)",
       color = NULL) +
  theme_bw(base_size = 11)

ggsave(file.path(out_dir, "04_obvious_vs_hard_features.png"),
       p_rna + p_dna, width = 8, h = 5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "05_RNA_vs_DNA_scatter_casetype.png"),
       p_scatter, width = 7, h = 5, dpi = 300, bg = "white")

# ── 4. Do RNA and DNA agree on which cases are obvious? ────
# If both modalities flag the same patients as obvious,
# they're redundant on easy cases. If they disagree, one
# modality is capturing a distinct subgroup.
df_model <- df_model %>%
  mutate(
    obvious_rna = as.integer(y == 1 & prob_rna >= jump_threshold),
    obvious_dna = as.integer(y == 1 & prob_dna >= jump_threshold)
  )

agreement <- df_model %>%
  filter(y == 1) %>%
  count(obvious_rna, obvious_dna) %>%
  mutate(label = paste0(
    ifelse(obvious_rna, "RNA obvious", "RNA hard"), " / ",
    ifelse(obvious_dna, "DNA obvious", "DNA hard")
  ))

cat("\nRNA vs DNA agreement on 'obvious' cases:\n")
print(agreement)
