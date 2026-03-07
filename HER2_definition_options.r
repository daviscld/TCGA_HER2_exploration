# ============================================================
# Multimodal HER2 Definition
# Layers: IHC | FISH | ERBB2 CNV | ERBB2 mRNA | Amplicon genes
# ============================================================

library(DESeq2); library(SummarizedExperiment)
library(dplyr); library(tibble); library(tidyr); library(stringr)
library(ggplot2); library(patchwork); library(ggpubr)
library(pheatmap); library(RColorBrewer); library(ComplexHeatmap); library(circlize)

out_dir <- "multimodal_plots";   dir.create(out_dir, showWarnings = FALSE)
res_dir <- "multimodal_results"; dir.create(res_dir, showWarnings = FALSE)

palette <- c("HER2pos" = "#E64B35", "HER2neg" = "#4DBBD5")
save_plot <- function(p, fname, w=8, h=6, dpi=300){
  ggsave(file.path(out_dir, fname), p, width=w, height=h, dpi=dpi, bg="white")}

# ═══════════════════════════════════════════════════════════
# 1. LOAD & ALIGN
# ═══════════════════════════════════════════════════════════
meta    <- read.csv("metadata_qc_passed.csv", row.names = 1)
vst     <- readRDS("vst_her2_brca.rds")
vst_mat <- assay(vst)

meta_ids <- if ("patient_id" %in% colnames(meta)) meta$patient_id else rownames(meta)

# ═══════════════════════════════════════════════════════════
# 2. BUILD FEATURE MATRIX
# ═══════════════════════════════════════════════════════════

# ── 2A. ERBB2 mRNA ────────────────────────────────────────
erbb2_row <- intersect("ERBB2", rownames(vst_mat))
if (length(erbb2_row) == 0) stop("ERBB2 not found in VST matrix")
erbb2_mRNA <- as.numeric(vst_mat[erbb2_row, ])
names(erbb2_mRNA) <- colnames(vst_mat)

# ── 2B. Amplicon neighbours (chr17q12) ────────────────────
# doi:10.1016/j.molonc.2012.10.012
larger_amplicon_genes <- c("ERBB2","TCAP","PNMT","PERLD1","C17orf37",
                           "NEUROD2","PPP1R1B","STARD3","GRB7","ZNFN1A3")
amp_rows <- intersect(larger_amplicon_genes, rownames(vst_mat))
if (length(amp_rows) == 0) stop("No amplicon genes found in VST matrix")
X_amplicon <- t(vst_mat[amp_rows, ]) %>% as.data.frame()   # samples × genes

# ── 2C. CNV ───────────────────────────────────────────────
# Pull directly from metadata; score: -2 deep del … 2 amplified
meta$ERBB2_CNV <- meta$erbb2_copy_number

# ── 2D. Assemble per-sample molecular data frame ──────────
mol_df <- data.frame(
  sample_id   = meta_ids,
  ERBB2_mRNA  = erbb2_mRNA[meta_ids],
  ERBB2_CNV   = meta$ERBB2_CNV,
  Fraction_Genome_Altered = meta$Fraction_Genome_Altered,
  Mutation_Count = meta$Mutation_Count,
  row.names   = meta_ids
) %>%
  bind_cols(X_amplicon[meta_ids, , drop = FALSE])

# ═══════════════════════════════════════════════════════════
# 3. CLINICAL LABELS
# ═══════════════════════════════════════════════════════════
# IHC_HER2: "Positive" / "Negative" / "Equivocal" / NA
# FISH_HER2: "Positive" / "Negative" / NA
# Build a three-tier label:
#   HER2+    : IHC Positive  OR  (IHC Equivocal AND FISH Positive)
#   HER2-    : IHC Negative  OR  (IHC Equivocal AND FISH Negative)
#   Equivocal: IHC Equivocal AND FISH missing/unavailable

mol_df$IHC_HER2  <- meta[meta_ids, "IHC_HER2"]
mol_df$FISH_HER2 <- meta[meta_ids, "HER2_fish_status"]

mol_df <- mol_df %>%
  mutate(
    clinical_label = case_when(
      IHC_HER2 == "Positive"                                      ~ "HER2+",
      IHC_HER2 == "Negative"                                      ~ "HER2-",
      IHC_HER2 == "Equivocal" & FISH_HER2 == "Positive"          ~ "HER2+",
      IHC_HER2 == "Equivocal" & FISH_HER2 == "Negative"          ~ "HER2-",
      IHC_HER2 == "Equivocal" & is.na(FISH_HER2)                 ~ "Equivocal",
      TRUE                                                         ~ NA_character_
    ),
    clinical_label = factor(clinical_label, levels = c("HER2+","HER2-","Equivocal"))
  )

cat("Clinical label breakdown:\n")
print(table(mol_df$clinical_label, useNA = "ifany"))

# ═══════════════════════════════════════════════════════════
# 4. MOLECULAR CALLS — binarise each layer
# ═══════════════════════════════════════════════════════════
# Thresholds are biologically motivated, not data-driven:
#   CNV >= 1   : gain or amplification (GISTIC2 convention)
#   mRNA       : top quartile within cohort (overexpression)
#   Amplicon   : sample's mean z-score across amplicon genes > 1 SD

mRNA_threshold   <- quantile(mol_df$ERBB2_mRNA, 0.75, na.rm = TRUE)
amplicon_zscore  <- scale(X_amplicon[meta_ids, amp_rows])   # gene-wise z-score
amplicon_mean_z  <- rowMeans(amplicon_zscore, na.rm = TRUE)

mol_df <- mol_df %>%
  mutate(
    CNV_positive      = as.integer(!is.na(ERBB2_CNV) & ERBB2_CNV >= 1),
    mRNA_positive     = as.integer(!is.na(ERBB2_mRNA) & ERBB2_mRNA >= mRNA_threshold),
    amplicon_positive = as.integer(amplicon_mean_z[sample_id] > 1),
    # Multimodal score: sum of positive molecular calls (0–3)
    molecular_score   = CNV_positive + mRNA_positive + amplicon_positive,
    # Multimodal label: positive if >= 2 of 3 layers agree
    molecular_label   = factor(
      ifelse(molecular_score >= 2, "HER2+", "HER2-"),
      levels = c("HER2+", "HER2-")
    )
  )

cat("\nMolecular label breakdown:\n")
print(table(mol_df$molecular_label, useNA = "ifany"))

cat("\nClinical vs molecular concordance:\n")
print(table(Clinical = mol_df$clinical_label,
            Molecular = mol_df$molecular_label, useNA = "ifany"))

# ═══════════════════════════════════════════════════════════
# 5. SAVE MULTIMODAL DEFINITION TABLE
# ═══════════════════════════════════════════════════════════
out_cols <- c("sample_id","IHC_HER2","FISH_HER2","clinical_label",
              "ERBB2_CNV","CNV_positive",
              "ERBB2_mRNA","mRNA_positive",
              "amplicon_positive","molecular_score","molecular_label")
write.csv(mol_df[, out_cols],
          file.path(res_dir, "multimodal_HER2_definition.csv"),
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════
# 6. PLOTS
# ═══════════════════════════════════════════════════════════

# ── 6A. Concordance bar: how many layers agree per sample ──
layer_conc <- mol_df %>%
  filter(!is.na(clinical_label)) %>%
  count(clinical_label, molecular_score) %>%
  mutate(molecular_score = factor(molecular_score, levels = 0:3))

p_score <- ggplot(layer_conc,
                  aes(x = molecular_score, y = n, fill = clinical_label)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("HER2+" = "#E64B35",
                                "HER2-" = "#4DBBD5",
                                "Equivocal" = "#F39B7F")) +
  labs(title    = "Molecular Score (0–3 layers positive) by Clinical Label",
       subtitle  = "Score = CNV + mRNA overexpression + amplicon z-score",
       x = "Number of positive molecular layers", y = "Samples", fill = NULL) +
  theme_bw(base_size = 12)
save_plot(p_score, "01_molecular_score_by_clinical.png", w = 7, h = 5)

# ── 6B. RNA vs CNV scatter, shape = IHC, colour = concordance ──
conc_pal <- c(
  "Concordant +"              = "#E64B35",
  "Concordant -"              = "#4DBBD5",
  "IHC+, molecular-"         = "#984EA3",
  "IHC-, molecular+"         = "#FF7F00",
  "Equivocal"                 = "#F39B7F"
)

conc_df <- mol_df %>%
  filter(!is.na(clinical_label), !is.na(ERBB2_CNV)) %>%
  mutate(
    concordance = case_when(
      clinical_label == "HER2+"     & molecular_label == "HER2+" ~ "Concordant +",
      clinical_label == "HER2-"     & molecular_label == "HER2-" ~ "Concordant -",
      clinical_label == "HER2+"     & molecular_label == "HER2-" ~ "IHC+, molecular-",
      clinical_label == "HER2-"     & molecular_label == "HER2+" ~ "IHC-, molecular+",
      clinical_label == "Equivocal"                               ~ "Equivocal"
    )
  )

p_scatter <- ggplot(conc_df,
                    aes(x = ERBB2_CNV, y = ERBB2_mRNA,
                        colour = concordance, shape = clinical_label)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.75) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = mRNA_threshold, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = conc_pal) +
  scale_shape_manual(values = c("HER2+" = 17, "HER2-" = 16, "Equivocal" = 15)) +
  labs(title    = "ERBB2 RNA vs CNV — Clinical/Molecular Concordance",
       x = "ERBB2 GISTIC2 Score", y = "ERBB2 mRNA (VST)",
       colour = "Concordance", shape = "Clinical label") +
  theme_bw(base_size = 11)
save_plot(p_scatter, "02_RNA_CNV_concordance_scatter.png", w = 9, h = 6)

# ── 6C. Amplicon heatmap — top/bottom molecular score samples ──
# Show amplicon expression for high-score (3) vs low-score (0) samples,
# annotated by clinical label and each molecular call
n_show   <- 50   # samples per group to display
top_samp <- mol_df %>% filter(molecular_score == 3) %>%
  slice_sample(n = n_show) %>% pull(sample_id)
bot_samp <- mol_df %>% filter(molecular_score == 0) %>%
  slice_sample(n = n_show) %>% pull(sample_id)
show_samp <- c(top_samp, bot_samp)

heat_mat <- t(scale(X_amplicon[show_samp, amp_rows]))   # gene-wise z-score
heat_mat[heat_mat >  3] <-  3
heat_mat[heat_mat < -3] <- -3

ann_df <- mol_df[show_samp, c("clinical_label","molecular_score",
                               "CNV_positive","mRNA_positive","amplicon_positive")] %>%
  mutate(across(c(CNV_positive, mRNA_positive, amplicon_positive), factor))

ha <- HeatmapAnnotation(
  df  = ann_df,
  col = list(
    clinical_label    = c("HER2+" = "#E64B35", "HER2-" = "#4DBBD5", "Equivocal" = "#F39B7F"),
    molecular_score   = circlize::colorRamp2(0:3, c("#FFFFFF","#FEE090","#FC8D59","#D73027")),
    CNV_positive      = c("0" = "grey85", "1" = "#3C5488"),
    mRNA_positive     = c("0" = "grey85", "1" = "#E64B35"),
    amplicon_positive = c("0" = "grey85", "1" = "#7E6148")
  )
)

ht <- Heatmap(heat_mat,
  name                      = "Z-score",
  col                       = colorRamp2(c(-3, 0, 3), c("#4DBBD5","white","#E64B35")),
  top_annotation            = ha,
  show_column_names         = FALSE,
  cluster_column_slices     = FALSE,
  column_split              = factor(mol_df[show_samp, "molecular_score"]),
  clustering_method_columns = "ward.D2",
  clustering_method_rows    = "ward.D2",
  column_title              = c("Score 0", "Score 3"),
  row_names_gp              = gpar(fontsize = 9),
  heatmap_legend_param      = list(title = "Z-score")
)

png(file.path(out_dir, "03_amplicon_heatmap_score0_vs_3.png"),
    width = 3000, height = 1000, res = 300)
draw(ht)
dev.off()

# ── 6D. ERBB2 mRNA by GISTIC tier, filled by clinical label ──
mol_df$GISTIC_tier <- factor(
  case_when(
    mol_df$ERBB2_CNV == -2 ~ "Deep del (-2)",
    mol_df$ERBB2_CNV == -1 ~ "Shallow del (-1)",
    mol_df$ERBB2_CNV ==  0 ~ "Neutral (0)",
    mol_df$ERBB2_CNV ==  1 ~ "Gain (1)",
    mol_df$ERBB2_CNV ==  2 ~ "Amplified (2)"
  ),
  levels = c("Deep del (-2)","Shallow del (-1)","Neutral (0)","Gain (1)","Amplified (2)")
)

p_violin <- ggplot(mol_df %>% filter(!is.na(GISTIC_tier), !is.na(clinical_label)),
                   aes(GISTIC_tier, ERBB2_mRNA, fill = clinical_label)) +
  geom_violin(alpha = 0.6, scale = "width") +
  geom_boxplot(width = 0.15, outlier.size = 0.5, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("HER2+" = "#E64B35", "HER2-" = "#4DBBD5",
                                "Equivocal" = "#F39B7F")) +
  labs(title = "ERBB2 mRNA by GISTIC2 Tier",
       x = "GISTIC2 Score", y = "ERBB2 VST Expression", fill = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p_violin, "04_mRNA_by_GISTIC_tier.png", w = 9, h = 5)

# ── 6E. FGA and mutation count ────────────────────────────
p_fga <- ggplot(mol_df %>% filter(!is.na(clinical_label)),
                aes(clinical_label, Fraction_Genome_Altered, fill = clinical_label)) +
  geom_violin(alpha = 0.6) + geom_boxplot(width = 0.12, outlier.size = 0.5) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("HER2+" = "#E64B35", "HER2-" = "#4DBBD5",
                                "Equivocal" = "#F39B7F")) +
  labs(title = "Fraction Genome Altered", x = NULL, y = "FGA") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_mut <- ggplot(mol_df %>% filter(!is.na(clinical_label)),
                aes(clinical_label, log10(Mutation_Count + 1), fill = clinical_label)) +
  geom_violin(alpha = 0.6) + geom_boxplot(width = 0.12, outlier.size = 0.5) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("HER2+" = "#E64B35", "HER2-" = "#4DBBD5",
                                "Equivocal" = "#F39B7F")) +
  labs(title = "Mutation Count (log10)", x = NULL, y = "log10(n+1)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

save_plot(p_fga + p_mut, "05_FGA_mutation_by_label.png", w = 9, h = 5)

cat("\n✓ Plots →", out_dir, "\n✓ Tables →", res_dir, "\n")
