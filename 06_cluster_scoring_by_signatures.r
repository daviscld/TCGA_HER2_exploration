# ============================================================
# Canonical Subtype Marker Annotation
# Maps unsupervised clusters to known breast cancer subtypes:
#   - PAM50 intrinsic subtypes (LumA, LumB, HER2-E, Basal, Normal)
#   - ALTTO HER2+ subtypes (Rediti et al., Nat Commun 2024):
#       Immune-enriched (IM), Proliferative/Metabolic (P/Met),
#       Mesenchymal/Stroma (Mes/S), Luminal (LUM), ERBB2-dependent (ERBB2-D)
#   - TNBC subtypes (Burstein: BLIS, BLIA, LAR, MES)
# Requires: vst_mat, umap_df (with umap_cluster), nmf_df (with nmf_cluster)
# ============================================================

library(DESeq2); library(SummarizedExperiment)
library(dplyr); library(tibble); library(tidyr); library(stringr)
library(ggplot2); library(patchwork)
library(GSVA)             # ssGSEA signature scores
library(ComplexHeatmap); library(circlize)
library(genefu)           # PAM50 calling
library(pheatmap); library(RColorBrewer)

out_dir <- "subtype_annotation"; dir.create(out_dir, showWarnings = FALSE)
res_dir <- "subtype_results";    dir.create(res_dir, showWarnings = FALSE)

save_plot <- function(p, fname, w=8, h=6, dpi=300){
  ggsave(file.path(out_dir, fname), p, width=w, height=h, dpi=dpi, bg="white")}

# ═══════════════════════════════════════════════════════════
# 1. MARKER GENE SIGNATURES
# ═══════════════════════════════════════════════════════════
# Sources:
#   PAM50 genes: Parker et al., JCO 2009
#   ALTTO subtypes: Rediti et al., Nat Commun 2024 (Fig 1)
#   TNBC subtypes: Burstein et al., Clin Cancer Research 2015; 
#                  Lehmann et al., PLoS ONE 2016
#   Immune markers: standard TIL/checkpoint literature
#   Proliferation: MKI67, PCNA, TOP2A (consensus proliferation genes)

subtype_markers <- list(

  # ── PAM50 canonical genes ─────────────────────────────────
  # High in Luminal A: low proliferation, strong ER signalling
  LuminalA = c("ESR1","FOXA1","GATA3","XBP1","TFF3","MLPH","NAT1",
               "SLC39A6","MAPT","PGR","BCL2","CCND1","SFRP1"),

  # Luminal B: ER+ but higher proliferation, some HER2 signal
  LuminalB = c("ESR1","FOXA1","MKI67","CCNE1","PCNA","CDC20",
               "AURKA","CDK1","BUB1","RRM2","BIRC5","MYBL2"),

  # HER2-Enriched (PAM50): high ERBB2 + amplicon + proliferation
  # ~50% of IHC HER2+ are HER2-E; also appears in IHC HER2- tumours
  HER2_Enriched = c("ERBB2","GRB7","STARD3","PERLD1","FGFR4",
                    "TMEM45B","GPR160","CDK12","PSMD3","MED24",
                    "MKI67","AURKA","CCNB1"),

  # Basal-like: TP53 mutant, EGFR+, low ER/PR/HER2, CK5/6+
  Basal = c("KRT5","KRT14","KRT17","EGFR","FOXC1","SFRP1",
            "ANXA8","FAR2","TRIM29","MIA","LAMC2","COL17A1",
            "CDH3","SERPINB5","ACTN1"),

  # ── ALTTO HER2+ subtypes (Rediti 2024) ────────────────────
  # Immune-enriched: TIL-high, checkpoint+, best prognosis
  ALTTO_Immune = c("CD8A","CD8B","GZMB","PRF1","CXCL9","CXCL10",
                   "CD274","PDCD1","PDCD1LG2","CD19","MS4A1",
                   "TIGIT","LAG3","HAVCR2","CTLA4","FOXP3",
                   "IFNG","IRF1","STAT1"),

  # Proliferative/Metabolic: glycolysis, cholesterol, MYC
  ALTTO_Proliferative = c("MKI67","TOP2A","AURKA","AURKB","BUB1",
                          "CDK1","CCNB1","CCNB2","MCM2","MCM6",
                          "LDHA","PKM","HK2","ACACA","FASN",
                          "SQLE","HMGCR","FDFT1","PCSK9","MYC"),

  # Mesenchymal/Stroma-enriched: EMT, TGF-β, angiogenesis, worst prognosis
  ALTTO_Mesenchymal = c("VIM","FN1","CDH2","ZEB1","ZEB2","TWIST1",
                        "SNAI1","SNAI2","TGFB1","TGFB2","TGFBR2",
                        "COL1A1","COL3A1","FAP","PDGFRB","ACTA2",
                        "VEGFA","ANGPT2","PDGFRA","THY1","S100A4"),

  # Luminal (ALTTO): HR+, good prognosis, low pCR
  ALTTO_Luminal = c("ESR1","PGR","FOXA1","GATA3","XBP1","TFF1",
                    "TFF3","AREG","PROM1","MYB","CCND1","BCL2",
                    "AGR2","ALCAM","CELSR2"),

  # ERBB2-Dependent: high ERBB2 + amplicon, high pCR with anti-HER2
  ALTTO_ERBB2dep = c("ERBB2","GRB7","STARD3","TCAP","PNMT",
                     "MIEN1","NEUROD2","PPP1R1B","ZNFN1A3",
                     "ERBB3","ERBB4","MUC4","EGFR","CDK12"),

  # ── TNBC subtypes (Burstein/Lehmann) ──────────────────────
  # Basal-like Immune-Activated (BLIA): immune high, best TNBC prognosis
  TNBC_BLIA = c("STAT1","IRF1","OAS1","MX1","IFIT1","ISG15",
                "CXCL10","CXCL9","CD274","JAK1","JAK2"),

  # Basal-like Immune-Suppressed (BLIS): low immune, FOXC1 high
  TNBC_BLIS = c("FOXC1","ANXA8","FAR2","COL17A1","CDKN2A",
                "RB1","BRCA1","BRCA2","PALB2","RAD51"),

  # Luminal Androgen Receptor (LAR): AR+, PI3K mutations
  TNBC_LAR  = c("AR","DHCR24","ALCAM","APOD","PIP","FSIP1",
                "LGR5","FKBP4","TMPRSS2","SPDEF","FOXA1"),

  # Mesenchymal (MES): EMT, stem-like, WNT/Notch
  TNBC_MES  = c("VIM","CDH2","FN1","TWIST2","AXL","PDGFRB",
                "WNT5A","NOTCH1","NOTCH3","SOX9","SPARC"),

  # ── Standalone biomarkers ─────────────────────────────────
  Proliferation = c("MKI67","PCNA","TOP2A","MCM2","CDK1",
                    "CCNB1","BUB1B","AURKA","BIRC5","RRM2"),

  ER_signalling = c("ESR1","PGR","FOXA1","GATA3","XBP1",
                    "TFF1","TFF3","CCND1","BCL2","AREG"),

  HER2_amplicon = c("ERBB2","GRB7","STARD3","TCAP","PNMT",
                    "PERLD1","C17orf37","NEUROD2","PPP1R1B"),

  Immune_general = c("CD8A","CD4","FOXP3","CD274","PDCD1",
                     "CD19","CD56","GZMB","PRF1","IFNG"),

  EMT           = c("VIM","CDH2","FN1","TWIST1","ZEB1",
                    "SNAI1","CDH1","OCLN","CLDN4","EPCAM"),

  Angiogenesis  = c("VEGFA","VEGFC","ANGPT2","PDGFRB",
                    "FLT1","KDR","TEK","NRP1","HIF1A")
)

# ═══════════════════════════════════════════════════════════
# 3. ssGSEA SIGNATURE SCORES
# ═══════════════════════════════════════════════════════════
# Score every sample on every marker gene set.
# ssGSEA is rank-based so it's robust to batch/normalisation differences.

# Get the samples used in clustering
all_cluster_samples <- nmf_df$patient_id

# Filter each signature to genes present in the expression matrix
sigs_present <- lapply(subtype_markers, function(g)
  intersect(g, rownames(vst_mat))
)
sigs_present <- Filter(function(g) length(g) >= 3, sigs_present)

cat(sprintf("\nRunning ssGSEA on %d signatures...\n", length(sigs_present)))
gsva_scores <- gsva(
  expr      = vst_mat[, all_cluster_samples],
  gset.idx.list = sigs_present,
  method    = "ssgsea",
  ssgsea.norm = TRUE,
  verbose   = FALSE
)
# gsva_scores: signatures × samples

# ═══════════════════════════════════════════════════════════
# 4. ATTACH SCORES + CLUSTER LABELS
# ═══════════════════════════════════════════════════════════

score_df <- as.data.frame(t(gsva_scores)) %>%
  rownames_to_column("patient_id") %>%
  left_join(nmf_df  %>% dplyr::select(patient_id, nmf_cluster, 
                clinical_label, ER_label, PR_label),
            by = "patient_id")

write.csv(score_df, file.path(res_dir, "signature_scores_all_samples.csv"),
          row.names = FALSE)

# ═══════════════════════════════════════════════════════════
# 5. PLOTS
# ═══════════════════════════════════════════════════════════

sig_names <- names(sigs_present)

# ── 5A. Signature score heatmap — UMAP clusters ───────────
plot_sig_heatmap <- function(score_df, cluster_col, title, fname) {
  df_filt <- score_df %>% filter(!is.na(.data[[cluster_col]]),
                                  .data[[cluster_col]] != "Noise")
  mat <- df_filt %>% dplyr::select(all_of(sig_names)) %>%
    as.matrix() %>% t()
  colnames(mat) <- df_filt$sample_id

  # Z-score across samples per signature
  mat_z <- t(scale(t(mat)))
  mat_z[mat_z >  2] <-  2; mat_z[mat_z < -2] <- -2

  cluster_vec <- df_filt[[cluster_col]]
  n_cl <- length(unique(cluster_vec))
  cl_cols <- setNames(
    RColorBrewer::brewer.pal(max(n_cl, 3), "Set2")[1:n_cl],
    sort(unique(cluster_vec))
  )

  ha <- HeatmapAnnotation(
    Cluster = cluster_vec,
    HER2    = df_filt$clinical_label,
    ER      = df_filt$ER_label,
    PR      = df_filt$PR_label,
    col     = list(
      Cluster = cl_cols,
      HER2    = c("HER2+" = "#E64B35", "HER2-" = "#4DBBD5")
    )
  )

  ht <- Heatmap(mat_z, name = "ssGSEA\nZ-score",
    col = colorRamp2(c(-2, 0, 2), c("#4DBBD5","white","#E64B35")),
    top_annotation    = ha,
    show_column_names = FALSE,
    column_split      = cluster_vec,
    cluster_column_slices = FALSE,
    clustering_method_rows = "ward.D2",
    row_names_gp      = gpar(fontsize = 8),
    column_title      = title,
    heatmap_legend_param = list(title = "Z-score")
  )

  png(file.path(out_dir, fname), width = 1100, height = 900, res = 300)
  draw(ht); dev.off()

  #also plot pdf for use in Illustrator to make final figure panels
    pdf(file.path(out_dir, sub("\\.png$", ".pdf", fname)), 
    width = 8, height = 9)
    draw(ht); dev.off()
}

plot_sig_heatmap(score_df, "nmf_cluster",
                 "Subtype Signatures — NMF Clusters",
                 "02_signatures_NMF_clusters.png")

# ── 5B. UMAP coloured by key individual signature scores ──
key_sigs <- c("HER2_amplicon","ER_signalling","Basal",
              "LuminalA","Proliferation","ALTTO_ERBB2dep",
              "LuminalB", "HER2_Enriched")
key_sigs <- intersect(key_sigs, sig_names)

umap_score_df <- umap_df %>%
  left_join(score_df %>% dplyr::select(patient_id, nmf_cluster, 
  all_of(key_sigs)),
            by = "patient_id")

plots_umap_sig <- lapply(key_sigs, function(sig) {
  ggplot(umap_score_df, aes(UMAP1, UMAP2, colour = .data[[sig]])) +
    geom_point(size = 1, alpha = 0.7) +
    scale_colour_distiller(palette = "RdYlBu", direction = -1) +
    labs(title = sig, colour = "ssGSEA") +
    theme_bw(base_size = 9) +
    theme(legend.key.height = unit(0.3, "cm"))
})
#add an additional plot showing the NMF clusters for reference
plots_umap_sig[[length(plots_umap_sig) + 1]] <- ggplot(umap_score_df, 
    aes(UMAP1, UMAP2, colour = nmf_cluster)) +
  geom_point(size = 1, alpha = 0.7) +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "NMF Clusters", colour = NULL) +
  theme_bw(base_size = 9) +
  theme(legend.key.height = unit(0.3, "cm"))

save_plot(
  wrap_plots(plots_umap_sig, ncol = 3),
  "04_UMAP_signature_scores.png", w = 14, h = 8
)

# ── 5D. Individual canonical marker genes — dot plot ──────
# Show mean expression + % expressing for key individual genes
canonical_genes <- c(
  "ERBB2","ESR1","PGR","MKI67",          # clinical biomarkers
  "KRT5","EGFR","FOXC1",                  # basal
  "FOXA1","GATA3","XBP1",                 # luminal
  "GRB7","STARD3",                        # HER2 amplicon
  "CD8A","CD274","FOXP3",                 # immune
  "VIM","CDH2","FN1",                     # EMT/mesenchymal
  "VEGFA","TGFB1"                         # angiogenesis/TGF-beta
)
canonical_genes <- intersect(canonical_genes, rownames(vst_mat))

make_dotplot <- function(score_df, cluster_col, method_label) {
  cl_sym <- sym(cluster_col)   # convert string to symbol once, reuse throughout

  samp_order <- score_df %>%
    filter(!is.na(!!cl_sym), !!cl_sym != "Noise") %>%
    pull(patient_id)

  expr_long <- vst_mat[canonical_genes, samp_order] %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "patient_id", values_to = "expr") %>%
    left_join(score_df %>% dplyr::select(patient_id, !!cl_sym),
              by = "patient_id") %>%
    rename(cluster = !!cl_sym) %>%
    group_by(gene, cluster) %>%
    summarise(
      mean_expr = mean(expr, na.rm = TRUE),
      pct_expr  = 100 * mean(expr > quantile(vst_mat[gene, ], 0.25)),
      .groups   = "drop"
    ) %>%
    mutate(gene = factor(gene, levels = canonical_genes))

  ggplot(expr_long, aes(x = cluster, y = gene,
                         size = pct_expr, colour = mean_expr)) +
    geom_point() +
    scale_colour_distiller(palette = "RdYlBu", direction = -1) +
    scale_size(range = c(1, 7)) +
    labs(title  = paste(method_label, "— Canonical Marker Dot Plot"),
         x = NULL, y = NULL,
         colour = "Mean VST", size = "% Above Q1") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

save_plot(make_dotplot(score_df, "nmf_cluster",  "NMF Clusters"),
          "08_dotplot_markers_NMF.png",  w = 8, h = 8)

# ── 5E. Proposed cluster label table ──────────────────────
# Print a summary to help you manually assign labels
cat("\n═══════════════════════════════════════════════════\n")
cat("  Cluster Annotation Summary\n")
cat("  Use this to assign biological labels to clusters\n")
cat("═══════════════════════════════════════════════════\n\n")

cat(sprintf("── NMF Clusters ──\n"))
summ <- score_df %>%
    filter(!is.na(nmf_cluster), nmf_cluster != "Noise") %>%
    group_by(nmf_cluster) %>%
    summarise(
        N = n(),
        pct_HER2pos = round(100 * mean(clinical_label == "HER2+", na.rm = TRUE), 1),
        ERBB2_score = round(mean(HER2_amplicon, na.rm = TRUE), 2),
        ER_score = round(mean(ER_signalling, na.rm = TRUE), 2),
        Immune_score = round(mean(ALTTO_Immune, na.rm = TRUE), 2),
        Mesenchymal_score = round(mean(ALTTO_Mesenchymal, na.rm = TRUE), 2),
        Prolif_score = round(mean(Proliferation, na.rm = TRUE), 2),
        .groups = "drop"
    )
print(as.data.frame(summ), row.names = FALSE)
cat("\n")
write.csv(summ, file.path(res_dir, "nmf_cluster_annotation_summary.csv"),
    row.names = FALSE)

cat("✓ Plots →", out_dir, "\n✓ Tables →", res_dir, "\n")
