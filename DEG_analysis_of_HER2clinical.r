#Differential expression analysis of HER2+ vs HER2- breast cancer 
#samples from TCGA, using DESeq2.
#This script loads the QC'd DESeq2 dataset and metadata, 
#performs DEG analysis, and saves results for downstream use.

# ============================================================
# DESeq2 Full Analysis Pipeline
# Input: dds_her2_brca.rds from previous QC script
# ============================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ashr)           # for lfcShrink method = "ashr"
library(scales)
library(patchwork)

# ── 0. CONFIG ──────────────────────────────────────────────
plot_dir   <- "deseq2_plots"
res_dir    <- "deseq2_results"
dir.create(plot_dir, showWarnings = FALSE)
dir.create(res_dir,  showWarnings = FALSE)

contrast   <- c("her2_clinical", "HER2pos", "HER2neg") 
p_thresh   <- 0.05
lfc_thresh <- 0.58        # log2 fold-change threshold for "significant" (~1.5x)
top_n_heat <- 50         # top N DEGs for heatmap
top_label  <- 20         # genes to label on volcano

# Helper: save plots consistently
save_plot <- function(p, fname, w=8, h=6, dpi=300) {
  ggsave(file.path(plot_dir, fname), p, width=w, height=h, dpi=dpi,
         bg = "white")
  invisible(p)
}

# ── 1. LOAD & RUN DESeq2 ───────────────────────────────────
dds <- readRDS("dds_her2_brca.rds")
vst <- readRDS("vst_her2_brca.rds")

# Run the full DESeq2 pipeline
# (estimateSizeFactors already called in QC script, so safe to re-run)
# Already assigned design (just her2_clinical), so no need to specify here
dds <- DESeq(dds, parallel = FALSE)
# ── 2. DISPERSION & FIT QC PLOTS ───────────────────────────

# 2a. Dispersion estimates (looking for decrease as mean increases)
png(file.path(plot_dir, "01_dispersion_estimates.png"),
    width = 800, height = 600, res = 120)
plotDispEsts(dds,
             main = "Dispersion Estimates",
             cex  = 0.4,
             genecol = "steelblue",
             fitcol  = "firebrick",
             finalcol = "black")
dev.off()

# 2b. Size factors — should be centered around 1 (looking for above 0.5, below 2)
sf_df <- data.frame(
  sample     = colnames(dds),
  size_factor = sizeFactors(dds),
  her2       = colData(dds)$her2_clinical
)

p_sf <- ggplot(sf_df, aes(x = reorder(sample, size_factor),
                           y = size_factor, fill = her2)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "firebrick") +
  scale_fill_manual(values = c("HER2pos" = "#E64B35", "HER2neg" = "#4DBBD5")) +
  labs(title = "DESeq2 Size Factors per Sample",
       x = "Samples", y = "Size Factor", fill = "HER2 Status") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_plot(p_sf, "02_size_factors.png", w=10, h=4)

# ── 3. SAMPLE-LEVEL QC (VST-based) ─────────────────────────

# 3a. PCA — PC1/PC2 and PC2/PC3
#plotPCA automatically uses the top 500 most variable genes
pca_data <- plotPCA(vst, intgroup = c("her2_clinical"), returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca12 <- ggplot(pca_data, aes(PC1, PC2, color = her2_clinical)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(group = her2_clinical), linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("HER2pos"="#E64B35","HER2neg"="#4DBBD5")) +
  xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  labs(title = "PCA — VST Counts", color = "HER2 Status") +
  theme_bw(base_size = 12)
save_plot(p_pca12, "03a_pca_PC1_PC2.png")

# PC2 vs PC3 - prcomp uses all genes, so percent variance will differ
pca_all <- prcomp(t(assay(vst)))
pca_df  <- as.data.frame(pca_all$x[, 1:5]) %>%
  rownames_to_column("sample") %>%
  left_join(as.data.frame(colData(vst)) %>% rownames_to_column("sample"),
            by = "sample")

pct_all <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)

p_pca23 <- ggplot(pca_df, aes(PC2, PC3, color = her2_clinical)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(group = her2_clinical), linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("HER2pos"="#E64B35","HER2neg"="#4DBBD5")) +
  xlab(paste0("PC2: ", pct_all[2], "% variance")) +
  ylab(paste0("PC3: ", pct_all[3], "% variance")) +
  labs(title = "PCA — PC2 vs PC3 (VST Counts)", color = "HER2 Status") +
  theme_bw(base_size = 12)
save_plot(p_pca23, "03b_pca_PC2_PC3.png")

# 3b. Scree plot (variance explained per PC)
scree_df <- data.frame(PC = 1:15, var_exp = pct_all[1:15])
p_scree  <- ggplot(scree_df, aes(PC, var_exp)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_line(color = "firebrick", linewidth = 0.8) +
  geom_point(color = "firebrick", size = 2) +
  scale_x_continuous(breaks = 1:15) +
  labs(title = "Scree Plot — Variance Explained per PC",
       x = "Principal Component", y = "% Variance Explained") +
  theme_bw(base_size = 12)
save_plot(p_scree, "04_scree_plot.png", w=7, h=4)

# 3c. Sample-sample distance heatmap (did this before HER2 clinical annotation,
#but now repeat it, and include PR+ and ER+ status in the annotation)
samp_dists <- dist(t(assay(vst)))
samp_mat   <- as.matrix(samp_dists)

anno_col <- as.data.frame(colData(vst)) %>%
    dplyr::select(her2_clinical, PR_status_by_ihc, ER_Status_By_IHC) %>%
    rename(HER2 = her2_clinical, PR = PR_status_by_ihc,
                 ER = ER_Status_By_IHC) %>%
    mutate(across(everything(), factor))
anno_colors <- list(HER2 = c("HER2pos" = "#E64B35", "HER2neg" = "#4DBBD5"),
                    PR = c("Positive" = "#E64B35", 
                            "Negative" = "#4DBBD5",
                            "Indeterminate" = "#999999"),
                    ER = c("Positive" = "#E64B35", 
                            "Negative" = "#4DBBD5",
                            "Indeterminate" = "#999999"))

png(file.path(plot_dir, "05_sample_distance_heatmap.png"),
    width = 1000, height = 900, res = 130)
pheatmap(samp_mat,
        clustering_distance_rows = samp_dists,
        clustering_distance_cols = samp_dists,
        annotation_col = anno_col,
        annotation_row = anno_col,
        annotation_colors = anno_colors,
        color = colorRampPalette(c("#2166AC", "white", "#D6604D"))(100),
        show_rownames = FALSE,
        show_colnames = FALSE,
        main = "Sample-Sample Euclidean Distance (VST)",
        fontsize = 9)
dev.off()

# 3d. Cook's distances — outlier samples
# High Cook's distance flags samples that strongly influence model fit
cook_mat  <- log10(assays(dds)[["cooks"]])
cook_max  <- apply(cook_mat, 2, max, na.rm = TRUE)
cook_df   <- data.frame(
  sample   = names(cook_max),
  max_cooks = cook_max,
  her2      = colData(dds)$her2_clinical
)

p_cooks <- ggplot(cook_df, aes(x = reorder(sample, max_cooks),
                                y = max_cooks, color = her2)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_hline(yintercept = log10(0.01 / nrow(dds)),  # rough threshold
             linetype = "dashed", color = "firebrick") +
  scale_color_manual(values = c("HER2pos"="#E64B35","HER2neg"="#4DBBD5")) +
  labs(title = "Max Cook's Distance per Sample",
       subtitle = "Dashed line: potential outlier threshold",
       x = NULL, y = "log10(max Cook's distance)", color = "HER2") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_plot(p_cooks, "06_cooks_distances.png", w=10, h=4)

# ── 4. EXTRACT & SHRINK RESULTS ────────────────────────────
# Raw results (for MA plot showing shrinkage effect)
res_raw <- results(dds,
                   contrast  = contrast,
                   alpha     = p_thresh,
                   pAdjustMethod = "BH")

summary(res_raw)

# Shrunken LFC (ashr — preserves sign, handles outliers well)
res_shr <- lfcShrink(dds,
                     coef = "her2_clinical_HER2pos_vs_HER2neg",
                     res      = res_raw,
                     type     = "ashr")

summary(res_shr)

# Convert to tidy data frame, add gene symbols
res_df <- as.data.frame(res_shr) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    sig    = case_when(
      padj < p_thresh & log2FoldChange >  lfc_thresh ~ "Up",
      padj < p_thresh & log2FoldChange < -lfc_thresh ~ "Down",
      TRUE ~ "NS"
    ),
    sig = factor(sig, levels = c("Up","Down","NS"))
  ) %>%
  arrange(padj)

# Save full results table
write.csv(res_df, file.path(res_dir, "DESeq2_HER2pos_vs_HER2neg_full.csv"),
          row.names = FALSE)

# Significant DEGs only
deg_df <- filter(res_df, sig != "NS")
write.csv(deg_df, file.path(res_dir, "DESeq2_HER2pos_vs_HER2neg_DEGs.csv"),
          row.names = FALSE)

cat(sprintf("\n── DEG Summary ──\n"))
cat(sprintf("  Upregulated   (HER2pos): %d\n", sum(res_df$sig=="Up",   na.rm=TRUE)))
cat(sprintf("  Downregulated (HER2pos): %d\n", sum(res_df$sig=="Down", na.rm=TRUE)))
cat(sprintf("  Total DEGs             : %d\n", nrow(deg_df)))

# ── 5. MA PLOTS (raw vs shrunken) ──────────────────────────

make_ma <- function(res, title) {
  df <- as.data.frame(res) %>%
    mutate(sig = !is.na(padj) & padj < p_thresh,
           baseMean_log = log10(baseMean + 1))
  ggplot(df, aes(baseMean_log, log2FoldChange, color = sig)) +
    geom_point(size = 0.5, alpha = 0.4) +
    geom_hline(yintercept = 0, color = "firebrick", linewidth = 0.7) +
    geom_hline(yintercept = c(-lfc_thresh, lfc_thresh),
               linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("TRUE"="#E64B35","FALSE"="grey70"),
                       labels = c("TRUE"="Significant","FALSE"="NS")) +
    labs(title = title,
         x = "log10(mean normalised count + 1)",
         y = "log2 Fold Change (HER2+ vs HER2–)",
         color = NULL) +
    ylim(-8, 8) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

#visualize the effect of shrinkage on the MA plot
p_ma_raw <- plotMA(res_raw)
p_ma_shr <- plotMA(res_shr)

png(file.path(plot_dir, "07_MA_plots_raw_vs_shrunken.png"),
    width = 1000, height = 2400, res = 300)
par(mfrow = c(2, 1))
plotMA(res_raw, main = "Raw LFC")
plotMA(res_shr, main = "ashr Shrunken LFC")
dev.off()

# ── 6. VOLCANO PLOT ────────────────────────────────────────
# Label top DEGs by |LFC| among significant hits
# Also label most significant DEGs by padj, to capture any with smaller effect 
# sizes but very strong evidence
top_labels <- res_df %>%
    filter(sig != "NS") %>%
    slice_max(abs(log2FoldChange), n = top_label) %>%
    bind_rows(
        res_df %>%
            filter(sig != "NS") %>%
            slice_min(padj, n = 10)
    ) %>%
    pull(gene_id) %>%
    unique()

p_volcano <- EnhancedVolcano(res_df,
  lab             = res_df$gene_id,
  x               = "log2FoldChange",
  y               = "padj",
  pCutoff         = p_thresh,
  FCcutoff        = lfc_thresh,
  selectLab       = top_labels,
  labSize         = 5,
  drawConnectors  = TRUE,
  widthConnectors = 0.4,
  colConnectors   = "grey50",
  col             = c("grey70","grey70","steelblue","#E64B35"),
  colAlpha        = 0.6,
  pointSize       = 2,
  title           = "HER2+ vs HER2– (DESeq2, ashr shrinkage)",
  subtitle        = paste0("padj < ", p_thresh, "  |LFC| > ", lfc_thresh),
  caption         = paste0("Total DEGs: ", nrow(deg_df)),
  legendPosition  = "bottom"
)
save_plot(p_volcano, "08_volcano.png", w=10, h=9)

# ── 7. P-VALUE HISTOGRAM ───────────────────────────────────
# Should show anti-conservative spike near 0; flat for null genes.
# A U-shape or spike at 1 indicates a model problem.
p_phist <- ggplot(res_df %>% filter(!is.na(pvalue)),
                  aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white", boundary = 0) +
  labs(title = "P-value Distribution (before adjustment)",
       subtitle = "Spike near 0 expected; uniform tail = model fits well",
       x = "Raw p-value", y = "Count") +
  theme_bw(base_size = 12)
save_plot(p_phist, "09_pvalue_histogram.png", w=6, h=4)

# ── 8. LFC DISTRIBUTION ────────────────────────────────────
p_lfc <- ggplot(res_df %>% filter(!is.na(log2FoldChange)),
                aes(x = log2FoldChange, fill = sig)) +
  geom_histogram(bins = 80, color = "white", boundary = 0) +
  scale_fill_manual(values = c("Up"="#E64B35","Down"="#4DBBD5","NS"="grey70")) +
  geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
             linetype = "dashed", color = "black") +
  labs(title = "Log2 Fold-Change Distribution",
       x = "log2 Fold Change (HER2+ vs HER2–)", y = "Count",
       fill = "DEG Status") +
  theme_bw(base_size = 12)
save_plot(p_lfc, "10_lfc_distribution.png", w=7, h=4)

# ── 9. TOP DEG HEATMAP ─────────────────────────────────────
top_degs <- res_df %>%
  filter(sig != "NS") %>%
  slice_max(abs(log2FoldChange), n = top_n_heat) %>%
  pull(gene_id)

# Z-score VST matrix for top DEGs
vst_mat <- assay(vst)[top_degs, ] %>%
  t() %>% scale() %>% t()    # gene-wise z-score

# Cap extreme z-scores for visual clarity
vst_mat[vst_mat >  3] <-  3
vst_mat[vst_mat < -3] <- -3

# Gene symbol labels
gene_labels <- res_df %>%
  filter(gene_id %in% top_degs) %>%
  column_to_rownames("gene_id")

# Sample annotation
ha_col <- HeatmapAnnotation(
  HER2    = colData(vst)$her2_clinical,
  col     = list(HER2 = c("HER2pos"="#E64B35","HER2neg"="#4DBBD5")),
  annotation_name_side = "left"
)

# Row (gene) annotation — LFC bar
ha_row <- rowAnnotation(
  LFC = anno_barplot(
    gene_labels[rownames(vst_mat), "log2FoldChange"],
    gp = gpar(fill = ifelse(
      gene_labels[rownames(vst_mat), "log2FoldChange"] > 0,
      "#E64B35", "#4DBBD5"))
  ),
  annotation_width = unit(2, "cm")
)

row_labels <- rownames(vst_mat)

ht <- Heatmap(vst_mat,
  name                  = "Z-score",
  col                   = colorRamp2(c(-3, 0, 3),
                                     c("#2166AC","white","#D6604D")),
  top_annotation        = ha_col,
  right_annotation      = ha_row,
  row_labels            = row_labels,
  row_names_gp          = gpar(fontsize = 7),
  show_column_names     = FALSE,
  column_split          = colData(vst)$her2_clinical,
  cluster_column_slices = FALSE,
  clustering_method_rows    = "ward.D2",
  clustering_method_columns = "ward.D2",
  column_title          = c("HER2–", "HER2+"),
  row_title             = paste0("Top ", top_n_heat, " DEGs by |LFC|"),
  heatmap_legend_param  = list(title = "Z-score\n(VST)")
)

png(file.path(plot_dir, "11_top_DEG_heatmap.png"),
    width = 1100, height = 1000, res = 300)
draw(ht)
dev.off()

# ── 10. INDIVIDUAL GENE BOXPLOTS (top 12 DEGs) ─────────────
top12 <- res_df %>%
  filter(sig != "NS") %>%
  slice_max(abs(log2FoldChange), n = 12) %>%
  pull(gene_id)

norm_counts <- counts(dds, normalized = TRUE)

box_df <- norm_counts[top12, ] %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    tidyr::pivot_longer(-gene_id, names_to = "sample", values_to = "norm_count") %>%
    left_join(res_df %>% dplyr::select(gene_id, log2FoldChange, padj),
                        by = "gene_id") %>%
    left_join(as.data.frame(colData(dds)) %>%
                            rownames_to_column("sample") %>%
                            dplyr::select(sample, her2_clinical),
                        by = "sample") %>%
    mutate(
        label = sprintf("%s\nLFC=%.2f\npadj=%.1e", gene_id, log2FoldChange, padj),
        log_count = log10(norm_count + 1)
    )

p_violin <- ggplot(box_df, aes(x = her2_clinical, y = log_count, 
                               color = her2_clinical)) +
  geom_violin(alpha = 0.7, size = 0.4) +
  geom_boxplot(width = 0.15, alpha = 0.8, size = 0.3, 
               outlier.size = 1) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.4) +
  scale_color_manual(values = c("HER2pos" = "#E64B35", 
                               "HER2neg" = "#4DBBD5")) +
  facet_wrap(~label, scales = "free_y", ncol = 4) +
  labs(title = "Top 12 DEGs — Normalised Counts",
       x = NULL, y = "log10(normalised count + 1)", 
       color = "HER2 Status") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_violin, "12_top12_gene_violins.png", w = 12, h = 8)

# ── 11. GENE ONTOLOGY ENRICHMENT ───────────────────────────
run_go <- function(genes_up, genes_dn, ont = "BP") {

  ego_up <- enrichGO(gene         = genes_up,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "SYMBOL",
                     ont          = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable     = TRUE)

  ego_dn <- enrichGO(gene         = genes_dn,
                     OrgDb        = org.Hs.eg.db,
                     keyType      = "SYMBOL",
                     ont          = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable     = TRUE)

  list(up = ego_up, down = ego_dn)
}
    
gene_up <- res_df %>% filter(sig == "Up",   
            !is.na(gene_id)) %>% pull(gene_id) %>% unique()
gene_dn <- res_df %>% filter(sig == "Down",  
            !is.na(gene_id)) %>% pull(gene_id) %>% unique()

#I don't generally love GO results but test anyways
#biological process and molecular function ontologies
go_bp <- run_go(gene_up, gene_dn, ont = "BP")
go_mf <- run_go(gene_up, gene_dn, ont = "MF")
#not currently using molecular function results

# GO dotplot — upregulated genes
if (!is.null(go_bp$up) && nrow(go_bp$up) > 0) {
  p_go_up <- dotplot(go_bp$up, showCategory = 20, font.size = 9) +
    labs(title = "GO Biological Process — Upregulated in HER2+") +
    theme_bw(base_size = 10)
  save_plot(p_go_up, "13a_GO_BP_upregulated.png", w=9, h=8)
}

# GO dotplot — downregulated genes
if (!is.null(go_bp$down) && nrow(go_bp$down) > 0) {
  p_go_dn <- dotplot(go_bp$down, showCategory = 20, font.size = 9) +
    labs(title = "GO Biological Process — Downregulated in HER2+") +
    theme_bw(base_size = 10)
  save_plot(p_go_dn, "13b_GO_BP_downregulated.png", w=9, h=8)
}

# ── 12. GSEA (pre-ranked on shrunken LFC × -log10 padj) ────
# Using Hallmark gene sets — best for clean biological interpretation
gsea_rank <- res_df %>%
  filter(!is.na(gene_id), !is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(stat = log2FoldChange) %>%
  arrange(desc(stat)) %>%
  dplyr::select(gene_id, stat) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  deframe()   # named vector

# Download Hallmark gene sets (MSigDB via msigdbr)
if (!requireNamespace("msigdbr", quietly = TRUE)) install.packages("msigdbr")
library(msigdbr)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol) %>%
    rename(entrez_gene = gene_symbol)

gsea_res <- GSEA(geneList    = gsea_rank,
                 TERM2GENE   = hallmark,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 minGSSize     = 15,
                 maxGSSize     = 500,
                 seed          = 42,
                 verbose       = FALSE)

# GSEA dotplot
if (nrow(as.data.frame(gsea_res)) > 0) {
    p_gsea <- dotplot(gsea_res, showCategory = 20, split = ".sign",
                                        font.size = 8) +
        facet_grid(~.sign, 
        labeller = labeller(.sign = c("activated" = "Activated in HER2+ samples",
                    "suppressed" = "Suppressed in HER2+ samples"))) +
        labs(title = "GSEA — Hallmark Gene Sets") +
        theme_bw(base_size = 9)
    save_plot(p_gsea, "14_GSEA_hallmark.png", w = 12, h = 8)

    # GSEA enrichment plot for top activated / suppressed pathway
    top_act  <- gsea_res@result %>% filter(NES > 0) %>% 
                slice_max(NES, n = 1) %>% pull(ID)
    top_supp <- gsea_res@result %>% filter(NES < 0) %>% 
                slice_min(NES, n = 1) %>% pull(ID)

    for (pw in c(top_act, top_supp)) {
        if (length(pw) > 0) {
            fname <- paste0("15_GSEA_enrichplot_",
                      str_replace_all(pw, "[^A-Za-z0-9]", "_"), ".png")
            p_ep  <- enrichplot::gseaplot2(gsea_res, geneSetID = pw,
                        title = pw, pvalue_table = TRUE)
            save_plot(p_ep, fname, w = 9, h = 6)
        }
    }

    write.csv(as.data.frame(gsea_res),
                        file.path(res_dir, "GSEA_hallmark_results.csv"),
                        row.names = FALSE)
}

# ── 13. SUMMARY FIGURE ─────────────────────────────────────
# Combine key numbers into a text-summary ggplot for reports
n_up   <- sum(res_df$sig == "Up",   na.rm = TRUE)
n_dn   <- sum(res_df$sig == "Down", na.rm = TRUE)
n_tot  <- nrow(dds)
n_test <- sum(!is.na(res_df$padj))

summary_text <- data.frame(
  label = c("Genes tested", "Significant (padj<0.05, |LFC|>0.58)",
            "  ↑ Upregulated in HER2+", "  ↓ Downregulated in HER2+"),
  value = c(n_test, n_up + n_dn, n_up, n_dn)
)

p_summary <- ggplot(summary_text, aes(x = 1, y = rev(seq_len(nrow(summary_text))),
                                       label = paste0(label, ":  ", value))) +
  geom_text(hjust = 0, size = 4.5) +
  xlim(1, 2.5) +
  labs(title = "DESeq2 Results Summary — HER2+ vs HER2–") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face="bold", size=13, hjust=0.1))
save_plot(p_summary, "00_results_summary.png", w=6, h=3)

cat("\n✓ All plots saved to:", plot_dir, "\n")
cat("✓ Results tables saved to:", res_dir, "\n")
