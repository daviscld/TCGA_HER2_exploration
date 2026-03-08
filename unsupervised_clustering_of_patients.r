#HER2+ patient groups

# ============================================================
# Unsupervised Clustering: UMAP + NMF
# Input: VST SummarizedExperiment + metadata with clinical HER2
# ============================================================

library(DESeq2); library(SummarizedExperiment)
library(dplyr); library(tibble); library(tidyr); library(stringr)
library(ggplot2); library(patchwork)
library(umap)             # UMAP
library(dbscan)           # HDBSCAN clustering on UMAP
library(NMF)              # NMF with rank survey
library(ComplexHeatmap); library(circlize)
library(clusterProfiler); library(org.Hs.eg.db); library(msigdbr)
library(DESeq2)           # re-used for within-cluster DE
library(EnhancedVolcano)
library(pheatmap)
library(XVector)

out_dir <- "clustering_plots";   dir.create(out_dir, showWarnings = FALSE)
res_dir <- "clustering_results"; dir.create(res_dir, showWarnings = FALSE)

random_seed  <- 42
n_hvg        <- 3000     # highly variable genes for UMAP + NMF input
min_cluster  <- 10       # minimum samples to retain a cluster
top_de_genes <- 50       # genes shown in DE heatmap
palette_her2 <- c("HER2+" = "#E64B35", "HER2-" = "#4DBBD5", 
                "Equivocal" = "#F39B7F")

save_plot <- function(p, fname, w=8, h=6, dpi=300){
  ggsave(file.path(out_dir, fname), p, width=w, height=h, dpi=dpi, bg="white")}

# ═══════════════════════════════════════════════════════════
# 1. SUBSET TO CLINICALLY ASSIGNED PATIENTS
# ═══════════════════════════════════════════════════════════
# Requires mol_df from multimodal definition script with
# clinical_label column ("HER2+" / "HER2-" / "Equivocal")

meta    <- read.csv("metadata_qc_passed.csv", row.names = 1)
vst     <- readRDS("vst_her2_brca.rds")
vst_mat <- assay(vst)

# rebuild clinical_label from meta

meta_ids <- if ("patient_id" %in% colnames(meta)) meta$patient_id else rownames(meta)
meta$FISH_HER2 <- meta[meta_ids, "HER2_fish_status"]
mol_df <- meta %>%
  mutate(
    clinical_label = case_when(
      FISH_HER2 == "Positive"                           ~ "HER2+",
      FISH_HER2 == "Negative"                           ~ "HER2-",
      is.na(FISH_HER2) & IHC_HER2 == "Positive"        ~ "HER2+",
      is.na(FISH_HER2) & IHC_HER2 == "Negative"        ~ "HER2-",
      is.na(FISH_HER2) & IHC_HER2 == "Equivocal"       ~ "Equivocal",
      TRUE                                              ~ NA_character_
    ),
    clinical_label = factor(
      clinical_label,
      levels = c("HER2+", "HER2-", "Equivocal")
    )
  )


#if ER and PR status are available, rename them to similar to HER2
if ("PR_status_by_ihc" %in% colnames(meta)) {
  mol_df <- mol_df %>%
    mutate(PR_label = case_when(
      PR_status_by_ihc == "Positive" ~ "PR+",
      PR_status_by_ihc == "Negative" ~ "PR-",
      TRUE                          ~ NA_character_
    ))
}
if ("ER_Status_By_IHC" %in% colnames(meta)) {
  mol_df <- mol_df %>%
    mutate(ER_label = case_when(
      ER_Status_By_IHC == "Positive" ~ "ER+",
      ER_Status_By_IHC == "Negative" ~ "ER-",
      TRUE                          ~ NA_character_
    ))
}

# Keep only samples with a definitive clinical label (drop Equivocal + NA)
assigned <- mol_df %>%
  filter(clinical_label %in% c("HER2+", "HER2-")) %>%
  filter(patient_id %in% colnames(vst_mat))

vst_sub <- vst_mat[, assigned$patient_id]
cat(sprintf("Clinically assigned samples: %d  (HER2+: %d  HER2-: %d)\n",
            nrow(assigned), sum(assigned$clinical_label == "HER2+"),
            sum(assigned$clinical_label == "HER2-")))

# ═══════════════════════════════════════════════════════════
# 2. HIGHLY VARIABLE GENE SELECTION
# ═══════════════════════════════════════════════════════════
gene_var <- apply(vst_sub, 1, var)
hvg      <- names(sort(gene_var, decreasing = TRUE))[1:n_hvg]
#check for ERBB2 and STARD3 in the selected HVGs
if (!"ERBB2" %in% hvg) warning("ERBB2 not in top HVGs — may not be biological")
if (!"STARD3" %in% hvg) warning("STARD3 not in top HVGs — may not be biological")
vst_hvg  <- vst_sub[hvg, ]   # HVG × samples

cat(sprintf("Using top %d highly variable genes\n", n_hvg))

# ═══════════════════════════════════════════════════════════
# 3. UMAP + CLUSTERING
# ═══════════════════════════════════════════════════════════

# ── 3A. PCA pre-reduction (standard before UMAP on bulk RNA) ─
set.seed(random_seed)
pca_res  <- prcomp(t(vst_hvg), scale. = TRUE)
n_pcs    <- min(50, ncol(pca_res$x))   # use up to 50 PCs as UMAP input
pca_scores <- pca_res$x[, 1:n_pcs]

# Scree to justify n_pcs choice
pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)
cum_var <- cumsum(pct_var)
n_80    <- which(cum_var >= 80)[1]
cat(sprintf("PCs explaining ≥80%% variance: %d (using %d for UMAP)\n", n_80, n_pcs))

# ── 3B. UMAP ──────────────────────────────────────────────
umap_config <- umap.defaults
umap_config$random_state  <- random_seed
umap_config$n_neighbors   <- 15
umap_config$min_dist      <- 0.1
umap_config$metric        <- "euclidean"

umap_res <- umap(pca_scores, config = umap_config)
umap_df  <- as.data.frame(umap_res$layout) %>%
  setNames(c("UMAP1","UMAP2")) %>%
  rownames_to_column("patient_id") %>%
  left_join(assigned %>% 
  dplyr::select(patient_id, clinical_label, ER_label, PR_label), 
  by = "patient_id")

# ── 3C. HDBSCAN clustering on UMAP coordinates ────────────
# minPts controls minimum cluster size — tune if clusters are too granular
set.seed(random_seed)
hdb <- hdbscan(umap_df[, c("UMAP1","UMAP2")], minPts = 15)
umap_df$umap_cluster <- factor(hdb$cluster)   # 0 = noise points

# Drop noise (cluster 0) and tiny clusters
cluster_sizes <- table(umap_df$umap_cluster)
keep_clusters <- names(cluster_sizes)[cluster_sizes >= min_cluster & names(cluster_sizes) != "0"]
umap_df <- umap_df %>%
  mutate(umap_cluster = ifelse(umap_cluster %in% keep_clusters,
                               paste0("C", umap_cluster), "Noise"))

cat("\nUMAP cluster sizes:\n"); print(table(umap_df$umap_cluster))

# ── 3D. HER2 enrichment per cluster ───────────────────────
her2_enrich <- umap_df %>%
  filter(umap_cluster != "Noise") %>%
  count(umap_cluster, clinical_label) %>%
  group_by(umap_cluster) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

# Fisher's exact test: HER2+ enrichment vs rest for each cluster
her2_fisher <- lapply(unique(umap_df$umap_cluster[umap_df$umap_cluster != "Noise"]),
                      function(cl) {
  tab <- table(
    in_cluster = umap_df$umap_cluster == cl,
    is_her2pos = umap_df$clinical_label == "HER2+"
  )
  ft  <- fisher.test(tab)
  data.frame(cluster = cl, OR = ft$estimate, p = ft$p.value)
}) %>% bind_rows() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

cat("\nHER2+ enrichment per UMAP cluster (Fisher):\n")
print(her2_fisher)
write.csv(her2_fisher, file.path(res_dir, "UMAP_HER2_enrichment.csv"), row.names = FALSE)

# ── 3E. UMAP plots ────────────────────────────────────────
p_umap_her2 <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = clinical_label)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = palette_her2) +
  labs(title = "UMAP — coloured by clinical HER2 label", colour = NULL) +
  theme_bw(base_size = 12)

#make PR and ER palette
palette_pr <- c("PR+" = "#E64B35", "PR-" = "#4DBBD5", "Equivocal" = "#F39B7F")
palette_er <- c("ER+" = "#E64B35", "ER-" = "#4DBBD5", "Equivocal" = "#F39B7F")

p_umap_PR <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = as.factor(PR_label))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = palette_pr) +
  labs(title = "UMAP — coloured by clinical PR label", colour = NULL) +
  theme_bw(base_size = 12)

p_umap_ER <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = as.factor(ER_label))) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_manual(values = palette_er) +
  labs(title = "UMAP — coloured by clinical ER label", colour = NULL) +
  theme_bw(base_size = 12)

p_umap_clust <- ggplot(umap_df, aes(UMAP1, UMAP2, colour = umap_cluster)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "UMAP — HDBSCAN clusters", colour = "Cluster") +
  theme_bw(base_size = 12)

save_plot(p_umap_her2 + p_umap_PR + p_umap_ER + p_umap_clust, 
    "01_UMAP_her2_pr_er_and_clusters.png", w=14, h=6)

# HER2+ fraction per cluster bar chart
p_her2_frac <- her2_enrich %>%
  filter(clinical_label == "HER2+") %>%
  left_join(her2_fisher %>% dplyr::select(cluster, p_adj), by = c("umap_cluster"="cluster")) %>%
  ggplot(aes(x = reorder(umap_cluster, -pct), y = pct,
             fill = p_adj < 0.05)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 100 * mean(assigned$clinical_label == "HER2+"),
             linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = c("TRUE"="#E64B35","FALSE"="grey70"),
                    labels = c("TRUE"="FDR<0.05","FALSE"="NS")) +
  labs(title = "HER2+ Fraction per UMAP Cluster",
       subtitle = "Dashed line = cohort prevalence",
       x = "Cluster", y = "% HER2+", fill = NULL) +
  theme_bw(base_size = 12)
save_plot(p_her2_frac, "02_HER2_fraction_per_cluster.png", w=7, h=5)

# ── 3F. Metadata enrichment per cluster ───────────────────
# Continuous: violin plots; categorical: stacked bar charts
meta_sub <- meta[assigned$patient_id, ] %>%
  left_join(umap_df %>% dplyr::select(patient_id, umap_cluster), by = "patient_id") %>%
  filter(umap_cluster != "Noise")

# Identify continuous and categorical metadata columns of interest
# Adjust these column names to match your metadata
continuous_cols    <- intersect(c("Fraction_Genome_Altered","Mutation_Count",
                                   "Diagnosis_Age", "Overall_Survival_Months"),
                                colnames(meta_sub))
categorical_cols   <- intersect(c("PR_status_by_ihc","ER_Status_By_IHC",
                                   "Race_Category","Menopause_Status",
                                   "Cancer_Stage"), colnames(meta_sub))

for (col in continuous_cols) {
  p <- ggplot(meta_sub %>% filter(!is.na(.data[[col]])),
              aes(x = umap_cluster, y = .data[[col]], fill = umap_cluster)) +
    geom_violin(alpha = 0.6) + geom_boxplot(width = 0.12, outlier.size = 0.4) +
    stat_compare_means(method = "wilcox.test", 
            label.y = max(meta_sub[[col]], na.rm = TRUE) * 1.05) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = paste(col, "by UMAP cluster"), x = NULL, y = col) +
    theme_bw(base_size = 11) + theme(legend.position = "none")
  save_plot(p, paste0("03_meta_", make.names(col), "_by_cluster.png"), w=7, h=4)
}

for (col in categorical_cols) {
  p <- meta_sub %>%
    filter(!is.na(.data[[col]])) %>%
    count(umap_cluster, .data[[col]]) %>%
    group_by(umap_cluster) %>% mutate(pct = 100 * n / sum(n)) %>%
    ggplot(aes(x = umap_cluster, y = pct, fill = .data[[col]])) +
    geom_col(width = 0.7) +
    stat_compare_means(method = "chisq.test", label.y = 75) +
    labs(title = paste(col, "composition by UMAP cluster"),
         x = NULL, y = "%", fill = col) +
    theme_bw(base_size = 11)
  save_plot(p, paste0("04_meta_", make.names(col), "_by_cluster.png"), w=8, h=5)
}

# ═══════════════════════════════════════════════════════════
# 4. DE + GSEA
# ═══════════════════════════════════════════════════════════
#define a function to run GSEA for clusters
run_de_gsea <- function(vst_counts, sample_meta, cluster_col,
                        prefix, top_n = top_de_genes) {

  clusters     <- sort(unique(sample_meta[[cluster_col]]))
  clusters     <- clusters[clusters != "Noise"]
  cluster_pairs <- combn(clusters, 2, simplify = FALSE)

  all_de   <- list()
  all_gsea <- list()

  for (pair in cluster_pairs) {
    cl_a <- pair[1]; cl_b <- pair[2]
    tag  <- paste0(cl_a, "_vs_", cl_b)

    samp_ab <- sample_meta %>%
      filter(.data[[cluster_col]] %in% c(cl_a, cl_b))

    counts_ab <- round(vst_counts[, samp_ab$patient_id])  # DESeq2 needs integers
    col_ab    <- samp_ab %>% column_to_rownames("patient_id")
    col_ab[[cluster_col]] <- factor(col_ab[[cluster_col]], levels = c(cl_b, cl_a))

    dds <- DESeqDataSetFromMatrix(counts_ab, col_ab,
                                  design = as.formula(paste("~", cluster_col)))
    dds <- DESeq(dds, quiet = TRUE)
    res <- lfcShrink(dds, coef = 2, type = "ashr") %>%
      as.data.frame() %>% rownames_to_column("gene") %>%
      arrange(padj) %>% filter(!is.na(padj))

    all_de[[tag]] <- res
    write.csv(res, file.path(res_dir, paste0(prefix, "_DE_", tag, ".csv")),
              row.names = FALSE)

    # Heatmap of top DE genes
    top_genes <- res %>%
      filter(padj < 0.05) %>%
      slice_max(abs(log2FoldChange), n = top_n) %>%
      pull(gene)

    if (length(top_genes) >= 5) {
      mat   <- vst_counts[top_genes, samp_ab$patient_id]
      mat_z <- t(scale(t(mat)))
      mat_z[mat_z >  3] <-  3; mat_z[mat_z < -3] <- -3

      ha <- HeatmapAnnotation(
        Cluster = col_ab[[cluster_col]],
        ER_label = sample_meta$ER_label[match(samp_ab$patient_id, sample_meta$patient_id)],
        PR_label = sample_meta$PR_label[match(samp_ab$patient_id, sample_meta$patient_id)],
        HER2    = sample_meta$clinical_label[match(samp_ab$patient_id, sample_meta$patient_id)],
        col = list(
          Cluster = setNames(RColorBrewer::brewer.pal(8,"Set2")[seq_along(clusters)], clusters),
            ER_label = palette_er,
            PR_label = palette_pr,
          HER2    = palette_her2
        )
      )
      ht <- Heatmap(mat_z, name = "Z-score",
        col               = colorRamp2(c(-3,0,3), c("#4DBBD5","white","#E64B35")),
        top_annotation    = ha,
        show_column_names = FALSE,
        column_split      = col_ab[[cluster_col]],
        cluster_column_slices = FALSE,
        row_names_gp      = gpar(fontsize = 7),
        column_title      = paste(prefix, tag),
        heatmap_legend_param = list(title = "Z-score")
      )
      png(file.path(out_dir, paste0(prefix, "_heatmap_", tag, ".png")),
          width=1000, height=900, res=130)
      draw(ht); dev.off()
    }

    # GSEA with Hallmark gene sets
    hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
      dplyr::select(gs_name, gene_symbol)

    rank_vec <- res %>%
      filter(!is.na(log2FoldChange), !is.na(padj)) %>%
      mutate(stat = log2FoldChange * -log10(padj + 1e-300)) %>%
      arrange(desc(stat)) %>%
      dplyr::select(gene, stat) %>%
      distinct(gene, .keep_all = TRUE) %>%
      deframe()

    gsea_res <- tryCatch(
      GSEA(rank_vec, TERM2GENE = hallmark, pvalueCutoff = 0.05,
           minGSSize = 15, maxGSSize = 500, seed = random_seed, verbose = FALSE),
      error = function(e) NULL
    )

    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      all_gsea[[tag]] <- gsea_res
      write.csv(as.data.frame(gsea_res),
                file.path(res_dir, paste0(prefix, "_GSEA_", tag, ".csv")),
                row.names = FALSE)

      p_gsea <- dotplot(gsea_res, showCategory = 20, split = ".sign",
                        font.size = 8) +
        facet_grid(~.sign) +
        labs(title = paste(prefix, "GSEA Hallmark —", tag)) +
        theme_bw(base_size = 9)
      save_plot(p_gsea, paste0(prefix, "_GSEA_", tag, ".png"), w=12, h=8)
    }
  }

  list(de = all_de, gsea = all_gsea)
}

# Reload raw counts for DESeq2 (VST is not suitable for DE testing)
#this will allow for DE testing between UMAP clusters
dds_full <- readRDS("dds_her2_brca.rds")
raw_counts <- counts(dds_full)[, assigned$patient_id]

umap_meta <- umap_df %>% filter(umap_cluster != "Noise")
umap_de_gsea <- run_de_gsea(raw_counts, umap_meta, 
                    "umap_cluster", "UMAP")

# ═══════════════════════════════════════════════════════════
# 5. NMF
# ═══════════════════════════════════════════════════════════

# ── 5A. Prepare NMF input ─────────────────────────────────
# NMF requires non-negative input — shift VST to be >= 0
nmf_mat <- vst_hvg - min(vst_hvg)   # genes × samples, non-negative

# ── 5B. Rank survey (k = 2 to 8) ─────────────────────────
# Cophenetic correlation and dispersion determine the best rank.
# This is the most computationally expensive step — reduce nrun if slow.
cat("\nRunning NMF rank survey (k=2:8, nrun=30)...\n")
set.seed(random_seed)
nmf_survey <- nmf(nmf_mat, rank = 2:8, method = "brunet",
                  nrun = 30, seed = random_seed, .opt = "vp4")
# v = verbose, p4 = 4 cores; adjust p to your core count or remove for single-core

# Plot rank selection metrics
png(file.path(out_dir, "05_NMF_rank_survey.png"), 
        width=4000, height=3000, res=300)
plot(nmf_survey)
dev.off()

# Cophenetic + dispersion summary
coph <- summary(nmf_survey)["cophenetic"]
disp <- summary(nmf_survey)["dispersion"]
rank_df <- data.frame(k = seq(2,8,1), cophenetic = coph, dispersion = disp)
cat("\nNMF rank survey:\n"); print(rank_df)
write.csv(rank_df, file.path(res_dir, "NMF_rank_survey.csv"), row.names = FALSE)

# ── 5C. Select optimal rank ───────────────────────────────
# Choose k where cophenetic is high and starts to decrease,
# and dispersion is high. The "elbow" in cophenetic is the key signal.

#chosen by looking for drop in dispersion, cophenetic, 
#silhouette score
best_k <- 4 

# ── 5D. Final NMF run at best_k ───────────────────────────
cat(sprintf("\nRunning final NMF at k=%d (nrun=100)...\n", best_k))
set.seed(random_seed)
nmf_final <- nmf(nmf_mat, rank = best_k, method = "brunet",
                 nrun = 100, seed = random_seed, .opt = "vp4")
saveRDS(nmf_final, file.path(res_dir, sprintf("NMF_final_k%d.rds", best_k)))

# Assign each sample to its dominant NMF factor (basis component)
nmf_clusters <- predict(nmf_final, what = "samples")   # named factor vector
nmf_df <- data.frame(
  patient_id   = names(nmf_clusters),
  nmf_cluster = paste0("NMF", as.integer(nmf_clusters))) %>% 
            left_join(assigned %>% 
        dplyr::select(patient_id, clinical_label, ER_label, PR_label), 
        by = "patient_id")

cat("\nNMF cluster sizes:\n"); print(table(nmf_df$nmf_cluster))

# HER2 enrichment in NMF clusters
nmf_fisher <- lapply(unique(nmf_df$nmf_cluster), function(cl) {
    tab <- table(in_cluster = nmf_df$nmf_cluster == cl,
               is_her2pos = nmf_df$clinical_label == "HER2+")
    ft  <- fisher.test(tab)
  data.frame(cluster = cl, OR = ft$estimate, p = ft$p.value)
}) 
%>% bind_rows() %>% mutate(p_adj = p.adjust(p, method = "BH"))

cat("\nHER2+ enrichment per NMF cluster (Fisher):\n"); print(nmf_fisher)
write.csv(nmf_fisher, file.path(res_dir, "NMF_HER2_enrichment.csv"), row.names = FALSE)

# ── 5E. NMF metagene heatmap ──────────────────────────────
# Top genes per NMF factor (high loading = characteristic of that program)
w_mat     <- basis(nmf_final)   # genes × k
top_genes_nmf <- apply(w_mat, 2, function(w)
  rownames(w_mat)[order(w, decreasing=TRUE)[1:30]]
)
all_nmf_genes <- unique(as.vector(top_genes_nmf))

ha_nmf <- HeatmapAnnotation(
  NMF_cluster = nmf_df$nmf_cluster,
  HER2        = nmf_df$clinical_label,
  PR_label    = nmf_df$PR_label,
  ER_label    = nmf_df$ER_label,
  col = list(
    NMF_cluster = setNames(
      RColorBrewer::brewer.pal(max(best_k, 3), "Set1")[1:best_k],
      paste0("NMF", 1:best_k)),
    HER2 = palette_her2,
    PR_label = palette_pr,
    ER_label = palette_er
  )
)

nmf_heat_mat <- vst_hvg[all_nmf_genes, nmf_df$patient_id]
nmf_heat_z   <- t(scale(t(nmf_heat_mat)))
nmf_heat_z[nmf_heat_z > 3] <- 3; nmf_heat_z[nmf_heat_z < -3] <- -3

ht_nmf <- Heatmap(nmf_heat_z, name = "Z-score",
  col               = colorRamp2(c(-3,0,3), 
        c("#4DBBD5","white","#E64B35")),
  top_annotation    = ha_nmf,
  show_column_names = FALSE,
  column_split      = nmf_df$nmf_cluster,
  cluster_column_slices = FALSE,
  row_names_gp      = gpar(fontsize = 7),
  column_title_gp   = gpar(fontsize = 10),
  heatmap_legend_param = list(title = "Z-score")
)
png(file.path(out_dir, "06_NMF_metagene_heatmap.png"),
    width=1100, height=1000, res=300)
draw(ht_nmf); dev.off()
#also make a pdf for use in Illustrator to make final figure panels
pdf(file.path(out_dir, "06_NMF_metagene_heatmap.pdf"),
    width=11, height=10)
draw(ht_nmf); dev.off()

# ── 5F. DE + GSEA between NMF clusters ────────────────────
nmf_de_gsea <- run_de_gsea(raw_counts, nmf_df, "nmf_cluster", "NMF")
saveRDS(list(umap_de_gsea = umap_de_gsea, nmf_de_gsea = nmf_de_gsea),
        file.path(res_dir, "DE_GSEA_results.rds"))

# ── 5G. NMF on UMAP ───────────────────────────────────────
umap_nmf_df <- umap_df %>%
  left_join(nmf_df %>% dplyr::select(patient_id, nmf_cluster), by = "patient_id")

p_umap_nmf <- ggplot(umap_nmf_df, aes(UMAP1, UMAP2, colour = nmf_cluster)) +
  geom_point(size = 1.5, alpha = 0.7) +
  scale_colour_brewer(palette = "Set1") +
  labs(title = "UMAP — coloured by NMF cluster", colour = "NMF") +
  theme_bw(base_size = 12)
save_plot(p_umap_nmf, "07_UMAP_NMF_overlay.png", w=7, h=6)

# Plot GSEA results for NMF clusters
plt <- lapply(names(nmf_de_gsea$gsea), function(tag) {
  gsea_res <- nmf_de_gsea$gsea[[tag]]
  if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
    p <- dotplot(gsea_res, showCategory = 20, split = ".sign",
                 font.size = 8) +
      facet_grid(~.sign) +
      labs(title = paste("NMF GSEA Hallmark —", tag)) +
      theme_bw(base_size = 9)
    save_plot(p, paste0("08_NMF_GSEA_", tag, ".png"), w=12, h=8)
    return(p)
  } else {
    return(NULL)
  }
}) %>% compact() %>% wrap_plots(ncol = 2)
save_plot(plt, "08_NMF_GSEA_all.png", w=14, h=14)
#save as pdf for figure creation
pdf(file.path(out_dir, "08_NMF_GSEA_all.pdf"), width=14, height=14)
print(plt); dev.off()