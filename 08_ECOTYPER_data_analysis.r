#Incomplete for TEMPUS assignment- ran out of time

#ECOTYPER data analysis
# ============================================================
# EcoTyper output parsing — Carcinoma Cell States
# Correlate cell state abundances with NMF cluster membership
# ============================================================

#read in VST matrix if necessary
vst <- readRDS("vst_her2_brca.rds")
vst_mat <- assay(vst)

#read in NMF cluster assignments
nmf_df <- read.csv("clustering_results/nmf_cluster_assignments.csv",
                    stringsAsFactors = FALSE)

#Generate outputs for ECOTYPER
#Download LM22 signature matrix from ECOTYPER website and read in gene list
lm22 <- read.table("LM22.txt", sep = "\t", header = TRUE)
lm22_genes <- lm22$Gene.symbol

# Filter VST matrix to LM22 genes only to reduce file size for ECOTYPER input
lm22_genes_present <- intersect(lm22_genes, rownames(vst_mat))

cat(sprintf("LM22 genes found in matrix: %d / %d\n", 
            length(lm22_genes_present), length(lm22_genes)))

mixture <- vst_mat[lm22_genes_present, nmf_df$patient_id] %>%
  as.data.frame() %>%
  rownames_to_column("Gene")

write.table(mixture, "cibersortx_mixture.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


library(tidyverse)
library(pheatmap)
library(ggplot2)

ecotyper_base <- "ecotyper_output/Carcinoma_Cell_States"
out_dir       <- "ecotyper_nmf"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Cell types to parse ───────────────────────────────────────
# Adjust this list to match exactly the folders present
cell_types <- c("B.cells", "CD4.T.cells", "CD8.T.cells", "Dendritic.cells",
                "Endothelial.cells", "Fibroblasts", "Monocytes.and.Macrophages",
                "Mast.cells", "NK.cells", "PCs", "PMNs", 
                "Epithelial.cells")

# ── Read and combine all state abundances ─────────────────────
abundance_list <- lapply(cell_types, function(ct) {
  f <- file.path(ecotyper_base, ct, paste0(ct, "_Cell_State_Abundance.txt"))
  if (!file.exists(f)) {
    message(sprintf("Skipping %s — file not found", ct))
    return(NULL)
  }
  df <- read.table(f, sep = "\t", header = TRUE, check.names = FALSE) %>%
    pivot_longer(-ID, names_to = "state", values_to = "abundance") %>%
    mutate(cell_type = ct,
           state_full = paste0(ct, "_", state))
  df
})

abundance_df <- bind_rows(abundance_list) %>%
  rename(ID = "patient_id") %>%
  left_join(nmf_df %>% dplyr::select("patient_id", "nmf_cluster"),
            by = "patient_id")

# ── Read and combine all state assignments ────────────────────
assignment_list <- lapply(cell_types, function(ct) {
  f <- file.path(ecotyper_base, ct, paste0(ct, "_Cell_State_Assignment.txt"))
  if (!file.exists(f)) return(NULL)
  read.table(f, sep = "\t", header = TRUE, check.names = FALSE) %>%
    mutate(cell_type = ct) %>%
    rename(ID = "patient_id")
})

assignment_df <- bind_rows(assignment_list) %>%
  left_join(nmf_df %>% dplyr::select(patient_id, nmf_cluster),
            by = "patient_id") %>%
  filter(!is.na(nmf_cluster))

# ── PLOT 1: Boxplot of state abundances by NMF cluster ────────
# One plot per cell type
for (ct in cell_types) {
  df_ct <- abundance_df %>% filter(cell_type == ct)
  if (nrow(df_ct) == 0) next

  p <- ggplot(df_ct, aes(x = nmf_cluster, y = abundance, fill = nmf_cluster)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~state, scales = "free_y") +
    theme_bw() +
    labs(title = paste("EcoTyper —", gsub("\\.", " ", ct)),
         x = "NMF Cluster", y = "State abundance") +
    theme(legend.position = "none")

  ggsave(file.path(out_dir, paste0("ecotyper_", ct, "_by_nmf.pdf")),
         p, width = 8, height = 6)
}

# ── PLOT 2: Heatmap — mean abundance per state per NMF cluster ─
# Rows = cell states, Columns = NMF clusters
mean_abundance <- abundance_df %>%
  group_by(nmf_cluster, state_full) %>%
  summarise(mean_abundance = mean(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = nmf_cluster, values_from = mean_abundance) %>%
  column_to_rownames("state_full")

# Remove states with zero variance across clusters
row_vars <- apply(mean_abundance, 1, var, na.rm = TRUE)
mean_abundance <- mean_abundance[row_vars > 0, ]

pheatmap(
  as.matrix(mean_abundance),
  scale        = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color        = colorRampPalette(c("steelblue", "white", "tomato"))(100),
  main         = "EcoTyper cell state abundance by NMF cluster",
  filename     = file.path(out_dir, "ecotyper_states_heatmap.pdf"),
  width        = 6,
  height       = max(6, nrow(mean_abundance) * 0.25)
)

# ── PLOT 3: State assignment frequencies by NMF cluster ───────
# What proportion of each NMF cluster is assigned to each state?
assignment_freq <- assignment_df %>%
  count(nmf_cluster, cell_type, State) %>%
  group_by(nmf_cluster, cell_type) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

for (ct in cell_types) {
  df_ct <- assignment_freq %>% filter(cell_type == ct)
  if (nrow(df_ct) == 0) next

  p <- ggplot(df_ct, aes(x = nmf_cluster, y = prop, fill = State)) +
    geom_col(position = "stack") +
    theme_bw() +
    labs(title = paste("State assignments —", gsub("\\.", " ", ct)),
         x = "NMF Cluster", y = "Proportion of samples",
         fill = "Cell state") +
    scale_y_continuous(labels = scales::percent)

  ggsave(file.path(out_dir, paste0("ecotyper_", ct, "_assignment_stack.pdf")),
         p, width = 5, height = 4)
}

# ── Save combined data ────────────────────────────────────────
saveRDS(list(
  abundance   = abundance_df,
  assignment  = assignment_df,
  mean_by_cluster = mean_abundance
), file.path(out_dir, "ecotyper_parsed.rds"))

cat("Done. Outputs saved to:", out_dir, "\n")