# ============================================================
# sPLS-DA + Random Forest to classify NMF clusters
# Requires:
#   nmf_df       – data.frame with patient_id, nmf_cluster
#   vst_mat      – VST expression matrix (genes × samples)
#   nmf_de_gsea  – DE results list with $de named list of DESeq2 result dfs
# ============================================================

library(mixOmics)
library(ranger)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(caret)       # for confusionMatrix

# ============================================================
# SHARED SETUP
# ============================================================

top_degs_per_cluster <- 100
padj_thresh          <- 0.05
lfc_thresh           <- 0.58
out_dir              <- "splsda_rf_nmf"
seed                 <- 42

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ── Load data ────────────────────────
nmf_results <- readRDS("clustering_results/NMF_final_k4.rds")
nmf_df <- read.csv("clustering_results/nmf_cluster_assignments.csv", 
                    stringsAsFactors = FALSE)
vst     <- readRDS("vst_her2_brca.rds")
vst_mat <- assay(vst)
vst_mat <- vst_mat[, nmf_df$patient_id]  # align columns to nmf_df order
nmf_de_gsea <- readRDS("clustering_results/DE_GSEA_results.rds")
#take only nmf results
nmf_de_gsea <- nmf_de_gsea[["nmf_de_gsea"]]

# ── Collect top DEGs  ────────────────────────
top_deg_genes <- lapply(nmf_de_gsea$de, function(res_df) {
  res_df %>%
    filter(padj < padj_thresh, abs(log2FoldChange) >= lfc_thresh) %>%
    slice_max(abs(log2FoldChange), n = top_degs_per_cluster, with_ties = FALSE) %>%
    pull(gene)
}) %>% unlist() %>% unique()

top_deg_genes <- intersect(top_deg_genes, rownames(vst_mat))
cat(sprintf("DEG feature set: %d genes\n", length(top_deg_genes)))

# ── Build X and y ─────────────────────────────────────────────
nmf_samples <- nmf_df %>% filter(!is.na(nmf_cluster)) %>% pull(patient_id)

X_full <- t(vst_mat[top_deg_genes, nmf_samples]) %>% as.data.frame()

y_full <- nmf_df %>%
  filter(!is.na(nmf_cluster)) %>%
  dplyr::select(patient_id, nmf_cluster) %>%
  deframe()

y_full <- as.factor(y_full[rownames(X_full)])

# ── Stratified holdout (15%) ───────────────────────────────────
set.seed(seed)
holdout_idx <- unlist(lapply(levels(y_full), function(cls) {
  idx <- which(y_full == cls)
  sample(idx, max(1, floor(length(idx) * 0.15)))
}))
train_idx <- setdiff(seq_len(nrow(X_full)), holdout_idx)

X_train <- X_full[train_idx, , drop = FALSE]
y_train <- y_full[train_idx]
X_hold  <- X_full[holdout_idx, , drop = FALSE]
y_hold  <- y_full[holdout_idx]

# ── Downsample training set to balance classes ────────────────
set.seed(seed)
min_n <- min(table(y_train))
ds_idx <- unlist(lapply(levels(y_train), function(cls) {
  idx <- which(y_train == cls)
  sample(idx, min_n)
}))
X_train <- X_train[ds_idx, , drop = FALSE]
y_train <- y_train[ds_idx]
cat(sprintf("After downsampling: %d training samples\n", nrow(X_train)))

# ── Scale (store params for holdout) ──────────────────────────
X_scaled    <- scale(X_train, center = TRUE, scale = TRUE)
scale_center <- attr(X_scaled, "scaled:center")
scale_sd     <- attr(X_scaled, "scaled:scale")

X_hold_scaled <- scale(X_hold,
                        center = scale_center,
                        scale  = scale_sd)

# ============================================================
# PART 1: SPARSE PLS-DA (mixOmics)
# ============================================================
cat("\n── Running sPLS-DA ─────────────────────────────────────\n")

# ── Tune number of components and keepX ───────────────────────
# keepX controls how many genes are retained per component
# Test a range and pick by BER (balanced error rate)

set.seed(seed)
tune_splsda <- mixOmics::tune.splsda(
  X          = X_scaled,
  Y          = y_train,
  ncomp      = 3,
  validation = "Mfold",
  folds      = 5,
  nrepeat    = 10,
  dist       = "max.dist",
  test.keepX = c(5, 10, 15, 20, 30, 50),
  progressBar = TRUE
)

optimal_keepX <- tune_splsda$choice.keepX
optimal_ncomp <- tune_splsda$choice.ncomp$ncomp
cat(sprintf("Optimal ncomp: %d, keepX: %s\n",
            optimal_ncomp,
            paste(optimal_keepX, collapse = ", ")))

# ── Fit final sPLS-DA model ───────────────────────────────────
splsda_model <- mixOmics::splsda(
  X      = X_scaled,
  Y      = y_train,
  ncomp  = optimal_ncomp,
  keepX  = optimal_keepX
)

# ── Cross-validation performance ─────────────────────────────
set.seed(seed)
splsda_perf <- mixOmics::perf(
  splsda_model,
  validation  = "Mfold",
  folds       = 5,
  nrepeat     = 10,
  dist        = "max.dist",
  progressBar = TRUE
)

cat("\nsPLS-DA CV Error Rates:\n")
print(splsda_perf$error.rate)

# ── Holdout prediction ────────────────────────────────────────
splsda_hold_pred <- predict(splsda_model,
                             newdata = X_hold_scaled,
                             dist    = "max.dist")
splsda_hold_class <- splsda_hold_pred$class$max.dist[, optimal_ncomp]
splsda_hold_acc   <- mean(splsda_hold_class == y_hold)
cat(sprintf("\nsPLS-DA holdout accuracy: %.3f\n", splsda_hold_acc))

print(caret::confusionMatrix(as.factor(splsda_hold_class), y_hold))

# ── Extract selected genes per component ──────────────────────
splsda_genes <- lapply(seq_len(optimal_ncomp), function(comp) {
  mixOmics::selectVar(splsda_model, comp = comp)$name
})
names(splsda_genes) <- paste0("Comp", seq_len(optimal_ncomp))

all_splsda_genes <- unique(unlist(splsda_genes))
cat(sprintf("\nsPLS-DA selected %d unique genes across %d components\n",
            length(all_splsda_genes), optimal_ncomp))

# ── Plots ─────────────────────────────────────────────────────

# Score plot comp 1 vs 2
pdf(file.path(out_dir, "splsda_scores_comp1v2.pdf"), width = 6, height = 5)
mixOmics::plotIndiv(splsda_model,
                    comp       = c(1, 2),
                    group      = y_train,
                    ind.names  = FALSE,
                    ellipse    = TRUE,
                    legend     = TRUE,
                    title      = "sPLS-DA — Comp 1 vs 2")
dev.off()

# Score plot comp 1 vs 3 (if ncomp >= 3)
if (optimal_ncomp >= 3) {
  pdf(file.path(out_dir, "splsda_scores_comp1v3.pdf"), width = 6, height = 5)
  mixOmics::plotIndiv(splsda_model,
                      comp       = c(1, 3),
                      group      = y_train,
                      ind.names  = FALSE,
                      ellipse    = TRUE,
                      legend     = TRUE,
                      title      = "sPLS-DA — Comp 1 vs 3")
  dev.off()
}

# Loading plot per component
for (comp in seq_len(optimal_ncomp)) {
  pdf(file.path(out_dir, sprintf("splsda_loadings_comp%d.pdf", comp)),
      width = 5, height = 4)
  mixOmics::plotLoadings(splsda_model,
                          comp   = comp,
                          method = "mean",
                          contrib = "max",
                          title  = sprintf("sPLS-DA Loadings — Comp %d", comp))
  dev.off()
}

# Heatmap of sPLS-DA components per sample
splsda_scores <- as.data.frame(splsda_model$variates$X)
colnames(splsda_scores) <- paste0("Comp", seq_len(ncol(splsda_scores)))

anno_row <- data.frame(NMF_cluster = y_train,
                        row.names   = rownames(splsda_scores))

pheatmap(
  as.matrix(splsda_scores),
  annotation_row = anno_row,
  cluster_cols   = FALSE,
  cluster_rows   = FALSE,
  show_rownames  = FALSE,
  scale          = "column",
  color          = colorRampPalette(c("steelblue", "white", "tomato"))(100),
  main           = "sPLS-DA component scores per sample",
  filename       = file.path(out_dir, "splsda_scores_heatmap.pdf"),
  width          = 6,
  height         = 7
)

# ============================================================
# PART 2: RANDOM FOREST (ranger)
# ============================================================
cat("\n── Running Random Forest ───────────────────────────────\n")

# ── Fit RF on training set ────────────────────────────────────
set.seed(seed)
rf_model <- ranger::ranger(
  x              = X_scaled,
  y              = y_train,
  num.trees      = 1000,
  importance     = "permutation",   # permutation importance is more reliable
  probability    = TRUE,            # output class probabilities
  num.threads    = parallel::detectCores() - 1,
  seed           = seed
)

# ── CV accuracy via OOB error ─────────────────────────────────
# ranger's OOB error is equivalent to leave-one-out CV for RF
oob_acc <- 1 - rf_model$prediction.error
cat(sprintf("RF OOB accuracy: %.3f\n", oob_acc))

# ── Holdout prediction ────────────────────────────────────────
rf_hold_pred  <- predict(rf_model, data = X_hold_scaled)
rf_hold_class <- colnames(rf_hold_pred$predictions)[
  apply(rf_hold_pred$predictions, 1, which.max)
]
rf_hold_acc   <- mean(rf_hold_class == as.character(y_hold))
cat(sprintf("RF holdout accuracy: %.3f\n", rf_hold_acc))

print(caret::confusionMatrix(as.factor(rf_hold_class), y_hold))

# ── Variable importance ───────────────────────────────────────
rf_imp <- sort(rf_model$variable.importance, decreasing = TRUE)

# Top 50 genes by permutation importance
rf_top50 <- data.frame(
  gene       = names(rf_imp)[1:min(50, length(rf_imp))],
  importance = rf_imp[1:min(50, length(rf_imp))]
)

plt_rf_imp <- ggplot(rf_top50,
                      aes(x = importance,
                          y = reorder(gene, importance))) +
  geom_col(fill = "steelblue") +
  labs(title = "Random Forest — Top 50 genes by permutation importance",
       x = "Permutation importance", y = NULL) +
  theme_bw()

ggsave(file.path(out_dir, "rf_variable_importance.pdf"),
       plt_rf_imp, width = 6, height = 8)

# ============================================================
# PART 3: OVERLAP BETWEEN sPLS-DA AND RF SIGNATURES
# ============================================================
cat("\n── Comparing signatures ────────────────────────────────\n")

# Top RF genes (importance > 0, ranked)
rf_genes_ranked <- names(rf_imp[rf_imp > 0])

# Overlap with sPLS-DA genes
overlap_genes <- intersect(all_splsda_genes, rf_genes_ranked)
cat(sprintf("sPLS-DA genes: %d\n", length(all_splsda_genes)))
cat(sprintf("RF genes (importance > 0): %d\n", length(rf_genes_ranked)))
cat(sprintf("Overlap: %d genes\n", length(overlap_genes)))
cat("Overlapping genes:\n")
print(overlap_genes)

# Venn-style overlap summary between the top 100 RF genes and all sPLS-DA genes
overlap_df <- data.frame(
  gene          = union(all_splsda_genes, rf_genes_ranked[1:100]),
  in_splsda     = union(all_splsda_genes, rf_genes_ranked[1:100]) %in% all_splsda_genes,
  in_rf_top100  = union(all_splsda_genes, rf_genes_ranked[1:100]) %in% rf_genes_ranked[1:100]
) %>%
  mutate(in_both = in_splsda & in_rf_top100)

write.csv(overlap_df,
          file.path(out_dir, "splsda_rf_gene_overlap.csv"),
          row.names = FALSE)

overlap_genes <- overlap_df %>%
  filter(in_both) %>%
  pull(gene)

# Heatmap of consensus signature (genes in both)
if (length(overlap_genes) >= 2) {
  X_consensus <- X_scaled[, overlap_genes, drop = FALSE]

  anno_row <- data.frame(NMF_cluster = y_train,
                          row.names   = rownames(X_consensus))

  pheatmap(
    t(X_consensus),
    annotation_col = anno_row,
    cluster_cols   = FALSE,
    cluster_rows   = TRUE,
    show_colnames  = FALSE,
    fontsize       = 12,
    fontsize_row   = 5,
    fontsize_col   = 8,
    scale          = "column",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    main           = sprintf("Consensus signature (%d genes)", length(overlap_genes)),
    filename       = file.path(out_dir, "consensus_signature_heatmap.pdf"),
    width          = 12,
    height         = 6
  )
}

# ============================================================
# PART 4: PERMUTATION TEST (shared, for both models)
# ============================================================
cat("\n── Running permutation test ────────────────────────────\n")

n_permutations <- 100
perm_acc_splsda <- numeric(n_permutations)
perm_acc_rf     <- numeric(n_permutations)

for (i in seq_len(n_permutations)) {
  set.seed(seed + i)
  y_perm <- sample(y_train)

  # RF permutation (fast — use OOB)
  rf_perm <- ranger::ranger(
    x          = X_scaled,
    y          = y_perm,
    num.trees  = 500,
    probability = TRUE,
    num.threads = parallel::detectCores() - 1,
    seed        = seed + i
  )
  perm_acc_rf[i] <- 1 - rf_perm$prediction.error

  # sPLS-DA permutation (single fit, no tuning for speed)
  splsda_perm <- try(mixOmics::splsda(
    X     = X_scaled,
    Y     = y_perm,
    ncomp = optimal_ncomp,
    keepX = optimal_keepX
  ), silent = TRUE)

  if (!inherits(splsda_perm, "try-error")) {
    perm_pred <- predict(splsda_perm,
                          newdata = X_scaled,
                          dist    = "max.dist")
    perm_class <- perm_pred$class$max.dist[, optimal_ncomp]
    perm_acc_splsda[i] <- mean(perm_class == y_perm)
  } else {
    perm_acc_splsda[i] <- NA
  }

  if (i %% 10 == 0) cat(sprintf("  Permutation %d / %d\n", i, n_permutations))
}

# ── Permutation plot ──────────────────────────────────────────
n_classes <- length(levels(y_train))
chance    <- 1 / n_classes

perm_plot_df <- bind_rows(
  data.frame(acc = perm_acc_rf,               method = "Random Forest"),
  data.frame(acc = perm_acc_splsda,           method = "sPLS-DA")
)

real_acc_df <- data.frame(
  method = c("Random Forest", "sPLS-DA"),
  acc    = c(oob_acc,          1 - min(splsda_perf$error.rate$overall))
)

plt_perm <- ggplot(perm_plot_df, aes(x = acc, fill = method)) +
  geom_histogram(bins = 30, alpha = 0.2, 
        position = "identity", colour = "white") +
  geom_vline(data = real_acc_df,
             aes(xintercept = acc, colour = method),
             linewidth = 1.2, linetype = "dashed") +
  geom_vline(xintercept = chance, linetype = "dashed",
             colour = "black", linewidth = 0.8) +
  facet_wrap(~method, ncol = 1) +
  labs(title    = "Permutation test — NMF cluster classification",
       subtitle = sprintf("Dashed black line = chance (1/%d = %.2f)\nDashed colored lines = observed accuracies", n_classes, chance),
       x = "Accuracy (permuted labels)", y = "Count") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "permutation_test.pdf"),
       plt_perm, width = 6, height = 6)

# Empirical p-values
p_rf     <- mean(perm_acc_rf     >= oob_acc,    na.rm = TRUE)
p_splsda <- mean(perm_acc_splsda >= (1 - min(splsda_perf$error.rate$overall)),
                  na.rm = TRUE)
cat(sprintf("\nRF permutation p-value:     %.4f\n", p_rf))
cat(sprintf("sPLS-DA permutation p-value: %.4f\n", p_splsda))

#create a confusion matrix for the holdout set predictions of both models
cat("\n── Holdout set confusion matrices ───────────────────────\n")
cat("\nRandom Forest:\n")
print(caret::confusionMatrix(as.factor(rf_hold_class), y_hold))
cat("\nsPLS-DA:\n")
print(caret::confusionMatrix(as.factor(splsda_hold_class), y_hold))

# ============================================================
# SAVE RESULTS
# ============================================================
saveRDS(list(
  splsda_model   = splsda_model,
  splsda_perf    = splsda_perf,
  splsda_genes   = splsda_genes,
  rf_model       = rf_model,
  rf_importance  = rf_imp,
  overlap_genes  = overlap_genes,
  perm_acc_rf    = perm_acc_rf,
  perm_acc_splsda = perm_acc_splsda,
  p_rf           = p_rf,
  p_splsda       = p_splsda,
  holdout_acc    = c(splsda = splsda_hold_acc, rf = rf_hold_acc)
), file.path(out_dir, "splsda_rf_results.rds"))

cat("\nDone. All outputs saved to:", out_dir, "\n")