library(tidyverse)

# Read all CSV files from the specified directory
data_dir <- "/Users/christinedavis/Documents/Professional/Post PhD Job Application/Tempus_Pharma_RD_Coding_Challenge_data.v2026.1"

files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

# Read each CSV file
clin_data <- read.csv(files[grepl("clin", files)])
erbb2_cn <- read.csv(files[grepl("erbb2", files)])
rsem <- read.csv(files[grepl("rsem", files)])

#standardize column names across all datasets
standardize_columns <- function(df) {
  rename_map <- c("sample_id" = "sample_id", 
                "Sample.ID" = "sample_id",
                "Sample ID" = "sample_id",
                "Patient ID" = "patient_id",
                  "patient_id" = "patient_id",
                  "Patient.ID" = "patient_id",
                  "Sample Type" = "sample_type",
                  "Sample.Type" = "sample_type",
                  "sample_type" = "sample_type")
  
  for (old_name in names(rename_map)) {
    if (old_name %in% colnames(df)) {
      colnames(df)[colnames(df) == old_name] <- rename_map[old_name]
    }
  }
  return(df)
}


clin_data <- standardize_columns(clin_data)
erbb2_cn <- standardize_columns(erbb2_cn)
rsem <- standardize_columns(rsem)

#check that at least sample_id OR patient_id column is present in all datasets
colnames(clin_data)[which(grepl("sample_id|patient_id", colnames(clin_data)))]
colnames(erbb2_cn)[which(grepl("sample_id|patient_id", colnames(erbb2_cn)))]
colnames(rsem)[which(grepl("sample_id|patient_id", colnames(rsem)))]

#how many sample_ids in clin_data are also in erbb2_cn?
table(clin_data$sample_id %in% erbb2_cn$sample_id)
table(erbb2_cn$sample_id %in% clin_data$sample_id)
#all erbb2_cn sample_ids are in clin_data, but not vice versa. 963 overlaps
#expect 145 clin_data samples with no erb22 data

#should combine clin_data and erbb2_cn by sample_id
combined_meta <- merge(clin_data, erbb2_cn, by.x = "sample_id",   
                       by.y = "sample_id", all = TRUE)

#make all . or .. in column names into _
colnames(combined_meta) <- gsub("\\.{1,2}", "_", 
                colnames(combined_meta)) %>%
  gsub("_$", "", .)

#check dimensions of datasets after merging
dim(clin_data)
dim(erbb2_cn)
dim(combined_meta)
#check that patient_id.x and patient_id.y are the same in combined_meta
table(combined_meta$patient_id_x == combined_meta$patient_id_y)
#which patient_id has all the NA values?
table(is.na(combined_meta$patient_id_x))
table(is.na(combined_meta$patient_id_y))
#looks like patient_id.y has all the NA values, 
#so we'll keep patient_id_x and drop patient_id_y
combined_meta$patient_id <- combined_meta$patient_id_x
#drop patient_id_x and patient_id_y
combined_meta <- combined_meta[, !colnames(combined_meta) %in% c("patient_id_x", "patient_id_y")]

#check replicate patient_id and sample_id columns
check_duplicate_patient_ids <- function(df, patient_id_col = "patient_id") {
    if (any(duplicated(df[[patient_id_col]]))) {
        warning("Duplicate patient_id found in dataset.")
        dup_ids <- df[[patient_id_col]][duplicated(df[[patient_id_col]])]
        duplicated_data <- df[df[[patient_id_col]] %in% dup_ids, ]
        return(duplicated_data)
    } else {
        message("No duplicate patient_ids found.")
        return(NULL)
    }
}

#write a table of the mean and stdev for demographic data
#RETURN TO THIS once I've defined her2_status
# Create HER2 grouping variable
combined_meta <- combined_meta %>%
    mutate(HER2_group = case_when(
        grepl("Positive|\\+", her2_status, ignore.case = TRUE) ~ "HER2+",
        grepl("Negative|\\-", her2_status, ignore.case = TRUE) ~ "HER2-",
        TRUE ~ NA_character_
    ))

# Summary statistics stratified by HER2 status
#include ethnicity category, menopause status, metastatic status, mutation count, 
#overall survival months, sex, erbb2 copy number
#Stage III/IV (%)			
#ER+ (%)			
#PR+ (%)
demographics_summary <- combined_meta %>%
    group_by(HER2_group) %>%
    summarise(
        n = n(),
        mean_age = mean(Diagnosis_Age, na.rm = TRUE),
        sd_age = sd(Diagnosis_Age, na.rm = TRUE),
        mean_survival = mean(Survival_Time, na.rm = TRUE),
        sd_survival = sd(Survival_Time, na.rm = TRUE),
        #mean_PFS = mean(days_to_patient_progression_free, na.rm = TRUE),
        #sd_PFS = sd(days_to_patient_progression_free, na.rm = TRUE),
        .groups = 'drop'
    )

#check metadata for % missingness, remove columns with >50% missing data
lapply(combined_meta, function(x) sum(is.na(x)/length(x))*100)

#decide what to keep out of the combined_meta dataset for downstream analyses
#add back in all HER2 variables, even if they have >50% missingness
#remove primarily sample description variables, metadata not directly
#related to question at hand, and variables with >50% missingness
final_meta <- combined_meta %>%
    select(patient_id, sample_id, erbb2_copy_number, 
    Diagnosis_Age, Sex, Cancer_Type_Detailed,
    Disease_Free_Months, Disease_Free_Status,
    ER_Status_By_IHC, PR_status_by_ihc,
    Fraction_Genome_Altered, sample_type,
    Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code,
    #everything HER2
    HER2_cent17_ratio, HER2_copy_number, HER2_fish_status,
    HER2_ihc_score, Prior_Cancer_Diagnosis_Occurence,
    IHC_HER2, IHC_Score, Mutation_Count, Oncotree_Code,
    Overall_Survival_Months, Overall_Survival_Status,
    Race_Category, Sex, Person_Neoplasm_Status,
    Menopause_Status
    ) 
#note: rejected a ton of made-up hallucinated variables

#rename Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code
final_meta <- final_meta %>%
    rename(Cancer_Stage = 
                 Neoplasm_Disease_Stage_American_Joint_Committee_on_Cancer_Code)
#make menopause status into pre, post, or unknown
final_meta <- final_meta %>%
    mutate(Menopause_Status = case_when(
        grepl("Indeterminate", Menopause_Status) ~ "Indeterminate",
        grepl("^Pre", Menopause_Status, ignore.case = TRUE) ~ "Pre",
        grepl("^Post", Menopause_Status, ignore.case = TRUE) ~ "Post",
        grepl("^Peri", Menopause_Status, ignore.case = TRUE) ~ "Peri",
        TRUE ~ "Unknown"
    ))

#iterate through final_meta and create plots to show distribution of each variable
#for numeric variables, create histograms; for categorical variables, create pie charts
numeric_vars <- c("Diagnosis_Age", "Disease_Free_Months", 
          "Fraction_Genome_Altered", 
            "HER2_cent17_ratio", "HER2_copy_number", "IHC_Score", 
            "Mutation_Count", 
            "erbb2_copy_number",
            "Overall_Survival_Months")
categorical_vars <- setdiff(colnames(final_meta), c("patient_id", "sample_id", 
            numeric_vars))
# Create histograms for numeric variables
for (var in numeric_vars) {
    pdf(paste0(var,"_histogram.pdf"))
    print(hist(as.numeric(final_meta[,var]),
    main = var,
    xlab = "Value"))
    dev.off()
}
# Create pie charts for categorical variables
for (var in categorical_vars) {
    p <- final_meta %>%
      count(.data[[var]]) %>%
      mutate(pct = n / sum(n) * 100) %>%
      ggplot(aes(x = "", y = n, fill = .data[[var]])) +
      geom_bar(stat = "identity", width = 1) +
      geom_text(aes(label = paste0(round(pct, 1), "%")), 
                position = position_stack(vjust = 0.5)) +
      coord_polar("y", start = 0) +
      labs(title = paste("Distribution of ", var)) +
      theme_minimal() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank())
    pdf(paste0(var,"_piechart.pdf"))
    print(p)
    dev.off()
}

#combine datasets by patient_id, sample_id not present in rsem
#only take primary tumors for now
combined_data_primary <- merge(final_meta[final_meta$sample_type == "Primary",], 
                            rsem[rsem$sample_type == "Primary Tumor", ], 
                            by.x = "patient_id", 
                            by.y = "patient_id", all = TRUE)

check_duplicate_patient_ids(combined_data_primary)

#keep sample_type.x, more complete
combined_data_primary$sample_type <- combined_data_primary$sample_type.x
#drop sample_type.x and sample_type.y
combined_data_primary <- combined_data_primary[, !colnames(combined_data_primary) %in% 
            c("sample_type.x", "sample_type.y")]

#look at distribution of erbb2 gene expression as well
pdf("erbb2_expression_histogram.pdf")
print(hist(as.numeric(combined_data_primary$ERBB2),
  main = "ERBB2 Expression",
  xlab = "Value"))
dev.off()

#add histogram zoomed in on 0-100000 range of erbb2 expression
pdf("erbb2_expression_histogram_zoomed.pdf")
print(hist(as.numeric(
  combined_data_primary$ERBB2[which(as.numeric(
    combined_data_primary$ERBB2) < 100000)]),
       breaks = 10,
       main = "ERBB2 Expression (0-100000 range)",
       xlab = "Value"))
dev.off()

#below generated by Claude

# ============================================================
# TCGA RSEM Counts: QC → DESeq2 Object
# ============================================================
# Expects: TCGA RSEM "expected_count" matrix (genes x samples)
# from TCGAbiolinks or GDC download. DESeq2 requires integer
# counts — RSEM expected_counts are floats and must be rounded.
# ============================================================

library(TCGAbiolinks)
library(DESeq2)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(pheatmap)
library(RColorBrewer)

# ── 0. CONFIG ──────────────────────────────────────────────
project      <- "TCGA-BRCA-HER2"
#DATA_DIR     <- "data/tcga_raw"
qc_dir       <- "qc_plots"
min_count    <- 10        # min total count across samples to keep gene
min_samples  <- 0.1      # gene must be expressed in at least X% of samples
min_genes    <- 5000     # flag samples with fewer detected genes than this
max_mito     <- 0.3      # max fraction mitochondrial reads (if annotated)
outlier_sd   <- 3        # SD cutoff for library-size outlier flagging

#dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir,   showWarnings = FALSE, recursive = TRUE)

#split combined_data_primary into metadata and counts matrix
#metadata should be all columns from final_meta, plus bcr_patient_barcode from rsem
#counts matrix should be all numeric columns from rsem
meta_cols <- c(colnames(final_meta), "bcr_patient_barcode")
meta <- combined_data_primary %>%
  dplyr::select(any_of(meta_cols))

counts_raw <- combined_data_primary %>%
  dplyr::select(-any_of(meta_cols)) %>%
    as.matrix() %>%
    round()  # RSEM expected_counts are floats; DESeq2 needs integers
    #double-checked validity of this on BioStars

# Transpose so genes are rows and samples are columns 
counts_raw <- t(counts_raw)
colnames(counts_raw) <- meta$patient_id  # Set sample IDs as column names

# Ensure metadata is ordered to match the sample order in counts_raw
# Just set this, so somewhat redundant 
meta <- meta[match(colnames(counts_raw), meta$patient_id), ]
rownames(meta) <- meta$patient_id  # Set sample IDs as rownames for easy indexing

#double-check counts_raw is numeric
if (!is.numeric(counts_raw)) {
  stop("Error: counts matrix must be numeric.")
}

#Check to make sure all samples have RNAseq data
if (any(is.na(counts_raw))) {
  warning("Warning: NA values found in counts matrix. Check for missing data.")
}
#remove any samples with all NA values in counts matrix
na_samples <- colnames(counts_raw)[colSums(is.na(counts_raw)) == nrow(counts_raw)]
if (length(na_samples) > 0) {
  warning(paste("Warning: Removing samples with all NA values in counts matrix:", paste(na_samples, collapse = ", ")))
  counts_raw <- counts_raw[, !colnames(counts_raw) %in% na_samples]
  meta <- meta[!rownames(meta) %in% na_samples, ]
}

# ── QC — STEP A: Gene-level filtering ───────────────────
cat("\n── Gene-level QC ──\n")
n_genes_raw <- nrow(counts_raw)

# b. Flag mitochondrial genes (optional — usually low in TCGA bulk)
mito_genes <- str_detect(rownames(counts_raw), "^MT-") #none detected
mito_frac  <- colSums(counts_raw[mito_genes, ]) / colSums(counts_raw)
meta$mito_frac <- mito_frac
cat(" Max mitochondrial fraction:", round(max(mito_frac), 4), "\n")

# c. Minimum count filter:
#     Keep genes with >= min_count total counts across >= min_samples% of samples
min_samp <- ceiling(min_samples * ncol(counts_raw))
expressed <- rowSums(counts_raw >= min_count) >= min_samp
counts_filtered <- counts_raw[expressed, ]

cat(sprintf("  Genes before filter : %d\n", n_genes_raw))
cat(sprintf("  Genes after filter  : %d (removed %d)\n",
            nrow(counts_filtered), n_genes_raw - nrow(counts_filtered)))

# ── QC — STEP B: Sample-level QC metrics ────────────────
cat("\n── Sample-level QC ──\n")

meta$lib_size       <- colSums(counts_raw)          # total mapped reads
meta$n_genes_expr   <- colSums(counts_raw > 0)      # genes with any count
meta$log10_lib_size <- log10(meta$lib_size)

# Flag low-library-size outliers (> outlier_sd SDs below mean on log scale)
lib_mean <- mean(meta$log10_lib_size)
lib_sd   <- sd(meta$log10_lib_size)
meta$low_lib_flag <- meta$log10_lib_size < (lib_mean - outlier_sd * lib_sd)

# Flag low gene-count samples
meta$low_gene_flag <- meta$n_genes_expr < min_genes

# Flag high mitochondrial fraction
meta$high_mito_flag <- meta$mito_frac > max_mito

# Combined QC fail flag
meta$qc_fail <- meta$low_lib_flag | meta$low_gene_flag | meta$high_mito_flag

n_fail <- sum(meta$qc_fail, na.rm = TRUE)
cat(sprintf("  Samples total       : %d\n", ncol(counts_raw)))
cat(sprintf("  Low library size    : %d\n", sum(meta$low_lib_flag, na.rm=TRUE)))
cat(sprintf("  Low gene count      : %d\n", sum(meta$low_gene_flag, na.rm=TRUE)))
cat(sprintf("  High mito fraction  : %d\n", sum(meta$high_mito_flag, na.rm=TRUE)))
cat(sprintf("  QC FAIL (total)     : %d\n", n_fail))

# ── QC PLOTS ────────────────────────────────────────────

# a. Library size distribution
p_lib <- ggplot(meta, aes(x = log10_lib_size)) +
  geom_histogram(bins = 60, alpha = 0.7, color = "white") +
  geom_vline(xintercept = lib_mean - outlier_sd * lib_sd,
             linetype = "dashed", color = "firebrick") +
  labs(title = "Library Size Distribution (log10 counts)",
       x = "log10(total counts)", y = "# Samples",
       caption = "Dashed line: 3 SD below mean (outlier cutoff)") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(qc_dir, "01_library_size.png"), p_lib, 
          width=8, height=4, dpi=300)

# b. Genes detected vs library size (scatter)
p_genes <- ggplot(meta, aes(x = log10_lib_size, y = n_genes_expr,
                             color = qc_fail)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                     labels = c("Pass", "Fail")) +
  geom_hline(yintercept = min_genes, linetype = "dashed", color = "firebrick") +
  labs(title = "Genes Detected vs Library Size",
       x = "log10(total counts)", y = "# Genes Detected",
       color = "QC status", shape = "Sample type") +
  theme_bw(base_size = 12)
ggsave(file.path(qc_dir, "02_genes_vs_libsize.png"), p_genes, 
          width=7, height=5, dpi=300  )


# c. Sample correlation heatmap on VST-normalised subset (top 500 variable genes)
#     Run on QC-passing tumor samples only (sample up to 200 for speed)
pass_tumors <- meta %>%
  filter(!qc_fail) %>%
  slice_sample(n = min(200, sum(!meta$qc_fail)))

counts_vst_input <- counts_filtered[, rownames(pass_tumors)]

# Quick VST without a full DESeq model (blind=TRUE for QC purposes)
#blind = TRUE ignores the design when estimating dispersion for transformation
vst_qc <- vst(counts_vst_input, blind = TRUE)
top500  <- head(order(rowVars(vst_qc), decreasing = TRUE), 500)
cor_mat <- cor(vst_qc[top500, ], method = "pearson")
plt <- pheatmap(cor_mat,
        annotation_col = pass_tumors %>%
          select(PR_status_by_ihc, ER_Status_By_IHC,
             Cancer_Stage, HER2_ihc_score, IHC_HER2,
             HER2_fish_status) %>%
          as.data.frame(),
        show_rownames = FALSE,
        show_colnames = FALSE,
        color = colorRampPalette(rev(brewer.pal(9,
                             "RdYlBu")))(100),
        main = "Sample-Sample Pearson Correlation (top 500 variable genes)",
        fontsize = 9)
ggsave(file.path(qc_dir, "03_sample_correlation_heatmap.png"), plt, 
      width=9, height=9)

# ── REMOVE QC-FAILING SAMPLES ───────────────────────────
pass_samples   <- rownames(meta)[!meta$qc_fail]
counts_qc      <- counts_filtered[, pass_samples]
meta_qc        <- meta[pass_samples, ]

cat(sprintf("\n── After QC: %d samples × %d genes retained\n",
            ncol(counts_qc), nrow(counts_qc)))

  #Use FISH when we have it, IHC otherwise, to create a HER2 clinical variable
  meta_qc <- meta_qc %>%
  mutate(
    her2_source = if_else(!is.na(HER2_fish_status), HER2_fish_status, IHC_HER2),
    her2_clinical = case_when(
      her2_source == "Positive" ~ "HER2pos",
      her2_source == "Negative" ~ "HER2neg",
      her2_source == "Indeterminate" ~ "HER2indeterminate",
      her2_source == "Equivocal" ~ "HER2equivocal",
      TRUE ~ "Unknown"
    )
  )

#Retain only HER2pos, HER2neg for now- unassessed/indeterminate/equivocal 
#samples will be excluded from DE testing but can be included later

# ── BUILD DESeq2 OBJECT ─────────────────────────────────
# Adjust the design formula to look at her2_clinical status

# Ensure design variable is a factor with a meaningful reference level
meta_qc$her2_clinical <- factor(meta_qc$her2_clinical,
                                 levels = c("HER2neg", "HER2pos"))

# Remove samples where design variable is NA
keep <- !is.na(meta_qc$her2_clinical) & meta_qc$her2_clinical != "Unknown"
counts_deseq <- counts_qc[, meta_qc$patient_id[keep]]
meta_deseq   <- meta_qc[keep, ] %>% as.data.frame()
rownames(meta_deseq) <- meta_deseq$patient_id

cat(sprintf("\n── After Retaining Only HER2-assigned: %d samples × %d genes retained\n",
            ncol(counts_deseq), nrow(counts_deseq)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_deseq,
  colData   = meta_deseq,
  design    = ~ her2_clinical       # contrast of interest
)

# ── PRE-FILTER DESeq OBJECT (optional redundant safety filter) ──
# DESeq2 docs recommend removing genes with very low counts before
# estimating dispersions — already done above but enforce here too.

#this removes a few more genes, now that we have subset down samples
#to those with her2_clinical status, but doesn't change sample count
keep_genes <- rowSums(counts(dds) >= 10) >= ceiling(0.1 * ncol(dds))
dds <- dds[keep_genes, ]

cat(sprintf("\n── DESeq2 dataset built: %d samples × %d genes\n",
            ncol(dds), nrow(dds)))
print(dds)

# ── SIZE FACTOR ESTIMATION & VST FOR DOWNSTREAM USE ────
dds <- estimateSizeFactors(dds) #add control for reading depth/library size 

# VST for PCA / clustering / heatmaps (NOT for DE testing — use dds for that)
vst_final <- vst(dds, blind = FALSE)   # blind=FALSE uses the design

# Quick PCA to confirm QC looks clean

pca_data <- prcomp(t(assay(vst_final)), scale. = FALSE)

pca_df <- data.frame(cbind(as.data.frame(pca_data[["x"]]), 
              patient_id = rownames(pca_data[["x"]])))

pca_df <- as.data.frame(pca_df) %>%
  left_join(meta_deseq, by = "patient_id")

pct_var <- round(100 * (pca_data$sdev^2 / sum(pca_data$sdev^2))[1:2])

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(shape = her2_clinical, 
                color = her2_clinical), 
                size = 2.5, alpha = 0.8, stroke = 1) +
xlab(paste0("PC1: ", pct_var[1], "% variance")) +
  ylab(paste0("PC2: ", pct_var[2], "% variance")) +
  labs(title = "PCA — VST normalised counts (post-QC)") +
  theme_bw(base_size = 12)
ggsave(file.path(qc_dir, "05_pca_post_qc.png"), p_pca, 
      width=7, height=5, dpi=300)

# ── SAVE OUTPUTS ───────────────────────────────────────
saveRDS(dds,       "dds_her2_brca.rds")
saveRDS(vst_final, "vst_her2_brca.rds")
write.csv(meta_deseq, "metadata_qc_passed.csv", row.names = TRUE)

cat("\n✓ Pipeline complete. Outputs:\n")
cat("  dds_her2_brca.rds        — DESeq2 dataset (use for DE testing)\n")
cat("  vst_her2_brca.rds        — VST-normalised SE (use for PCA/clustering)\n")
cat("  metadata_qc_passed.csv   — cleaned sample metadata\n")
cat("  qc_plots/                — all QC figures\n")