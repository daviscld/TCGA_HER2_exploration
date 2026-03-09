# TCGA HER2 Analysis

This repository contains code for an exploratory analysis of the TCGA breast cancer cohort focusing on HER2-positive tumors. The analysis integrates RNA-seq expression, copy-number variation, and clinical metadata.

The repository accompanies a written analysis submission; therefore the README is intentionally brief and assumes familiarity with the analysis prompts.

Additional plots and tables not included in the written report are included here for completeness.

---

## Repository structure

Scripts are intended to be run sequentially:

01_Load_and_QC.r

02_DEG_analysis_of_HER2clinical.r

03_HER2_clinical_multimodal_definitions.r

04_rna_vs_dna.r

05_unsupervised_clustering_UMAP_NMF.r

06_cluster_scoring_by_signatures.r

07_sPLSDA_and_rf.r

---

## Environment

Scripts were developed and executed using:

R 4.3.1 (2023-06-16)
macOS Sonoma 14.2.1

---
