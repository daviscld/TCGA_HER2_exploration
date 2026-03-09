This repository provides the code for a sample TCGA analysis looking at HER2+ samples, among others. Little information is provided here, as this repository is meant for only people who already have access to the analysis prompts.

Plots and tables not included in the official writeup are included here for completeness.

Run all scripts in order:

01_Load_and_QC.r

02_DEG_analysis_of_HER2clinical.r

03_HER2_clinical_multimodal_definitions.r

04_rna_vs_dna.r

05_unsupervised_clustering_UMAP_NMF.r

06_cluster_scoring_by_signatures.r

07_sPLSDA_and_rf.r

These scripts were developed and run in Visual Studio Code running R version 4.3.1 (2023-06-16). Running under: macOS Sonoma 14.2.1
