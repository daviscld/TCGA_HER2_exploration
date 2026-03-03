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
#check dimensions of datasets after merging
dim(clin_data)
dim(erbb2_cn)
dim(combined_meta)
#check that patient_id.x and patient_id.y are the same in combined_meta
table(combined_meta$patient_id.x == combined_meta$patient_id.y)
#which patient_id has all the NA values?
table(is.na(combined_meta$patient_id.x))
table(is.na(combined_meta$patient_id.y))
#looks like patient_id.y has all the NA values, 
#so we'll keep patient_id.x and drop patient_id.y
combined_meta$patient_id <- combined_meta$patient_id.x
#drop patient_id.x and patient_id.y
combined_meta <- combined_meta[, !colnames(combined_meta) %in% c("patient_id.x", "patient_id.y")]

#combine datasets by patient_id, sample_id not present in rsem
#only take primary tumors for now
combined_data_primary <- merge(combined_meta[combined_meta$sample_type == "Primary",], 
                            rsem[rsem$sample_type == "Primary Tumor", ], 
                            by.x = "patient_id", 
                            by.y = "patient_id", all = TRUE)


#check replicate patient_id and sample_id columns
check_duplicate_patient_ids <- function(df, patient_id_col = "Patient.ID") {
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

check_duplicate_patient_ids(combined_data_primary)

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
        mean_age = mean(Diagnosis.Age, na.rm = TRUE),
        sd_age = sd(Diagnosis.Age, na.rm = TRUE),
        mean_survival = mean(Survival.Time, na.rm = TRUE),
        sd_survival = sd(Survival.Time, na.rm = TRUE),
        #mean_PFS = mean(days_to_patient_progression_free, na.rm = TRUE),
        #sd_PFS = sd(days_to_patient_progression_free, na.rm = TRUE),
        .groups = 'drop'
    )

#check metadata for % missingness, remove columns with >50% missing data
lapply(combined_meta, function(x) sum(is.na(x)/length(x))*100)
colnames(combined_meta)[which(sapply(combined_meta, function(x) sum(is.na(x)/length(x))*100 > 50))]


#decide what to keep out of the combined_meta dataset for downstream analyses
#add back in all HER2 variables, even if they have >50% missingness
#remove primarily sample description variables, metadata not directly
#related to question at hand, and variables with >50% missingness
final_meta <- combined_meta %>%
    select(patient_id, sample_id, erbb2_copy_number, 
    Diagnosis.Age, Sex, Cancer.Type.Detailed,
    Disease.Free..Months., Disease.Free.Status,
    ER.Status.By.IHC, PR.status.by.ihc,
    Fraction.Genome.Altered, 
    Neoplasm.Disease.Stage.American.Joint.Committee.on.Cancer.Code,
    #everything HER2
    HER2.cent17.ratio, HER2.copy.number, HER2.fish.status,
    HER2.ihc.score, Prior.Cancer.Diagnosis.Occurence,
    IHC.HER2, IHC.Score, Mutation.Count, Oncotree.Code,
    Overall.Survival..Months., Overall.Survival.Status,
    Race.Category, Sex, Person.Neoplasm.Status,
    Menopause.Status
    ) 
#note: rejected a ton of made-up hallucinated variables

#NOTE: DO NOT YET HAVE HER2 STATUS DEFINED

#iterate through final_meta and create plots to show distribution of each variable, stratified by HER2 status
#for numeric variables, create boxplots; for categorical variables, create bar plots
numeric_vars <- c("Diagnosis.Age", "Disease.Free..Months.", "Fraction.Genome.Altered", 
            "HER2.cent17.ratio", "HER2.copy.number", "IHC.Score", "Mutation.Count", 
            "Overall.Survival..Months.")
categorical_vars <- setdiff(colnames(final_meta), c("patient_id", "sample_id", 
            numeric_vars))
# Create boxplots for numeric variables
for (var in numeric_vars) { 
    p <- ggplot(final_meta, aes_string(x = "HER2_group", y = var)) +
        geom_boxplot() +
        labs(title = paste("Distribution of", var, "by HER2 status"),
             x = "HER2 Status", y = var) +
        theme_minimal()
    print(p)
}   
# Create bar plots for categorical variables
for (var in categorical_vars) {
    p <- ggplot(final_meta, aes_string(x = var, fill = "HER2_group")) +
        geom_bar(position = "dodge") +
        labs(title = paste("Distribution of", var, "by HER2 status"),
             x = var, y = "Count") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
}