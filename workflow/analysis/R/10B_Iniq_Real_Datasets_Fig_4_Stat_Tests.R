library(data.table)
library(tidyverse)
library(reshape2)

### NOTE THAT THIS RUN DOES NOT USE LIGER ###

# Helper functions 
`%ni%` <- Negate(`%in%`)

# Set seed for any sampling done 
set.seed(42)

# Change to results dir for lowcap modified
setwd("../../../results/lowcap_modified/")

# Load in and concatenate celltype imbalance summary files
setwd("celltype_imbalance_summaries")
cimba_files <- list.files()
pbmc_2_batch_imba_cimba_files <- grep(
  "pbmc_2_batch",
  cimba_files,
  value = TRUE
)
pbmc_2_batch_imba_cimba_loaded <- lapply(
  pbmc_2_batch_imba_cimba_files, 
  fread
)
pbmc_2_batch_imba_cimba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_imba_cimba_loaded
)
gc()

# Load in and concatenate full imbalance summary files 
setwd("../imbalance_summaries/")
imba_files <- list.files()
pbmc_2_batch_imba_imba_files <- grep(
  "pbmc_2_batch",
  imba_files,
  value = TRUE
)
pbmc_2_batch_imba_imba_loaded <- lapply(
  pbmc_2_batch_imba_imba_files, 
  fread
)
pbmc_2_batch_imba_imba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_imba_imba_loaded
)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
pbmc_2_batch_imba_clus_files <- grep(
  "pbmc_2_batch",
  clus_files,
  value = TRUE
)
pbmc_2_batch_imba_clus_loaded <- lapply(
  pbmc_2_batch_imba_clus_files, 
  fread
)
pbmc_2_batch_imba_clus_concat <- Reduce(
  rbind, 
  pbmc_2_batch_imba_clus_loaded
)
pbmc_2_batch_imba_clus_concat <- pbmc_2_batch_imba_clus_concat[
  pbmc_2_batch_imba_clus_concat$Method != "liger"
] 
gc()

# Load in relatedness metric summary 
setwd("../relatedness_results/")
relate_files <- list.files()
pbmc_2_batch_imba_relate_file <- grep(
  "pbmc_2_batch",
  relate_files,
  value = TRUE
)
pbmc_2_batch_imba_relate_loaded <- fread(pbmc_2_batch_imba_relate_file)

# Change to results dir for control data 
setwd("../../control/")

# Load in and concatenate celltype imbalance summary files
setwd("celltype_imbalance_summaries")
cimba_files <- list.files()
pbmc_2_batch_bal_cimba_files <- grep(
  "pbmc_2_batch_base_balanced",
  cimba_files,
  value = TRUE
)
pbmc_2_batch_bal_cimba_loaded <- lapply(
  pbmc_2_batch_bal_cimba_files, 
  fread
)
pbmc_2_batch_bal_cimba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_bal_cimba_loaded
)
gc()

# Load in and concatenate full imbalance summary files 
setwd("../imbalance_summaries/")
imba_files <- list.files()
pbmc_2_batch_bal_imba_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_files,
  value = TRUE
)
pbmc_2_batch_bal_imba_loaded <- lapply(
  pbmc_2_batch_bal_imba_files, 
  fread
)
pbmc_2_batch_bal_imba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_bal_imba_loaded
)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
pbmc_2_batch_bal_clus_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_files,
  value = TRUE
)
pbmc_2_batch_bal_clus_loaded <- lapply(
  pbmc_2_batch_bal_clus_files, 
  fread
)
pbmc_2_batch_bal_clus_concat <- Reduce(
  rbind, 
  pbmc_2_batch_bal_clus_loaded
)
pbmc_2_batch_bal_clus_concat <- pbmc_2_batch_bal_clus_concat[
  pbmc_2_batch_bal_clus_concat$Method != "liger"
]
gc()

# Load in relatedness metric summary 
setwd("../relatedness_results/")
relate_files <- list.files()
pbmc_2_batch_bal_relate_file <- grep(
  "pbmc_2_batch_base_balanced",
  relate_files,
  value = TRUE
)
pbmc_2_batch_bal_relate_loaded <- fread(pbmc_2_batch_bal_relate_file)

# Change to top level dir 
setwd("../../..")

# Create directory for lowcap modified stats test if it doesn't exist
if (!dir.exists("outs/lowcap_modified/results")) {
  dir.create("outs/lowcap_modified/results", recursive = TRUE)
}

### TEST 1 - Difference in ARI (batch and celltype) between the balanced- 
# and imbalanced PBMC 2 batch cohorts ###

# Format the relatedness metric results for both the balanced and imbalanced
# PBMC 2 batch datasets 
pbmc_2_batch_bal_relate_formatted <- pbmc_2_batch_bal_relate_loaded %>%
  group_by(`Celltype 1`, `Celltype 2`) %>%
  summarize(`Average PCA cosine dist` = mean(`PCA cosine dist`)) %>%
  as.data.frame() 

pbmc_2_batch_imba_relate_formatted <- pbmc_2_batch_imba_relate_loaded %>%
  group_by(`Celltype 1`, `Celltype 2`) %>%
  summarize(`Average PCA cosine dist` = mean(`PCA cosine dist`)) %>%
  as.data.frame() 

# Merge and subset balanced clustering data to only include control experiments 
imba_clus_merged_bal <- merge(
  pbmc_2_batch_bal_clus_concat,
  pbmc_2_batch_bal_imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_clus_merged_bal$type <- ifelse(
  imba_clus_merged_bal$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged_bal$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)
imba_clus_merged_bal_control <- imba_clus_merged_bal[
  imba_clus_merged_bal$type %in% c("Control")
]

# Merge imbalanced clustering data
imba_clus_merged_imba <- merge(
  pbmc_2_batch_imba_clus_concat,
  pbmc_2_batch_imba_imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Create dataframe of ARI celltype and batch values for both datasets 
imba_clus_imba_sub <- imba_clus_merged_imba[,
  c("Method", "Celltype ARI Imbalanced", "Batch ARI")
]
imba_clus_imba_sub$Dataset <- "PBMC full dataset imbalanced"

imba_clus_bal_sub <- imba_clus_merged_bal_control[,
  c("Method", "Celltype ARI Imbalanced", "Batch ARI")
]
imba_clus_bal_sub$Dataset <- "PBMC control dataset balanced"

imba_clus_imba_bal_sub_merged <- rbind(
  imba_clus_imba_sub, imba_clus_bal_sub
)
colnames(imba_clus_imba_bal_sub_merged) <- c(
  "method", "celltype_ari", "batch_ari", "dataset"
)

# Perform ANOVA tests for both the batch and celltype ARI 

# Celltype ARI ANOVA 
model_fit_celltype_ari_pbmc_2 <- lm(
  as.formula(
    paste0(
      "celltype_ari",
      "~",
      "method+",
      "dataset"
    )
  ),
  data = imba_clus_imba_bal_sub_merged
)
anova_result_celltype_ari_pbmc_2 <- anova(
  model_fit_celltype_ari_pbmc_2, test = "F"
)
anova_result_dt_celltype_ari_pbmc_2 <- as.data.table(
  anova_result_celltype_ari_pbmc_2,
  keep.rownames = TRUE
)
colnames(anova_result_dt_celltype_ari_pbmc_2)[1] <- "Covariate"
anova_result_dt_celltype_ari_pbmc_2$dataset_name <- "PBMC 2 batch bal/imbal"
anova_result_dt_celltype_ari_pbmc_2$metric <- "Celltype ARI"
anova_result_dt_celltype_ari_pbmc_2$last_covariate <- "Dataset"

# Batch ARI ANOVA 
model_fit_batch_ari_pbmc_2 <- lm(
  as.formula(
    paste0(
      "batch_ari",
      "~",
      "method+",
      "dataset"
    )
  ),
  data = imba_clus_imba_bal_sub_merged
)
anova_result_batch_ari_pbmc_2 <- anova(
  model_fit_batch_ari_pbmc_2, test = "F"
)
anova_result_dt_batch_ari_pbmc_2 <- as.data.table(
  anova_result_batch_ari_pbmc_2,
  keep.rownames = TRUE
)
colnames(anova_result_dt_batch_ari_pbmc_2)[1] <- "Covariate"
anova_result_dt_batch_ari_pbmc_2$dataset_name <- "PBMC 2 batch bal/imbal"
anova_result_dt_batch_ari_pbmc_2$metric <- "Batch ARI"
anova_result_dt_batch_ari_pbmc_2$last_covariate <- "Dataset"


# Combine the two and do a p-value correction 
anova_result_dt_batch_celltype_ari_pbmc_2 <- rbind(
  anova_result_dt_celltype_ari_pbmc_2,
  anova_result_dt_batch_ari_pbmc_2
)
anova_result_dt_batch_celltype_ari_pbmc_2$q_val <- p.adjust(
  anova_result_dt_batch_celltype_ari_pbmc_2$`Pr(>F)`,
  method = "BH"
)

# Save the results
fwrite(
  anova_result_dt_batch_celltype_ari_pbmc_2,
  "outs/lowcap_modified/results/10B_pbmc_2_bal_imbal_c_b_ari_anova_tests.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# 