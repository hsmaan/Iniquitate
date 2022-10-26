library(data.table)
library(tidyverse)
library(reshape2)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Change to results dir for control data 
setwd("../../../results/control/")

# Set seed for any sampling done 
set.seed(42)

##### Statistical analysis of PBMC 2 batch balanced data - baseline Fig 2 #####

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_files,
  value = TRUE
)
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
clus_concord_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_concord_files,
  value = TRUE
)
clus_concord_loaded <- lapply(clus_concord_files, fread)
clus_concord_concat <- Reduce(rbind, clus_concord_loaded)

# Load in and concatenate dge concordance summaries
setwd("../dge_concord_stats/")
dge_files <- list.files()
dge_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_files,
  value = TRUE
)
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_files <- grep(
  "pbmc_2_batch_base_balanced",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)

# Change to top level dir 
setwd("../../..")

### Statistical tests for ARI/AMI/etc. ### 

# Merge imbalance and clustering summary results
imba_clus_merged <- merge(
  clus_concat,
  imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_clus_unique_celltypes <- unique(imba_clus_merged$`Downsampled celltypes`)
imba_clus_unique_celltypes <- imba_clus_unique_celltypes[
  imba_clus_unique_celltypes %ni% "None"
]
imba_clus_none_celltype_len <- length(
  which(
    imba_clus_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_clus_none_celltype_draw <- sample(
  imba_clus_unique_celltypes,
  imba_clus_none_celltype_len,
  replace = TRUE
)
imba_clus_merged$`Downsampled celltypes`[
  which(
    imba_clus_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_clus_none_celltype_draw

# Format celltype names 
imba_clus_merged$`Downsampled celltypes` <- plyr::mapvalues(
  imba_clus_merged$`Downsampled celltypes`,
  from = c(
    "Monocyte_CD14",
    "Monocyte_FCGR3A",
    "CD4 T cell",
    "CD8 T cell"
  ),
  to = c(
    "CD14+ Monocyte",
    "FCGR3A+ Monocyte",
    "CD4+ T cell",
    "CD8+ T cell"
  )
)

# Format column names
colnames(imba_clus_merged) <- plyr::mapvalues(
  colnames(imba_clus_merged),
  from = c(
    "Celltype ARI Imbalanced",
    "Celltype AMI Imbalanced",
    "Celltype Homogeneity Imbalanced",
    "Celltype Completeness Imbalanced"
  ),
  to = c(
    "Celltype ARI",
    "Celltype AMI",
    "Celltype Homogeneity",
    "Celltype Completeness"
  )
)

# Indicate which samples are controls and which are real runs
imba_clus_merged$type <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
) 

# Create a function to do a wilcoxon rank sum test across the given
# subsets for the indicated metric 
wilcox_metric_test <- function(
    dataset, 
    dataset_name,
    subset_col,
    subset_1, 
    subset_2, 
    metric
  ) {
  subset_1_sub <- dataset[dataset[[subset_col]] %in% subset_1][[metric]]
  subset_2_sub <- dataset[dataset[[subset_col]] %in% subset_2][[metric]]
  wrs_p_value <- wilcox.test(subset_1_sub, subset_2_sub)$p.value
  wrs_summary_df <- data.frame(
    "Dataset" = dataset_name,
    "Subset_name" = subset_col,
    "Subset_1" = subset_1,
    "Subset_2" = subset_2,
    "Metric" = metric,
    "Wilcoxon_Rank_Sum_p" = wrs_p_value,
    "Samples_subset_1" = length(subset_1_sub),
    "Samples_subset_2" = length(subset_2_sub)
  )
  return(wrs_summary_df)
}

# Do a wilcoxin rank sum test between all control, ablated, and downsampled
# experiment results for all metrics 
metrics <- list(
  "Celltype ARI",
  "Celltype AMI",
  "Celltype Homogeneity",
  "Celltype Completeness",
  "Batch ARI",
  "Batch AMI",
  "Batch Homogeneity",
  "Batch Completeness"
)
comparison_subsets <- list(
  c("Control", "Downsampled"),
  c("Control", "Ablated")
)

# Get all results for batch and celltype metrics (globally) for the two given
# comparisons 
global_metric_results <- lapply(metrics, function(x) {
  lapply(comparison_subsets, function(y) {
    wrs_summary_df <- wilcox_metric_test(
      dataset = imba_clus_merged,
      dataset_name = "PBMC 2 batch base balanced",
      subset_col = "type",
      subset_1 = y[[1]],
      subset_2 = y[[2]],
      metric = x
    )
  })
})

# Concatenate all results for type subset
global_metric_results_concat_1 <- lapply(
  global_metric_results,
  function(x) {
    x_red <- Reduce(rbind, x)
    return(x_red)
  }
)
global_metric_results_concat_2 <- Reduce(rbind, global_metric_results_concat_1)

# Reperform for pooled (downsampled, ablated) comparison

# Create new pooled type column
imba_clus_merged$type_pooled <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  "Pooled downsampled/ablated"
) 

# Get all results for batch and celltype metrics (globally) for the given
# comparison
global_metric_results_pooled <- lapply(metrics, function(x) {
  wrs_summary_df <- wilcox_metric_test(
    dataset = imba_clus_merged,
    dataset_name = "PBMC 2 batch base balanced",
    subset_col = "type_pooled",
    subset_1 = "Control",
    subset_2 = "Pooled downsampled/ablated",
    metric = x
  )
})

# Concatenate all results for type_pooled subset
global_metric_results_pooled_concat <- Reduce(
  rbind,
  global_metric_results_pooled
)

# Concatenate individual and pooled subset results and save
global_metric_results_indi_pooled <- rbind(
  global_metric_results_concat_2,
  global_metric_results_pooled_concat
)
fwrite(
  global_metric_results_indi_pooled,
  "outs/control/results/06_pbmc_base_global_metric_wilcoxon_tests.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Create function to perform ANOVA test, given either type or type pooled,
# a metric, and taking into account method and specific downsampled 
# celltype 
anova_metric_test <- function(
  dataset, 
  dataset_name,
  metric,
  last_covariate
){
  model_fit <- lm(
    as.formula(
      paste0(
        metric, 
        "~",
        "method+",
        "downsampled_celltypes+",
        last_covariate
      )
    ),
    data = dataset
  )
  anova_result <- anova(model_fit, test = "F")
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$metric <- metric
  anova_result_dt$last_covariate <- last_covariate
  return(anova_result_dt)
}

# Rename columns for use in formulas (no spaces)
colnames(imba_clus_merged) <- plyr::mapvalues(
  colnames(imba_clus_merged),
  from = c(
    "Celltype ARI",
    "Celltype AMI",
    "Celltype Homogeneity",
    "Celltype Completeness",
    "Batch ARI",
    "Batch AMI",
    "Batch Homogeneity",
    "Batch Completeness",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "celltype_ari",
    "celltype_ami",
    "celltype_homogeneity",
    "celltype_completeness",
    "batch_ari",
    "batch_ami",
    "batch_homogeneity",
    "batch_completeness",
    "method",
    "downsampled_celltypes"
  )
)

# Iterate over metrics and the important last covariate (type)
metrics_aov <- list(
  "celltype_ari",
  "celltype_ami",
  "celltype_homogeneity",
  "celltype_completeness",
  "batch_ari",
  "batch_ami",
  "batch_homogeneity",
  "batch_completeness"
)
metric_aov_results <- lapply(metrics_aov, function(x) {
  aov_summary_df <- anova_metric_test(
    dataset = imba_clus_merged,
    dataset_name = "PBMC 2 batch base balanced",
    metric = x,
    last_covariate = "type"
  )
})

# Concatenate results, perform FDR correction, and save 
metric_aov_results_concat <- Reduce(rbind, metric_aov_results)
metric_aov_results_concat$`FDR_q` <- p.adjust(
  metric_aov_results_concat$`Pr(>F)`
)
fwrite(
  metric_aov_results_concat,
  "outs/control/results/06_pbmc_base_metric_aov_results_types_ct_method_ctrl.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

### Statistical tests for KNN classification metrics for baseline data ###

# Merge imbalance and knn classification results together
imba_knn_merged <- merge(
  imba_concat,
  knn_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_knn_merged <- distinct(imba_knn_merged)

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_knn_unique_celltypes <- unique(imba_knn_merged$`Downsampled celltypes`)
imba_knn_unique_celltypes <- imba_knn_unique_celltypes[
  imba_knn_unique_celltypes %ni% "None"
]
imba_knn_none_celltype_len <- length(
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_knn_none_celltype_draw <- sample(
  imba_knn_unique_celltypes,
  imba_knn_none_celltype_len,
  replace = TRUE
)
imba_knn_merged$`Downsampled celltypes`[
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_knn_none_celltype_draw

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype$type <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Create another column for pooling the ablated and downsampled data together
imba_knn_merged_celltype$type_pooled <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  "Pooled downsampled/ablated"
)

# Create a function to do a wilcoxon rank sum test across the given
# subsets for the KNN classification F1 scores where the F1 score 
# is for the celltype that was downsampled
wilcox_knn_test <- function(
  dataset, 
  dataset_name,
  subset_col,
  subset_1, 
  subset_2, 
  metric = "F1-score"
) {
  subset_1_sub <- dataset[dataset[[subset_col]] %in% subset_1][[metric]]
  subset_2_sub <- dataset[dataset[[subset_col]] %in% subset_2][[metric]]
  wrs_p_value <- wilcox.test(subset_1_sub, subset_2_sub)$p.value
  wrs_summary_df <- data.frame(
    "Dataset" = dataset_name,
    "Subset_name" = subset_col,
    "Subset_1" = subset_1,
    "Subset_2" = subset_2,
    "Metric" = metric,
    "Wilcoxon_Rank_Sum_p" = wrs_p_value,
    "Samples_subset_1" = length(subset_1_sub),
    "Samples_subset_2" = length(subset_2_sub)
  )
  return(wrs_summary_df)
}

# Get global KNN classification wilcoxon result across type subsets
# and concatenate
global_celltype_knn_results_type <- lapply(comparison_subsets, function(x) {
    wrs_summary_df <- wilcox_knn_test(
      dataset = imba_knn_merged_celltype,
      dataset_name = "PBMC 2 batch base balanced",
      subset_col = "type",
      subset_1 = x[[1]],
      subset_2 = x[[2]],
      metric = "F1-score"
    )
})
global_celltype_knn_results_type_concat <- Reduce(
  rbind,
  global_celltype_knn_results_type
)

# Get global KNN classification wilcoxon result across type_pooled subsets
# and concatenate 
global_celltype_knn_results_type_pooled <- wilcox_knn_test(
    dataset = imba_knn_merged_celltype,
    dataset_name = "PBMC 2 batch base balanced",
    subset_col = "type_pooled",
    subset_1 = "Control",
    subset_2 = "Pooled downsampled/ablated",
    metric = "F1-score"
)

# Concatenate type and type pooled results together and save
global_celltype_knn_results_type_and_pooled <- rbind(
  global_celltype_knn_results_type_concat,
  global_celltype_knn_results_type_pooled
)
fwrite(
  global_celltype_knn_results_type_and_pooled,
  "outs/control/results/06_pbmc_base_knn_cell_ds_cell_test_wilcoxon_tests.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Create function to perform ANOVA test, given either type or type pooled,
# and taking into account method used and downsampled celltype (same as 
# the examined celltype)
anova_knn_test <- function(
  dataset, 
  dataset_name,
  last_covariate,
  metric = "f1_score"
){
  model_fit <- lm(
    as.formula(
      paste0(
        metric, 
        "~",
        "method+",
        "downsampled_celltypes+",
        last_covariate
      )
    ),
    data = dataset
  )
  anova_result <- anova(model_fit, test = "F")
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$metric <- metric
  anova_result_dt$last_covariate <- last_covariate
  return(anova_result_dt)
}

# Rename columns for use in formulas 
colnames(imba_knn_merged_celltype) <- plyr::mapvalues(
  colnames(imba_knn_merged_celltype),
  from = c(
    "F1-score",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "f1_score",
    "method",
    "downsampled_celltypes"
  )
)

# Get the ANOVA KNN classification F1 score for type after taking into account
# method and downsampled celltype, and save
knn_aov_results <- anova_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch base balanced",
  last_covariate = "type",
  metric = "f1_score"
)
fwrite(
  knn_aov_results,
  paste0(
    "outs/control/results/",
    "06_pbmc_base_knn_cell_ds_cell_aov_results_ctrl_method_celltype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

##### Statistical analysis of PBMC 2 batch balanced - hierarchical Fig 2 #####

# Clear all previous files and reload KNN functions 
rm(list = ls())
gc()

wilcox_knn_test <- function(
  dataset, 
  dataset_name,
  subset_col,
  subset_1, 
  subset_2, 
  metric = "F1-score"
) {
  subset_1_sub <- dataset[dataset[[subset_col]] %in% subset_1][[metric]]
  subset_2_sub <- dataset[dataset[[subset_col]] %in% subset_2][[metric]]
  wrs_p_value <- wilcox.test(subset_1_sub, subset_2_sub)$p.value
  wrs_summary_df <- data.frame(
    "Dataset" = dataset_name,
    "Subset_name" = subset_col,
    "Subset_1" = subset_1,
    "Subset_2" = subset_2,
    "Metric" = metric,
    "Wilcoxon_Rank_Sum_p" = wrs_p_value,
    "Samples_subset_1" = length(subset_1_sub),
    "Samples_subset_2" = length(subset_2_sub)
  )
  return(wrs_summary_df)
}

anova_knn_test <- function(
  dataset, 
  dataset_name,
  last_covariate,
  metric = "f1_score"
){
  model_fit <- lm(
    as.formula(
      paste0(
        metric, 
        "~",
        "method+",
        "downsampled_celltypes+",
        last_covariate
      )
    ),
    data = dataset
  )
  anova_result <- anova(model_fit, test = "F")
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$metric <- metric
  anova_result_dt$last_covariate <- last_covariate
  return(anova_result_dt)
}

# Reload helper functions and seed
`%ni%` <- Negate(`%in%`)

# Set seed for any sampling done 
set.seed(42)

# Load in and concatenate summary files - just imbalance and KNN results
setwd("results/control/imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)

# Change to top level dir
setwd("../../..")

# Merge imbalance and knn classification results together
imba_knn_merged <- merge(
  imba_concat,
  knn_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_knn_merged <- distinct(imba_knn_merged)

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_knn_unique_celltypes <- unique(imba_knn_merged$`Downsampled celltypes`)
imba_knn_unique_celltypes <- imba_knn_unique_celltypes[
  imba_knn_unique_celltypes %ni% "None"
]
imba_knn_none_celltype_len <- length(
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_knn_none_celltype_draw <- sample(
  imba_knn_unique_celltypes,
  imba_knn_none_celltype_len,
  replace = TRUE
)
imba_knn_merged$`Downsampled celltypes`[
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_knn_none_celltype_draw

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype$type <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Create another column for pooling the ablated and downsampled data together
imba_knn_merged_celltype$type_pooled <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  "Pooled downsampled/ablated"
)

# Get global KNN classification wilcoxon result across type subsets
# and concatenate
comparison_subsets <- list(
  c("Control", "Downsampled"),
  c("Control", "Ablated")
)

global_celltype_knn_results_type <- lapply(comparison_subsets, function(x) {
  wrs_summary_df <- wilcox_knn_test(
    dataset = imba_knn_merged_celltype,
    dataset_name = "PBMC 2 batch hierarchical balanced",
    subset_col = "type",
    subset_1 = x[[1]],
    subset_2 = x[[2]],
    metric = "F1-score"
  )
})
global_celltype_knn_results_type_concat <- Reduce(
  rbind,
  global_celltype_knn_results_type
)

# Get global KNN classification Wilcoxin result across type_pooled subsets
# and concatenate 
global_celltype_knn_results_type_pooled <- wilcox_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch hierarchical balanced",
  subset_col = "type_pooled",
  subset_1 = "Control",
  subset_2 = "Pooled downsampled/ablated",
  metric = "F1-score"
)

# Concatenate type and type pooled results together and save
global_celltype_knn_results_type_and_pooled <- rbind(
  global_celltype_knn_results_type_concat,
  global_celltype_knn_results_type_pooled
)
fwrite(
  global_celltype_knn_results_type_and_pooled,
  paste0(
    "outs/control/results/",
    "06_pbmc_hierarchical_knn_cell_ds_cell_test_wilcoxon_tests.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Rename columns for use in formulas 
colnames(imba_knn_merged_celltype) <- plyr::mapvalues(
  colnames(imba_knn_merged_celltype),
  from = c(
    "F1-score",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "f1_score",
    "method",
    "downsampled_celltypes"
  )
)

# Get the ANOVA KNN classification F1 score for type after taking into account
# method and downsampled celltype, and save
knn_aov_results <- anova_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch hierarchical balanced",
  last_covariate = "type",
  metric = "f1_score"
)
fwrite(
  knn_aov_results,
  paste0(
    "outs/control/results/",
    "06_pbmc_hierarchical_knn_cell_ds_cell_aov_results_ctrl_method_ds_celltype.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#############################################################################
### REDO ALL STATISTICAL TESTS ABOVE BUT WITHOUT LIGER ###
#############################################################################

rm(list = ls())
gc()


`%ni%` <- Negate(`%in%`)
setwd("results/control/")

# Set seed for any sampling done 
set.seed(42)

##### Statistical analysis of PBMC 2 batch balanced data - baseline Fig 2 #####

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_files,
  value = TRUE
)
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)
clus_concat <- clus_concat[clus_concat$Method != "liger"]

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
clus_concord_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_concord_files,
  value = TRUE
)
clus_concord_loaded <- lapply(clus_concord_files, fread)
clus_concord_concat <- Reduce(rbind, clus_concord_loaded)
clus_concord_concat <- clus_concord_concat[
  clus_concord_concat$`Method 1` != "liger"
] 
clus_concord_concat <- clus_concord_concat[
  clus_concord_concat$`Method 2` != "liger"
]

# Load in and concatenate dge concordance summaries
setwd("../dge_concord_stats/")
dge_files <- list.files()
dge_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_files,
  value = TRUE
)
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)
dge_concat <- dge_concat[dge_concat$`Method 1` != "liger"]
dge_concat <- dge_concat[dge_concat$`Method 2` != "liger"]

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_files <- grep(
  "pbmc_2_batch_base_balanced",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
knn_concat <- knn_concat[knn_concat$Method != "liger"]

# Change to top level dir 
setwd("../../..")

### Statistical tests for ARI/AMI/etc. ### 

# Merge imbalance and clustering summary results
imba_clus_merged <- merge(
  clus_concat,
  imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_clus_unique_celltypes <- unique(imba_clus_merged$`Downsampled celltypes`)
imba_clus_unique_celltypes <- imba_clus_unique_celltypes[
  imba_clus_unique_celltypes %ni% "None"
]
imba_clus_none_celltype_len <- length(
  which(
    imba_clus_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_clus_none_celltype_draw <- sample(
  imba_clus_unique_celltypes,
  imba_clus_none_celltype_len,
  replace = TRUE
)
imba_clus_merged$`Downsampled celltypes`[
  which(
    imba_clus_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_clus_none_celltype_draw

# Format celltype names 
imba_clus_merged$`Downsampled celltypes` <- plyr::mapvalues(
  imba_clus_merged$`Downsampled celltypes`,
  from = c(
    "Monocyte_CD14",
    "Monocyte_FCGR3A",
    "CD4 T cell",
    "CD8 T cell"
  ),
  to = c(
    "CD14+ Monocyte",
    "FCGR3A+ Monocyte",
    "CD4+ T cell",
    "CD8+ T cell"
  )
)

# Format column names
colnames(imba_clus_merged) <- plyr::mapvalues(
  colnames(imba_clus_merged),
  from = c(
    "Celltype ARI Imbalanced",
    "Celltype AMI Imbalanced",
    "Celltype Homogeneity Imbalanced",
    "Celltype Completeness Imbalanced"
  ),
  to = c(
    "Celltype ARI",
    "Celltype AMI",
    "Celltype Homogeneity",
    "Celltype Completeness"
  )
)

# Indicate which samples are controls and which are real runs
imba_clus_merged$type <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
) 

# Create a function to do a wilcoxon rank sum test across the given
# subsets for the indicated metric 
wilcox_metric_test <- function(
    dataset, 
    dataset_name,
    subset_col,
    subset_1, 
    subset_2, 
    metric
) {
  subset_1_sub <- dataset[dataset[[subset_col]] %in% subset_1][[metric]]
  subset_2_sub <- dataset[dataset[[subset_col]] %in% subset_2][[metric]]
  wrs_p_value <- wilcox.test(subset_1_sub, subset_2_sub)$p.value
  wrs_summary_df <- data.frame(
    "Dataset" = dataset_name,
    "Subset_name" = subset_col,
    "Subset_1" = subset_1,
    "Subset_2" = subset_2,
    "Metric" = metric,
    "Wilcoxon_Rank_Sum_p" = wrs_p_value,
    "Samples_subset_1" = length(subset_1_sub),
    "Samples_subset_2" = length(subset_2_sub)
  )
  return(wrs_summary_df)
}

# Do a wilcoxin rank sum test between all control, ablated, and downsampled
# experiment results for all metrics 
metrics <- list(
  "Celltype ARI",
  "Celltype AMI",
  "Celltype Homogeneity",
  "Celltype Completeness",
  "Batch ARI",
  "Batch AMI",
  "Batch Homogeneity",
  "Batch Completeness"
)
comparison_subsets <- list(
  c("Control", "Downsampled"),
  c("Control", "Ablated")
)

# Get all results for batch and celltype metrics (globally) for the two given
# comparisons 
global_metric_results <- lapply(metrics, function(x) {
  lapply(comparison_subsets, function(y) {
    wrs_summary_df <- wilcox_metric_test(
      dataset = imba_clus_merged,
      dataset_name = "PBMC 2 batch base balanced",
      subset_col = "type",
      subset_1 = y[[1]],
      subset_2 = y[[2]],
      metric = x
    )
  })
})

# Concatenate all results for type subset
global_metric_results_concat_1 <- lapply(
  global_metric_results,
  function(x) {
    x_red <- Reduce(rbind, x)
    return(x_red)
  }
)
global_metric_results_concat_2 <- Reduce(rbind, global_metric_results_concat_1)

# Reperform for pooled (downsampled, ablated) comparison

# Create new pooled type column
imba_clus_merged$type_pooled <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  "Pooled downsampled/ablated"
) 

# Get all results for batch and celltype metrics (globally) for the given
# comparison
global_metric_results_pooled <- lapply(metrics, function(x) {
  wrs_summary_df <- wilcox_metric_test(
    dataset = imba_clus_merged,
    dataset_name = "PBMC 2 batch base balanced",
    subset_col = "type_pooled",
    subset_1 = "Control",
    subset_2 = "Pooled downsampled/ablated",
    metric = x
  )
})

# Concatenate all results for type_pooled subset
global_metric_results_pooled_concat <- Reduce(
  rbind,
  global_metric_results_pooled
)

# Concatenate individual and pooled subset results and save
global_metric_results_indi_pooled <- rbind(
  global_metric_results_concat_2,
  global_metric_results_pooled_concat
)
fwrite(
  global_metric_results_indi_pooled,
  "outs/control/results/06_pbmc_base_global_metric_wilcoxon_tests_no_liger.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Create function to perform ANOVA test, given either type or type pooled,
# a metric, and taking into account method and specific downsampled 
# celltype 
anova_metric_test <- function(
    dataset, 
    dataset_name,
    metric,
    last_covariate
){
  model_fit <- lm(
    as.formula(
      paste0(
        metric, 
        "~",
        "method+",
        "downsampled_celltypes+",
        last_covariate
      )
    ),
    data = dataset
  )
  anova_result <- anova(model_fit, test = "F")
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$metric <- metric
  anova_result_dt$last_covariate <- last_covariate
  return(anova_result_dt)
}

# Rename columns for use in formulas (no spaces)
colnames(imba_clus_merged) <- plyr::mapvalues(
  colnames(imba_clus_merged),
  from = c(
    "Celltype ARI",
    "Celltype AMI",
    "Celltype Homogeneity",
    "Celltype Completeness",
    "Batch ARI",
    "Batch AMI",
    "Batch Homogeneity",
    "Batch Completeness",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "celltype_ari",
    "celltype_ami",
    "celltype_homogeneity",
    "celltype_completeness",
    "batch_ari",
    "batch_ami",
    "batch_homogeneity",
    "batch_completeness",
    "method",
    "downsampled_celltypes"
  )
)

# Iterate over metrics and the important last covariate (type)
metrics_aov <- list(
  "celltype_ari",
  "celltype_ami",
  "celltype_homogeneity",
  "celltype_completeness",
  "batch_ari",
  "batch_ami",
  "batch_homogeneity",
  "batch_completeness"
)
metric_aov_results <- lapply(metrics_aov, function(x) {
  aov_summary_df <- anova_metric_test(
    dataset = imba_clus_merged,
    dataset_name = "PBMC 2 batch base balanced",
    metric = x,
    last_covariate = "type"
  )
})

# Concatenate results, perform FDR correction, and save 
metric_aov_results_concat <- Reduce(rbind, metric_aov_results)
metric_aov_results_concat$`FDR_q` <- p.adjust(
  metric_aov_results_concat$`Pr(>F)`
)
fwrite(
  metric_aov_results_concat,
  "outs/control/results/06_pbmc_base_metric_aov_results_types_ct_method_ctrl_no_liger.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

### Statistical tests for KNN classification metrics for baseline data ###

# Merge imbalance and knn classification results together
imba_knn_merged <- merge(
  imba_concat,
  knn_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_knn_merged <- distinct(imba_knn_merged)

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_knn_unique_celltypes <- unique(imba_knn_merged$`Downsampled celltypes`)
imba_knn_unique_celltypes <- imba_knn_unique_celltypes[
  imba_knn_unique_celltypes %ni% "None"
]
imba_knn_none_celltype_len <- length(
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_knn_none_celltype_draw <- sample(
  imba_knn_unique_celltypes,
  imba_knn_none_celltype_len,
  replace = TRUE
)
imba_knn_merged$`Downsampled celltypes`[
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_knn_none_celltype_draw

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype$type <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Create another column for pooling the ablated and downsampled data together
imba_knn_merged_celltype$type_pooled <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  "Pooled downsampled/ablated"
)

# Create a function to do a wilcoxon rank sum test across the given
# subsets for the KNN classification F1 scores where the F1 score 
# is for the celltype that was downsampled
wilcox_knn_test <- function(
    dataset, 
    dataset_name,
    subset_col,
    subset_1, 
    subset_2, 
    metric = "F1-score"
) {
  subset_1_sub <- dataset[dataset[[subset_col]] %in% subset_1][[metric]]
  subset_2_sub <- dataset[dataset[[subset_col]] %in% subset_2][[metric]]
  wrs_p_value <- wilcox.test(subset_1_sub, subset_2_sub)$p.value
  wrs_summary_df <- data.frame(
    "Dataset" = dataset_name,
    "Subset_name" = subset_col,
    "Subset_1" = subset_1,
    "Subset_2" = subset_2,
    "Metric" = metric,
    "Wilcoxon_Rank_Sum_p" = wrs_p_value,
    "Samples_subset_1" = length(subset_1_sub),
    "Samples_subset_2" = length(subset_2_sub)
  )
  return(wrs_summary_df)
}

# Get global KNN classification wilcoxon result across type subsets
# and concatenate
global_celltype_knn_results_type <- lapply(comparison_subsets, function(x) {
  wrs_summary_df <- wilcox_knn_test(
    dataset = imba_knn_merged_celltype,
    dataset_name = "PBMC 2 batch base balanced",
    subset_col = "type",
    subset_1 = x[[1]],
    subset_2 = x[[2]],
    metric = "F1-score"
  )
})
global_celltype_knn_results_type_concat <- Reduce(
  rbind,
  global_celltype_knn_results_type
)

# Get global KNN classification wilcoxon result across type_pooled subsets
# and concatenate 
global_celltype_knn_results_type_pooled <- wilcox_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch base balanced",
  subset_col = "type_pooled",
  subset_1 = "Control",
  subset_2 = "Pooled downsampled/ablated",
  metric = "F1-score"
)

# Concatenate type and type pooled results together and save
global_celltype_knn_results_type_and_pooled <- rbind(
  global_celltype_knn_results_type_concat,
  global_celltype_knn_results_type_pooled
)
fwrite(
  global_celltype_knn_results_type_and_pooled,
  "outs/control/results/06_pbmc_base_knn_cell_ds_cell_test_wilcoxon_tests_no_liger.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Create function to perform ANOVA test, given either type or type pooled,
# and taking into account method used and downsampled celltype (same as 
# the examined celltype)
anova_knn_test <- function(
    dataset, 
    dataset_name,
    last_covariate,
    metric = "f1_score"
){
  model_fit <- lm(
    as.formula(
      paste0(
        metric, 
        "~",
        "method+",
        "downsampled_celltypes+",
        last_covariate
      )
    ),
    data = dataset
  )
  anova_result <- anova(model_fit, test = "F")
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$metric <- metric
  anova_result_dt$last_covariate <- last_covariate
  return(anova_result_dt)
}

# Rename columns for use in formulas 
colnames(imba_knn_merged_celltype) <- plyr::mapvalues(
  colnames(imba_knn_merged_celltype),
  from = c(
    "F1-score",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "f1_score",
    "method",
    "downsampled_celltypes"
  )
)

# Get the ANOVA KNN classification F1 score for type after taking into account
# method and downsampled celltype, and save
knn_aov_results <- anova_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch base balanced",
  last_covariate = "type",
  metric = "f1_score"
)
fwrite(
  knn_aov_results,
  paste0(
    "outs/control/results/",
    "06_pbmc_base_knn_cell_ds_cell_aov_results_ctrl_method_celltype_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

##### Statistical analysis of PBMC 2 batch balanced - hierarchical Fig 2 #####

# Clear all previous files and reload KNN functions 
rm(list = ls())
gc()

wilcox_knn_test <- function(
    dataset, 
    dataset_name,
    subset_col,
    subset_1, 
    subset_2, 
    metric = "F1-score"
) {
  subset_1_sub <- dataset[dataset[[subset_col]] %in% subset_1][[metric]]
  subset_2_sub <- dataset[dataset[[subset_col]] %in% subset_2][[metric]]
  wrs_p_value <- wilcox.test(subset_1_sub, subset_2_sub)$p.value
  wrs_summary_df <- data.frame(
    "Dataset" = dataset_name,
    "Subset_name" = subset_col,
    "Subset_1" = subset_1,
    "Subset_2" = subset_2,
    "Metric" = metric,
    "Wilcoxon_Rank_Sum_p" = wrs_p_value,
    "Samples_subset_1" = length(subset_1_sub),
    "Samples_subset_2" = length(subset_2_sub)
  )
  return(wrs_summary_df)
}

anova_knn_test <- function(
    dataset, 
    dataset_name,
    last_covariate,
    metric = "f1_score"
){
  model_fit <- lm(
    as.formula(
      paste0(
        metric, 
        "~",
        "method+",
        "downsampled_celltypes+",
        last_covariate
      )
    ),
    data = dataset
  )
  anova_result <- anova(model_fit, test = "F")
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$metric <- metric
  anova_result_dt$last_covariate <- last_covariate
  return(anova_result_dt)
}

# Reload helper functions and seed
`%ni%` <- Negate(`%in%`)

# Set seed for any sampling done 
set.seed(42)

# Load in and concatenate summary files - just imbalance and KNN results
setwd("results/control/imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
knn_concat <- knn_concat[knn_concat$Method != "liger"]

# Change to top level dir
setwd("../../..")

# Merge imbalance and knn classification results together
imba_knn_merged <- merge(
  imba_concat,
  knn_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_knn_merged <- distinct(imba_knn_merged)

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_knn_unique_celltypes <- unique(imba_knn_merged$`Downsampled celltypes`)
imba_knn_unique_celltypes <- imba_knn_unique_celltypes[
  imba_knn_unique_celltypes %ni% "None"
]
imba_knn_none_celltype_len <- length(
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_knn_none_celltype_draw <- sample(
  imba_knn_unique_celltypes,
  imba_knn_none_celltype_len,
  replace = TRUE
)
imba_knn_merged$`Downsampled celltypes`[
  which(
    imba_knn_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_knn_none_celltype_draw

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype$type <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Create another column for pooling the ablated and downsampled data together
imba_knn_merged_celltype$type_pooled <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  "Pooled downsampled/ablated"
)

# Get global KNN classification wilcoxon result across type subsets
# and concatenate
comparison_subsets <- list(
  c("Control", "Downsampled"),
  c("Control", "Ablated")
)

global_celltype_knn_results_type <- lapply(comparison_subsets, function(x) {
  wrs_summary_df <- wilcox_knn_test(
    dataset = imba_knn_merged_celltype,
    dataset_name = "PBMC 2 batch hierarchical balanced",
    subset_col = "type",
    subset_1 = x[[1]],
    subset_2 = x[[2]],
    metric = "F1-score"
  )
})
global_celltype_knn_results_type_concat <- Reduce(
  rbind,
  global_celltype_knn_results_type
)

# Get global KNN classification wilcoxon result across type_pooled subsets
# and concatenate 
global_celltype_knn_results_type_pooled <- wilcox_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch hierarchical balanced",
  subset_col = "type_pooled",
  subset_1 = "Control",
  subset_2 = "Pooled downsampled/ablated",
  metric = "F1-score"
)

# Concatenate type and type pooled results together and save
global_celltype_knn_results_type_and_pooled <- rbind(
  global_celltype_knn_results_type_concat,
  global_celltype_knn_results_type_pooled
)
fwrite(
  global_celltype_knn_results_type_and_pooled,
  paste0(
    "outs/control/results/",
    "06_pbmc_hierarchical_knn_cell_ds_cell_test_wilcoxon_tests_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Rename columns for use in formulas 
colnames(imba_knn_merged_celltype) <- plyr::mapvalues(
  colnames(imba_knn_merged_celltype),
  from = c(
    "F1-score",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "f1_score",
    "method",
    "downsampled_celltypes"
  )
)

# Get the ANOVA KNN classification F1 score for type after taking into account
# method and downsampled celltype, and save
knn_aov_results <- anova_knn_test(
  dataset = imba_knn_merged_celltype,
  dataset_name = "PBMC 2 batch hierarchical balanced",
  last_covariate = "type",
  metric = "f1_score"
)
fwrite(
  knn_aov_results,
  paste0(
    "outs/control/results/",
    "06_pbmc_hierarchical_knn_cell_ds_cell_aov_results_ctrl_method_ds_celltype_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

#### Comparison of ANOVA KNN classification F1 scores for baseline and #### 
#### hierarchical setups - not including LIGER ####

# Rename hierarchical AOV results and reload in baseline
knn_aov_results_hierarchical <- knn_aov_results
knn_aov_results_baseline <- fread(
  paste0(
    "outs/control/results/",
    "06_pbmc_base_knn_cell_ds_cell_aov_results_ctrl_method_celltype_no_liger.tsv"
  )
)

# Plot the F-values for the three covariates considered in the ANOVA test 
f_vals_baseline <- knn_aov_results_baseline$`F value`[1:3]
f_vals_hierarchical <- knn_aov_results_hierarchical$`F value`[1:3]
covars_baseline <- knn_aov_results_baseline$Covariate[1:3]
covars_hierarchical <- knn_aov_results_hierarchical$Covariate[1:3]

knn_aov_comp_df <- data.frame(
  "Covariates" = c(covars_baseline, covars_hierarchical),
  "F_values" = c(f_vals_baseline, f_vals_hierarchical),
  "Subset" = c(
    "Baseline", 
    "Baseline", 
    "Baseline", 
    "Hierarchical", 
    "Hierarchical", 
    "Hierarchical"
  )
)

knn_aov_comp_df_melted <- reshape2::melt(
  knn_aov_comp_df,
  id.vars = c("Subset", "Covariates"),
  measure.vars = "F_values"
)
knn_aov_comp_df_melted$Covariates <- plyr::mapvalues(
  knn_aov_comp_df_melted$Covariates,
  from = c(
    "type",
    "method",
    "downsampled_celltypes"
  ),
  to = c(
    "Balanced vs perturbed",
    "Integration method",
    "Celltype downsampled"
  )
) 

knn_aov_comp_df_melted$Subset <- plyr::mapvalues(
  x = knn_aov_comp_df_melted$Subset,
  from = c("Baseline", "Hierarchical"),
  to = c(
    "Baseline balanced \nPBMC data", 
    "Hierarchically clustered \nPBMC data"
  )
)

ggplot(data = knn_aov_comp_df_melted, aes(Covariates, value)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(),
    aes(
      fill = Subset
    )
  ) + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic() + 
  coord_flip () +
  labs(x = "Covariate", y = "ANOVA F-statistic") +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "06_pbmc_base_vs_hierarchical_knn_aov_results_comp_no_liger.pdf"
  ),
  width = 12,
  height = 8,
  device = cairo_pdf
)
