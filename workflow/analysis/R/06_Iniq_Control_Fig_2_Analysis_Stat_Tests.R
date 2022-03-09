library(data.table)
library(tidyverse)
library(reshape2)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Change to results dir for control data 
setwd("../../../results/control/")

# Set seed for any sampling done 
set.seed(42)

###Statistical analysis of PBMC 2 batch balanced data - baseline ### 

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
    metric,
    samples_per_subset
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
    "Samples_subset_2" = length(subset)
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
      metric = x,
      samples_per_subset = 1200
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
    metric = x,
    samples_per_subset = 1200
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
  "outs/control/results/06_global_metric_wilcoxon_tests.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# # Create function to perform ANOVA test, given either type or type pooled,
# # a metric, and taking into account method and specific downsampled 
# # celltype 
# anova_metric_test <- function(
#   dataset, 
#   dataset_name,
#   metric,
#   total_samples,
# )
#   
# # Rename columns for use in formulas 
# colnames(imba_clus_merged) <- plyr::mapvalues(
#   colnames(imba_clus_merged),
#   from = c(
#     "Celltype ARI",
#     "Celltype AMI",
#     "Celltype Homogeneity",
#     "Celltype Completeness",
#     "Batch ARI",
#     "Batch AMI",
#     "Batch Homoge"
#   )
# )
# 
# # Subset for and rename relevant columns in imba_clus_merged
# imba_clus_merged_sub_cari <- imba_clus_merged[
#   ,c(
#     "Celltype ARI Imbalanced",
#     "Method",
#     "Downsampled celltypes",
#     "type"
#   )
# ]
# colnames(imba_clus_merged_sub_cari) <- c(
#   "celltype_ari_imbalanced",
#   "method",
#   "downsampled_celltypes",
#   "type"
# )
# 
# # Perform ANOVA test for celltype ari and type, taking into account the 
# # celltype that was downsampled and the method used - save result
# celltype_ari_model <- lm(
#   as.formula(
#     paste0(
#       "celltype_ari_imbalanced", 
#       "~",
#       "method+",
#       "downsampled_celltypes+",
#       "type"
#     )
#   ),
#   data = imba_clus_merged_sub_cari
# )
# cari_anova_result <- anova(celltype_ari_model, test = "F")
# fwrite(
#   cari_anova_result,
#   "outs/control/results/06_celltype_ari_global_anova.tsv",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )
# 
# # Perform the same as above for batch ARI through an ANOVA test, after taking 
# # into account the method and celltype that was downsampled 
# 
# # Subset for and rename relevant columns in imba_clus_merged
# imba_clus_merged_sub_bari <- imba_clus_merged[
#   ,c(
#     "Batch ARI",
#     "Method",
#     "Downsampled celltypes",
#     "type"
#   )
# ]
# colnames(imba_clus_merged_sub_bari) <- c(
#   "batch_ari_imbalanced",
#   "method",
#   "downsampled_celltypes",
#   "type"
# )
# 
# # Perform ANOVA test for batch ari and type, taking into account the 
# # celltype that was downsampled and the method used - save result
# batch_ari_model <- lm(
#   as.formula(
#     paste0(
#       "batch_ari_imbalanced", 
#       "~",
#       "method+",
#       "downsampled_celltypes+",
#       "type"
#     )
#   ),
#   data = imba_clus_merged_sub_bari
# )
# bari_anova_result <- anova(batch_ari_model, test = "F")
# fwrite(
#   bari_anova_result,
#   "outs/control/results/06_batch_ari_global_anova.tsv",
#   sep = "\t",
#   quote = FALSE,
#   row.names = FALSE,
#   col.names = TRUE
# )
# 
# ### Statistical tests for KNN classification results ###
# 
# 
