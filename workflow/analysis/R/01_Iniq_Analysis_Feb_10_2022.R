library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(dotwhisker)

# Change to results dir
setwd("../../../results/")

# Helper functions
`%ni%` <- Negate(`%in%`)

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)

# Change to top level dir 
setwd("..")

# Create output directory for figures if not present
if (!dir.exists("../outs/figures")) {
  dir.create("../outs/figures", recursive = TRUE)
}

### Figure 1 - broad patterns of effects of imbalance on integration

# Merge imbalance and clustering summary results
colnames(clus_concat)[2] <- "Number of batches downsampled"
clus_imba_merged <- merge(
  imba_concat,
  clus_concat,
  by = c(
    "Dataset",
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Replicate"
  ),
  allow.cartesian = TRUE
)

## Determine ANOVA p-values for batch and celltype ARI based on other 
# characteristics 
lm_celltype_ari <- lm(
  `Celltype ARI` ~ `Total batches` + `Cell number` + `Method` + 
  `Number of batches downsampled` + `Number of celltypes downsampled` + 
  `Proportion downsampled.x`, 
  data = clus_imba_merged
)

dwplot(lm_celltype_ari) +
  theme_bw() +
  labs(
    x = "Coefficient estimate for Celltype ARI",
    y = "Covariate"
  )
ggsave("../outs/figures/01_fig_1_anova_celltype_ari.pdf", height = 8, width = 6)

### Caveats here are that there are not enough (n = 0) points, and because these 
### methods are deterministic (most of them), they might not be independant runs
### Should I still run all of these iterations regardless? To add more power 
### to the case of n = 0 (0 celltypes downsampled?) - all others have 10 
### replicates


## Plot celltype and batch statistics by dataset and based on summary 
## statistics previously indicated 
dataset_list = list(
  "mouse_hindbrain_6_batch",
  "pbmc_2_batch",            
  "pbmc_4_batch",            
  "peng_pdac_8_batch"
)

# Keep only relevant columns and clean up - both batch and dataset
clus_imba_merged_celltype <- clus_imba_merged[,
  colnames(clus_imba_merged) %ni% c(
    "Proportion downsampled.x",
    "Batch ARI",
    "Batch AMI",
    "Batch Homogeneity",
    "Batch Completeness"
  ),
  with = FALSE
]
colnames(clus_imba_merged_celltype)[11] <- "Proportion downsampled"

clus_imba_merged_batch <- clus_imba_merged[,
  colnames(clus_imba_merged) %ni% c(
    "Proportion downsampled.x",
    "Celltype ARI",
    "Celltype AMI",
    "Celltype Homogeneity",
    "Celltype Completeness"
  ),
  with = FALSE
]
colnames(clus_imba_merged_batch)[11] <- "Proportion downsampled"

# Melt data based on celltype/batch ARI/AMI/etc
clus_imba_merged_celltype_melt <- reshape2::melt(
  clus_imba_merged_celltype,
  id = c(1:14)
)
colnames(clus_imba_merged_celltype_melt)[15] <- "Metric"
colnames(clus_imba_merged_celltype_melt)[16] <- "Value"

clus_imba_merged_batch_melt <- reshape2::melt(
  clus_imba_merged_batch,
  id = c(1:14)
)
colnames(clus_imba_merged_batch_melt)[15] <- "Metric"
colnames(clus_imba_merged_batch_melt)[16] <- "Value"

# Plot celltype results based on dataset using 
# celltype intersection ratio 
lapply(dataset_list, function(x) {
  # Subset based on dataset 
  clus_imba_celltype_subset <- clus_imba_merged_celltype_melt[
    clus_imba_merged_celltype_melt$Dataset %in% x,
  ]
  ggplot(data = clus_imba_celltype_subset, aes(
    x = factor(round(`Celltype intersection ratio`, 3)), 
    y = Value,
    groups = factor(round(`Celltype intersection ratio`, 3))
  )) + 
    facet_wrap(.~Metric, scales = "free_y") +
    geom_boxplot(fill = "dodgerblue2", alpha = 0.8, notch = TRUE) +
    theme_few() +
    labs(
      title = x,
      x = "Celltype intersection ratio (Jaccard distance)",
      y = ""
    )
  ggsave(
    paste0(
      "../outs/figures/", 
      "01_fig_1_celltype_int_ratio_celltype_stats_", 
      x, 
      ".pdf"
    ),
    height = 8, 
    width = 14
  )
})

# Plot celltype results based on dataset using 
# cosine proportion distance 
lapply(dataset_list, function(x) {
  # Subset based on dataset 
  clus_imba_celltype_subset <- clus_imba_merged_celltype_melt[
    clus_imba_merged_celltype_melt$Dataset %in% x,
  ]
  ggplot(data = clus_imba_celltype_subset, aes(
    x = factor(round(`Mean proportion cosine distance`, 3)), 
    y = Value,
    groups = factor(round(`Mean proportion cosine distance`, 3))
  )) + 
    facet_wrap(.~Metric, scales = "free_y") +
    geom_boxplot(fill = "dodgerblue2", alpha = 0.8, notch = TRUE) +
    theme_few() +
    labs(
      title = x,
      x = "Mean cosine proportion distance",
      y = ""
    )
  ggsave(
    paste0(
      "../outs/figures/", 
      "01_fig_1_mean_prop_cosine_distance_celltype_stats_", 
      x, 
      ".pdf"
    ),
    height = 8, 
    width = 16
  )
})

# Plot celltype results based on dataset using 
# cell number coefficient of var 
lapply(dataset_list, function(x) {
  # Subset based on dataset 
  clus_imba_celltype_subset <- clus_imba_merged_celltype_melt[
    clus_imba_merged_celltype_melt$Dataset %in% x,
  ]
  ggplot(data = clus_imba_celltype_subset, aes(
    x = factor(round(`Length coeff var`, 3)), 
    y = Value,
    groups = factor(round(`Length coeff var`, 3))
  )) + 
    facet_wrap(.~Metric, scales = "free_y") +
    geom_boxplot(fill = "dodgerblue2", alpha = 0.8, notch = TRUE) +
    theme_few() +
    labs(
      title = x,
      x = "Cell number coefficient of variation",
      y = ""
    )
  ggsave(
    paste0(
      "../outs/figures/", 
      "01_fig_1_cell_num_coeff_var_celltype_stats_", 
      x, 
      ".pdf"
    ),
    height = 8, 
    width = 16
  )
})

# Plot batch results based on dataset using
# celltype intersection ratio 
lapply(dataset_list, function(x) {
  # Subset based on dataset 
  clus_imba_batch_subset <- clus_imba_merged_batch_melt[
    clus_imba_merged_batch_melt$Dataset %in% x,
  ]
  ggplot(data = clus_imba_batch_subset, aes(
    x = factor(round(`Celltype intersection ratio`, 3)), 
    y = Value,
    groups = factor(round(`Celltype intersection ratio`, 3))
  )) + 
    facet_wrap(.~Metric, scales = "free_y") +
    geom_boxplot(fill = "dodgerblue2", alpha = 0.8, notch = TRUE) +
    theme_few() +
    labs(
      title = x,
      x = "Celltype intersection ratio (Jaccard distance)",
      y = ""
    )
  ggsave(
    paste0(
      "../outs/figures/", 
      "01_fig_1_celltype_int_ratio_batch_stats_", 
      x, 
      ".pdf"
    ),
    height = 8, 
    width = 14
  )
})

# Plot batch results based on dataset using
# mean proportion cosine distance distance 
lapply(dataset_list, function(x) {
  # Subset based on dataset 
  clus_imba_batch_subset <- clus_imba_merged_batch_melt[
    clus_imba_merged_batch_melt$Dataset %in% x,
  ]
  ggplot(data = clus_imba_batch_subset, aes(
    x = factor(round(`Mean proportion cosine distance`, 3)), 
    y = Value,
    groups = factor(round(`Mean proportion cosine distance`, 3))
  )) + 
    facet_wrap(.~Metric, scales = "free_y") +
    geom_boxplot(fill = "dodgerblue2", alpha = 0.8, notch = TRUE) +
    theme_few() +
    labs(
      title = x,
      x = "Mean cosine proportion distance",
      y = ""
    )
  ggsave(
    paste0(
      "../outs/figures/", 
      "01_fig_1_cosine_prop_distance_batch_stats_", 
      x, 
      ".pdf"
    ),
    height = 8, 
    width = 16
  )
})

### Figure 2 - broad patterns of effects of imbalance on integration

