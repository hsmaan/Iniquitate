library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(dotwhisker)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Change to results dir for control data 
setwd("../../../results/control/")

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

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
clus_concord_loaded <- lapply(clus_concord_files, fread)
clus_concord_concat <- Reduce(rbind, clus_concord_loaded)

# Load in and concatenate dge concordance summaries
setwd("../dge_concord_stats/")
dge_files <- list.files()
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)

# Change to top level dir 
setwd("../../..")

### Celltype-specific clustering index and cluster number 
### results for integration pre and post downsampling

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

# Subset data by just results for Seurat
imba_clus_merged_seurat <- imba_clus_merged[
  imba_clus_merged$Method == "seurat"
]

# Plot results for Celltype ARI based on downsampling, ablation and control,
# subset by celltype
ggplot(data = imba_clus_merged_seurat, aes(
  x = factor(`Downsampled celltypes`),
  y = `Celltype ARI Imbalanced`,
  color = `type`
)) +
  geom_jitter() +
  theme_few() +
  labs(
    color = "Datapoint type",
    title = paste0(
      "PBMC 2 Batch Balanced Downsampling/Ablation (One Batch Randomized) ",
      "- Seurat results"
    ),
    x = "Celltype downsampled",
    y = "Celltype Adjusted Rand Index Post-integration"
  ) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))
ggsave(
  "outs/control/figures/04_pbmc_ds_ablate_seurat_celltype_ari.pdf",
  width = 14,
  height = 7
)

# Plot results for Batch ARI based on downsampled, ablation and control,
# subset by celltype 
ggplot(data = imba_clus_merged_seurat, aes(
  x = factor(`Downsampled celltypes`),
  y = `Batch ARI`,
  color = `type`
)) +
  geom_jitter() +
  theme_few() +
  labs(
    color = "Datapoint type",
    title = paste0(
      "PBMC 2 Batch Balanced Downsampling/Ablation (One Batch Randomized) ",
      "- Seurat results"
    ),
    x = "Celltype downsampled",
    y = "(1 - Batch) Adjusted Rand Index Post-integration"
  ) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))
ggsave(
  "outs/control/figures/04_pbmc_ds_ablate_seurat_batch_ari.pdf",
  width = 14,
  height = 7
)

# Plot batch and celltype ARI values in the same panes side by side 
imba_clus_merged_seurat_melt <- melt(
  imba_clus_merged_seurat,
  id.vars = colnames(imba_clus_merged_seurat)[
    which(colnames(imba_clus_merged_seurat) %ni% c(
      "Celltype ARI Imbalanced",
      "Celltype AMI Imbalanced",
      "Celltype Completeness Imbalanced",
      "Celltype Homogeneity Imbalanced",
      "Celltype ARI Balanced",
      "Celltype AMI Balanced",
      "Celltype Completeness Balanced",
      "Celltype Homogeneity Balanced",
      "Batch ARI",
      "Batch AMI",
      "Batch Homogeneity",
      "Batch Completeness"
  ))]
)
colnames(imba_clus_merged_seurat_melt)[
  ncol(imba_clus_merged_seurat_melt)
] <- "Score"
colnames(imba_clus_merged_seurat_melt)[
  ncol(imba_clus_merged_seurat_melt) - 1
] <- "Metric"
imba_clus_merged_seurat_melt_ari_sub <- imba_clus_merged_seurat_melt[
  imba_clus_merged_seurat_melt$Metric %in% c(
    "Celltype ARI Imbalanced",
    "Batch ARI"
  ),
]   
imba_clus_merged_seurat_melt_ari_sub$Metric <- gsub(
  "Celltype ARI Imbalanced",
  "Celltype ARI",
  imba_clus_merged_seurat_melt_ari_sub$Metric, 
)

ggplot(data = imba_clus_merged_seurat_melt_ari_sub, aes(
  x = factor(`Downsampled celltypes`),
  y = `Score`,
  color = `type`
)) +
  facet_wrap(.~Metric) +
  geom_jitter() +
  theme_few() +
  labs(
    color = "Datapoint type",
    title = paste0(
      "PBMC 2 Batch Balanced Downsampling/Ablation (One Batch Randomized) ",
      "- Seurat results"
    ),
    x = "Celltype downsampled",
    y = "(1 - Batch) Adjusted Rand Index Post-integration"
  ) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.6)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))
ggsave(
  "outs/control/figures/04_pbmc_ds_ablate_seurat_batch_celltype_ari.pdf",
  width = 14,
  height = 7
)

# Determine the change in cluster number based on downsampling and/or 
# ablation of different celltypes
ggplot(data = imba_clus_merged_seurat, aes(
  x = factor(`Downsampled celltypes`),
  y = `Cluster number`,
  color = `type`
)) +
  geom_jitter() +
  theme_few() +
  labs(
    color = "Datapoint type",
    title = paste0(
      "PBMC 2 Batch Balanced Downsampling/Ablation (One Batch Randomized) ",
      "- Seurat results"
    ),
    x = "Celltype downsampled",
    y = "Number of Leiden clusters post integration"
  ) +
  theme(axis.title.x = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16)) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))
ggsave(
  "outs/control/figures/04_pbmc_ds_ablate_seurat_num_clusters.pdf",
  width = 14,
  height = 7
)

### Celltype-specific KNN classification results for different celltypes 
### integration pre and post downsampling

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

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Subset for only the Seurat method 
imba_knn_merged_celltype_seurat <- imba_knn_merged_celltype[
  imba_knn_merged_celltype$Method %in% c("seurat")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype_seurat$type <- ifelse(
  imba_knn_merged_celltype_seurat$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype_seurat$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Plot results for KNN classification score by celltype 
ggplot(data = imba_knn_merged_celltype_seurat, aes(
  x = factor(
    `type`,
    levels = c(
      "Control",
      "Downsampled",
      "Ablated"
    )
  ),
  y = `F1-score`,
  color = factor(
    `type`,
    levels = c(
      "Control",
      "Downsampled",
      "Ablated"
    )
  ),
)) +
  facet_wrap(.~Celltype) +
  geom_jitter() +
  theme_few() +
  labs(
    color = "Datapoint type",
    title = paste0(
      "PBMC 2 Batch Balanced Downsampling/Ablation (One Batch Randomized) ",
      "- Seurat results"
    ),
    x = "",
    y = "KNN classification F1-score - post integration"
  ) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_text(size = 14)) +
  theme(strip.text.x = element_text(size = 14)) +
  theme(plot.title = element_text(size = 16)) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  theme(axis.ticks.x = element_blank())
ggsave(
    "outs/control/figures/04_pbmc_ds_ablate_seurat_knn_class_score.pdf",
    width = 12,
    height = 7
)
