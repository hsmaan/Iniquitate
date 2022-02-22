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
setwd("../dge_concord_summaries/")
dge_files <- list.files()
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)

# Change to top level dir 
setwd("../../..")

### Celltype-specific results for integration pre and post downsampling

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
      "PBMC 2 Batch Balanced Downsampling/Ablation Ablation (One Batch) ",
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
