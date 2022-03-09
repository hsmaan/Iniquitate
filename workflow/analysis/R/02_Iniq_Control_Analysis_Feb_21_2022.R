library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(dotwhisker)

# Change to results dir for control data 
setwd("../../../results/control/") 

# Helper functions
`%ni%` <- Negate(`%in%`)

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
setwd("../dge_concord_stats//")
dge_files <- list.files()
dge_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_files,
  value = TRUE
)
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)

# Change to top level dir 
setwd("../../..")

# Create output directory for figures if not present
if (!dir.exists("outs/control/figures")) {
  dir.create("outs/control/figures", recursive = TRUE)
}

### Analysis of clustering results based on method 

# Merge clustering data with imbalance data 
clus_concat_imba <- merge(
  clus_concat, 
  imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Subset data for comparison of 1 celltype downsampled to none, full ablation
clus_concat_imba_sub_ablation <- clus_concat_imba[
  clus_concat_imba$`Proportion downsampled` %in% c(0, 1)
]
clus_concat_imba_sub_ablation$`Ablation` <- ifelse(
  clus_concat_imba_sub_ablation$`Number of batches downsampled` %in% 1 &
    clus_concat_imba_sub_ablation$`Number of celltypes downsampled` %in% 1 &
    clus_concat_imba_sub_ablation$`Proportion downsampled` %in% 0,
  "Yes",
  "No"
)

# Subset data for only those that had FCGR3A monocytes downsampled 
clus_concat_sub_ablation_mfcg <- clus_concat_imba_sub_ablation[
  clus_concat_imba_sub_ablation$`Downsampled celltypes` %in% c(
    "Monocyte_FCGR3A",
    "None"
  )
]

# Plot Celltype ARI values based on method (balanced) for ablation vs non
# ablation cases 
ggplot(data = clus_concat_sub_ablation_mfcg, aes(
  x = factor(`Ablation`),
  y = `Celltype ARI Imbalanced`,
  color = `Method`
)) +
  facet_wrap(.~Method) +
  ylim(c(0, 1)) +
  geom_jitter() +
  theme_few() +
  labs(
    title = paste0(
      "PBMC 2 Batch Balanced FCGR3A+ Monocyte Ablation ",
      "(One Batch Randomized)"
    ),
    x = "Ablation of FCGR3A+ Monocytes",
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
  "outs/control/figures/02_monocyte_ablation_results_celltype_ari.pdf",
  width = 12, 
  height = 8
)

# Plot Batch ARI values based on method (balanced) for ablation vs non
# ablation cases 
ggplot(data = clus_concat_sub_ablation_mfcg, aes(
  x = factor(`Ablation`),
  y = `Batch ARI`,
  color = `Method`
)) +
  facet_wrap(.~Method) +
  ylim(c(0, 1.05)) +
  geom_jitter() +
  theme_few() +
  labs(
    title = paste0(
      "PBMC 2 Batch Balanced FCGR3A+ Monocyte Ablation ",
      "(One Batch Randomized)"
    ),
    x = "Ablation of FCGR3A+ Monocytes",
    y = "Batch Adjusted Rand Index Post-integration"
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
  "outs/control/figures/02_monocyte_ablation_results_batch_ari.pdf",
  width = 12, 
  height = 8
)

# Plot the results for the number of clusters based on method for ablation vs
# non-ablation cases 
ggplot(data = clus_concat_sub_ablation_mfcg, aes(
  x = factor(`Ablation`),
  y = `Cluster number`,
  color = `Method`
)) +
  facet_wrap(.~Method) +
  geom_jitter() +
  theme_few() +
  labs(
    title = paste0(
      "PBMC 2 Batch Balanced FCGR3A+ Monocyte Ablation ",
      "(One Batch Randomized)"
    ),
    x = "Ablation of FCGR3A+ Monocytes",
    y = "Number of clusters post-integration"
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
  "outs/control/figures/02_monocyte_ablation_results_clus_number.pdf",
  width = 12, 
  height = 8
)

# Subset data for comparison of 1 celltype downsampled to none, 0.1 proportion
clus_concat_imba_sub_ds <- clus_concat_imba[
  clus_concat_imba$`Proportion downsampled` %in% c(1, 0.1)
]
clus_concat_imba_sub_ds$`Downsample` <- ifelse(
  clus_concat_imba_sub_ds$`Number of batches downsampled` %in% 1 &
    clus_concat_imba_sub_ds$`Number of celltypes downsampled` %in% 1 &
    clus_concat_imba_sub_ds$`Proportion downsampled` %in% 0.1,
  "Yes",
  "No"
)

# Subset data for only those that had FCGR3A monocytes downsampled 
clus_concat_imba_sub_ds_mfcg <- clus_concat_imba_sub_ds[
  clus_concat_imba_sub_ds$`Downsampled celltypes` %in% c(
    "Monocyte_FCGR3A",
    "None"
  )
]

# Plot Celltype ARI values based on method (balanced) for downsample vs non
# downsample cases 
ggplot(data = clus_concat_imba_sub_ds_mfcg, aes(
  x = factor(`Downsample`),
  y = `Celltype ARI Imbalanced`,
  color = `Method`
)) +
  facet_wrap(.~Method) +
  ylim(c(0, 1)) +
  geom_jitter() +
  theme_few() +
  labs(
    title = paste0(
      "PBMC 2 Batch Balanced FCGR3A+ Monocyte Downsampling ",
      "(One Batch Randomized)"
    ),
    x = "Downsampling (to 10% of cells) of FCGR3A+ Monocytes",
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
  "outs/control/figures/02_monocyte_downsample_results_celltype_ari.pdf",
  width = 12, 
  height = 8
)

# Plot Batch ARI values based on method (balanced) for downsample vs non
# downsample cases 
ggplot(data = clus_concat_imba_sub_ds_mfcg, aes(
  x = factor(`Downsample`),
  y = `Batch ARI`,
  color = `Method`
)) +
  facet_wrap(.~Method) +
  ylim(c(0, 1.05)) +
  geom_jitter() +
  theme_few() +
  labs(
    title = paste0(
      "PBMC 2 Batch Balanced FCGR3A+ Monocyte Downsampling ",
      "(One Batch Randomized)"
    ),
    x = "Downsampling (to 10% of cells) of FCGR3A+ Monocytes",
    y = "Batch Adjusted Rand Index Post-integration"
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
  "outs/control/figures/02_monocyte_downsample_results_batch_ari.pdf",
  width = 12, 
  height = 8
)

# Plot the results for the number of clusters based on method for downsampling 
# vs non-downsampling cases 
ggplot(data = clus_concat_imba_sub_ds_mfcg, aes(
  x = factor(`Downsample`),
  y = `Cluster number`,
  color = `Method`
)) +
  facet_wrap(.~Method) +
  geom_jitter() +
  theme_few() +
  labs(
    title = paste0(
      "PBMC 2 Batch Balanced FCGR3A+ Monocyte Downsampling ",
      "(One Batch Randomized)"
    ),
    x = "Downsampling (to 10% of cells) of FCGR3A+ Monocytes",
    y = "Number of clusters post-integration"
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
  "outs/control/figures/02_monocyte_downsample_results_clus_number.pdf",
  width = 12, 
  height = 8
)



