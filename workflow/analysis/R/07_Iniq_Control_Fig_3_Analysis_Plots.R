library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(ggpubr)
library(dotwhisker)
library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)
library(le)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Change to results dir for control data 
setwd("../../../results/control/")

##### Analysis of PBMC 2 batch balanced data - baseline #####

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
gc()

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
gc()

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
gc()

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
gc()

# Load in and concatenate dge ranking summaries
setwd("../dge_ranking_stats/")
dge_rank_files <- list.files()
dge_rank_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_rank_files,
  value = TRUE
)
dge_rank_loaded <- lapply(dge_rank_files, fread)
dge_rank_concat <- Reduce(rbind, dge_rank_loaded)
gc()

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
gc()

# Change to top level dir 
setwd("../../..")

### Fig 3A) - Variability in number of clusters after downsampling the control
### dataset using different celltypes 

# Merge together imbalance and cluster number results
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

# Indicate which panels are control and which ones are ablations or downsampling
imba_clus_merged$type <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
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

# Create plot of effects on cluster number, facetted by method, filled by
# type (control, ablation, downsampling), and plotted based on downsampled
# celltype 
ggplot(data = imba_clus_merged, aes(
    x = `Downsampled celltypes`,
    y = `Cluster number`
  )
) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
) + 
  geom_jitter(color = "black", size = 0.4, alpha = 0.8) +
  facet_wrap(.~Method, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "Number of Leiden clusters post-integration"
  ) +
  coord_flip() +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  "outs/control/figures/07_pbmc_ds_ablate_allmethod_cluster_number.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)  

### Fig 3B) - Correlation between the number of clusters per method and the 
### adjusted rand index for celltype and batch ARI, on a per-method level
ggscatter(imba_clus_merged, 
          x = "Celltype ARI Imbalanced", 
          y = "Cluster number", 
          size = 0.4,
          combine = TRUE,
          xlab = "Celltype ARI post-integration",
          ylab = "Number of Leiden clusters post-integration",
          palette = "jco",
          add = "reg.line", 
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE) +
  facet_wrap(.~`Method`, scales = "fixed") +
  stat_cor(
    method = "spearman", 
    label.x = 0.1, 
    label.y = 0.1, 
    p.digits = 2,
    size = 5
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  "outs/control/figures/07_pbmc_ds_ablate_clus_num_celltype_ari_corr.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

ggscatter(imba_clus_merged, 
          x = "Batch ARI", 
          y = "Cluster number", 
          size = 0.4,
          combine = TRUE,
          xlab = "Batch ARI post-integration",
          ylab = "Number of Leiden clusters post-integration",
          palette = "jco",
          add = "reg.line", 
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE) +
  facet_wrap(.~`Method`, scales = "fixed") +
  stat_cor(
    method = "spearman", 
    label.x = 0.9, 
    label.y = 0.1, 
    p.digits = 2,
    size = 5
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  "outs/control/figures/07_pbmc_ds_ablate_clus_num_batch_ari_corr.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

### Fig 3C) - concordance of DGE 
