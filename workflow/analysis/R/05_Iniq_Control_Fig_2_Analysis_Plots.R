library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(dotwhisker)
library(Seurat)
library(SeuratDisk)

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

# Change to PBMC data dir  
setwd("resources/h5ad_files/int_datasets/pbmc_2_batch_balanced/")

# Convert balanced pbmc data files to h5seurat format
Convert(
  "tran_exp5_pbmc_batch1_balanced.h5ad", 
  dest = "h5seurat", 
  overwrite = TRUE
)
Convert(
  "tran_exp5_pbmc_batch2_balanced.h5ad",
  dest = "h5seurat",
  overwrite = TRUE
)

# Load h5ad files for both balanced pbmc datasets 
pbmc_1 <- SeuratDisk::LoadH5Seurat("tran_exp5_pbmc_batch1_balanced.h5seurat")
pbmc_2 <- SeuratDisk::LoadH5Seurat("tran_exp5_pbmc_batch2_balanced.h5seurat")

# Change to top level dir 
setwd("../../../../")

# Process both datasets independantly and together
pbmc_combined <- merge(pbmc_1, pbmc_2)

pbmc_1 <- NormalizeData(pbmc_1)
pbmc_1 <- FindVariableFeatures(
  pbmc_1, selection.method = "vst", nfeatures = 2000
)
all.genes <- rownames(pbmc_1)
pbmc_1 <- ScaleData(pbmc_1, features = all.genes)
pbmc_1 <- RunPCA(pbmc_1, features = VariableFeatures(object = pbmc_1))
pbmc_1 <- FindNeighbors(pbmc_1, dims = 1:20)
pbmc_1 <- FindClusters(pbmc_1, resolution = 0.5)
pbmc_1 <- RunUMAP(pbmc_1, dims = 1:20)

pbmc_2 <- NormalizeData(pbmc_2)
pbmc_2 <- FindVariableFeatures(
  pbmc_2, selection.method = "vst", nfeatures = 2000
)
all.genes <- rownames(pbmc_2)
pbmc_2 <- ScaleData(pbmc_2, features = all.genes)
pbmc_2 <- RunPCA(pbmc_2, features = VariableFeatures(object = pbmc_2))
pbmc_2 <- FindNeighbors(pbmc_2, dims = 1:20)
pbmc_2 <- FindClusters(pbmc_2, resolution = 0.5)
pbmc_2 <- RunUMAP(pbmc_2, dims = 1:20)

pbmc_combined <- NormalizeData(pbmc_combined)
pbmc_combined <- FindVariableFeatures(
  pbmc_combined, selection.method = "vst", nfeatures = 2000
)
all.genes <- rownames(pbmc_combined)
pbmc_combined <- ScaleData(pbmc_combined, features = all.genes)
pbmc_combined <- RunPCA(pbmc_combined, features = VariableFeatures(
  object = pbmc_combined
))
pbmc_combined <- FindNeighbors(pbmc_combined, dims = 1:20)
pbmc_combined <- FindClusters(pbmc_combined, resolution = 0.5)
pbmc_combined <- RunUMAP(pbmc_combined, dims = 1:20)

### Create plots of celltype embeddings for the two datasets combined

# First format the metadata a bit
pbmc_combined@meta.data$celltype <- plyr::mapvalues(
  pbmc_combined@meta.data$celltype,
  from = c(
    "Monocyte_CD14",
    "Monocyte_FCGR3A"
  ),
  to = c(
    "CD14+ Monocyte",
    "FCGR3A+ Monocyte"
  )
)
pbmc_combined@meta.data$batch <- plyr::mapvalues(
  pbmc_combined@meta.data$batch, 
  from = c(
    "batch_1",
    "batch_2"
  ),
  to = c(
    "Batch 1",
    "Batch 2"
  )
)

# Plot combined data results as they are, for both celltype and batch 
DimPlot(pbmc_combined, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  "outs/control/figures/05_pbmc_balanced_combined_celltypes.pdf",
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  "outs/control/figures/05_pbmc_balanced_combined_batch.pdf",
  width = 6,
  height = 6
)

# Remove (ablate) CD14+ Monocytes from batch 2 and replot figures 
pbmc_combined_cd14_ablate <- pbmc_combined[, 
  -(which(pbmc_combined@meta.data$celltype %in% ("CD14+ Monocyte") &
    pbmc_combined@meta.data$batch %in% ("Batch 2")))
]
DimPlot(pbmc_combined_cd14_ablate, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  "outs/control/figures/05_pbmc_balanced_combined_celltypes_cd14_ablate.pdf",
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined_cd14_ablate, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  "outs/control/figures/05_pbmc_balanced_combined_batch_cd14_ablate.pdf",
  width = 6,
  height = 6
)

# Downsample CD14+ Monocytes from batch 2 and replot figures
cd_14_batch_2_indices <- (
  which(pbmc_combined@meta.data$celltype %in% ("CD14+ Monocyte") &
          pbmc_combined@meta.data$batch %in% ("Batch 2"))
)
cd_14_batch_2_indices_sample_90 <- sample(
  x = cd_14_batch_2_indices,
  size = round(length(cd_14_batch_2_indices)*0.9, 0)
)
pbmc_combined_cd14_ds <- pbmc_combined[, -(cd_14_batch_2_indices_sample_90)]

DimPlot(pbmc_combined_cd14_ds, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  "outs/control/figures/05_pbmc_balanced_combined_celltypes_cd14_10_pct_ds.pdf",
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined_cd14_ds, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  "outs/control/figures/05_pbmc_balanced_combined_batch_cd14_10_pct_ds.pdf",
  width = 6,
  height = 6
)

