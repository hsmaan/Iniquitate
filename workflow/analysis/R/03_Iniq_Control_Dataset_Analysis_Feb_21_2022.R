library(ggplot2)
library(ggthemes)
library(Seurat)
library(SeuratDisk)

# Change to PBMC data dir  
setwd("../../../resources/h5ad_files/int_datasets/pbmc_2_batch_base_balanced")

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

# Create plots of celltype embeddings for the two datasets independantly
DimPlot(pbmc_1, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank())
ggsave(
  "../../../../outs/control/figures/03_pbmc_balanced_batch_1_celltypes.pdf",
  width = 7,
  height = 6
)

DimPlot(pbmc_2, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank())
ggsave(
  "../../../../outs/control/figures/03_pbmc_balanced_batch_2_celltypes.pdf",
  width = 7,
  height = 6
)

# Create plots of celltype embeddings for the two datasets combined
DimPlot(pbmc_combined, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank())
ggsave(
  "../../../../outs/control/figures/03_pbmc_balanced_combined_celltypes.pdf",
  width = 7,
  height = 6
)
DimPlot(pbmc_combined, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank())
ggsave(
  "../../../../outs/control/figures/03_pbmc_balanced_combined_batch.pdf",
  width = 7,
  height = 6
)
