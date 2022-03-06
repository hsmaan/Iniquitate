library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(dotwhisker)
library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)

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

### Fig 1A) - plotting of pbmc balanced dataset and examples of downsampling
### and ablation on CD14 Monocyte cells 

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

### Fig 2A) - summary of ablation and downsampling effects on batch and 
### celltype ARI values (base metrics), dependant on method/technique 

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

# Get median celltype ARI based on each method and whether or not
# it's a control, downsampling, or ablation, and by celltype 
median_celltype_ari_results <- imba_clus_merged %>% 
  group_by(Method, type, `Downsampled celltypes`) %>% 
  summarize(
    `Median celltype ARI` = median(`Celltype ARI Imbalanced`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Melt and format for ComplexHeatMap plotting 
median_celltype_ari_results_vals_long <- reshape2::dcast(
  median_celltype_ari_results,
  formula = type + `Downsampled celltypes` ~ `Method`,
  value.var = "Median celltype ARI"
)
median_celltype_ari_results_vals_long$type <- factor(
  median_celltype_ari_results_vals_long$type,
  levels = c("Control", "Downsampled", "Ablated")
)
median_celltype_ari_results_vals_long <- median_celltype_ari_results_vals_long[
  order(
    median_celltype_ari_results_vals_long$type,
    median_celltype_ari_results_vals_long$`Downsampled celltypes`
  ),
]
rownames(median_celltype_ari_results_vals_long) <- c(
  paste0("R_", seq(1, nrow(median_celltype_ari_results_vals_long)))
)
colnames(median_celltype_ari_results_vals_long)[1] <- c(
  "Type"
)
median_celltype_ari_long_vals_only <- median_celltype_ari_results_vals_long[
  ,-c(1,2)
]
median_celltype_ari_long_vals_only_scaled <- scale(
  median_celltype_ari_long_vals_only,
  center = TRUE,
  scale = TRUE
)
median_celltype_ari_long_type <- median_celltype_ari_results_vals_long[
  ,1, drop = FALSE
]
median_celltype_ari_long_celltype <- median_celltype_ari_results_vals_long[
  ,2, drop = FALSE
]

# Plot the three heatmaps together for median celltype ARI post integration
dark_2_cols = palette.colors(n = 8, "Dark2")
col_type = c(
  "Control" = "forestgreen",
  "Downsampled" = "darkorchid3",
  "Ablated" = "firebrick2"
)
col_celltype = c(
  "B cell" = dark_2_cols[1], 
  "CD14+ Monocyte" = dark_2_cols[2], 
  "CD4+ T cell" = dark_2_cols[3],
  "CD8+ T cell" = dark_2_cols[4],
  "FCGR3A+ Monocyte" = dark_2_cols[5],
  "NK cell" = dark_2_cols[6],
  "None" = "black"
)

ht1 = Heatmap(
  as.matrix(median_celltype_ari_long_vals_only_scaled), 
  name = "Scaled median \ncelltype ARI", 
  width = unit(5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)
ht2 = Heatmap(
  as.matrix(median_celltype_ari_long_type), 
  name = "Type",
  col = col_type,
  width = unit(0.5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
ht3 = Heatmap(
  as.matrix(median_celltype_ari_long_celltype), 
  name = "Affected celltype",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
celltype_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/05_celltype_ari_ds_effects_heatmap.pdf", 
  width = 8, 
  height = 6
)
draw(
  celltype_ari_hm,
  column_title = "Integration method",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

# Perform the exact same analysis/heatmap as above, but now for Batch ARI

# Get median batch ARI based on each method and whether or not
# it's a control, downsampling, or ablation, and by celltype 
median_batch_ari_results <- imba_clus_merged %>% 
  group_by(Method, type, `Downsampled celltypes`) %>% 
  summarize(
    `Median batch ARI` = median(`Batch ARI`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Melt and format for ComplexHeatMap plotting 
median_batch_ari_results_vals_long <- reshape2::dcast(
  median_batch_ari_results,
  formula = type + `Downsampled celltypes` ~ `Method`,
  value.var = "Median batch ARI"
)
median_batch_ari_results_vals_long$type <- factor(
  median_batch_ari_results_vals_long$type,
  levels = c("Control", "Downsampled", "Ablated")
)
median_batch_ari_results_vals_long <- median_batch_ari_results_vals_long[
  order(
    median_batch_ari_results_vals_long$type,
    median_batch_ari_results_vals_long$`Downsampled celltypes`
  ),
]
rownames(median_batch_ari_results_vals_long) <- c(
  paste0("R_", seq(1, nrow(median_batch_ari_results_vals_long)))
)
colnames(median_batch_ari_results_vals_long)[1] <- c(
  "Type"
)
median_batch_ari_long_vals_only <- median_batch_ari_results_vals_long[
  ,-c(1,2)
]
median_batch_ari_long_vals_only_scaled <- scale(
  median_batch_ari_long_vals_only,
  center = TRUE,
  scale = TRUE
)
median_batch_ari_long_type <- median_batch_ari_results_vals_long[
  ,1, drop = FALSE
]
median_batch_ari_long_celltype <- median_batch_ari_results_vals_long[
  ,2, drop = FALSE
]

# Plot the three heatmaps together for median batch ARI post integration
dark_2_cols = palette.colors(n = 8, "Dark2")
col_type = c(
  "Control" = "forestgreen",
  "Downsampled" = "darkorchid3",
  "Ablated" = "firebrick2"
)
col_celltype = c(
  "B cell" = dark_2_cols[1], 
  "CD14+ Monocyte" = dark_2_cols[2], 
  "CD4+ T cell" = dark_2_cols[3],
  "CD8+ T cell" = dark_2_cols[4],
  "FCGR3A+ Monocyte" = dark_2_cols[5],
  "NK cell" = dark_2_cols[6],
  "None" = "black"
)

ht1 = Heatmap(
  as.matrix(median_batch_ari_long_vals_only_scaled), 
  name = "Scaled median \nbatch ARI", 
  width = unit(5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)
ht2 = Heatmap(
  as.matrix(median_batch_ari_long_type), 
  name = "Type",
  col = col_type,
  width = unit(0.5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
ht3 = Heatmap(
  as.matrix(median_batch_ari_long_celltype), 
  name = "Affected celltype",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
batch_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/05_batch_ari_ds_effects_heatmap.pdf", 
  width = 8, 
  height = 6
)
draw(
  batch_ari_hm,
  column_title = "Integration method",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

