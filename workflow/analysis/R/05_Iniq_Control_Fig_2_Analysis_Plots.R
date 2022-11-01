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

# Change to PBMC data dir  
setwd("resources/h5ad_files/int_datasets/pbmc_2_batch_base_balanced/")

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

### Fig 2A) - plotting of pbmc balanced dataset and examples of downsampling
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

### Fig 2B) - summary of ablation and downsampling effects on batch and 
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

### Fig 2C) - results of celltype downsampling and ablation on  
### KNN classification scores 

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

# Format celltype names
imba_knn_merged_celltype$Celltype <- plyr::mapvalues(
  imba_knn_merged_celltype$Celltype,
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

# Create function to format facet labels (downsampled celltypes)
ds_celltype_labelled <- function(variable,value){
  return(paste0("Cell-type affected = ", value))
}

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  ylim(0, 1) +
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 2
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected celltype F1-classification score post-integration"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/05_pbmc_ds_ablate_allmethod_knn_f1_score.pdf",
  width = 12,
  height = 14,
  device = cairo_pdf
)

### Fig 2D) - correlation of batch and celltype ARI values with KNN 
### classification scores to show discordance of these results 

# Merge imbalance, clustering, and knn classification results together
imba_knn_cluster_merged <- merge(
  imba_knn_merged,
  clus_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Method"
  )
)
imba_knn_cluster_merged <- distinct(imba_knn_cluster_merged)

# Format celltype names
imba_knn_cluster_merged$Celltype <- plyr::mapvalues(
  imba_knn_cluster_merged$Celltype,
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

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_cluster_merged$type <- ifelse(
  imba_knn_cluster_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_cluster_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Get the median celltype classification score for each of these,
# grouped by method, type and celltype ARI score
median_imba_knn_cluster_results_cari <- imba_knn_cluster_merged %>% 
  group_by(Method, type, `Celltype ARI Imbalanced`) %>% 
  summarize(
    `Median F1-score` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Plot the result for correlation with the Celltype ARI values 
# THIS IS THE MEDIAN SCORE - must be indicated in figure legend
ggscatter(median_imba_knn_cluster_results_cari, 
          x = "Celltype ARI Imbalanced", 
          y = "Median F1-score", 
          size = 0.4,
          combine = TRUE,
          xlab = "Celltype ARI post-integration",
          ylab = "Celltype F1-classification score post-integration",
          palette = "jco",
          add = "reg.line", 
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE) +
  facet_wrap(.~`Method`, scales = "fixed") +
  stat_cor(
    aes(label = ..r.label..),
    method = "spearman", 
    label.x = 0.6, 
    label.y = 0.4, 
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
  "outs/control/figures/05_pbmc_ds_ablate_f1_score_celltype_ari_corr.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Get the median celltype classification score for each of these,
# grouped by method, type and Batch ARI score
median_imba_knn_cluster_results_bari <- imba_knn_cluster_merged %>% 
  group_by(Method, type, `Batch ARI`) %>% 
  summarize(
    `Median F1-score` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Plot the result for correlation with the Batch ARI values 
# THIS IS THE MEDIAN SCORE - must be indicated in figure legend
ggscatter(median_imba_knn_cluster_results_bari, 
          x = "Batch ARI", 
          y = "Median F1-score", 
          size = 0.4,
          combine = TRUE,
          xlab = "Batch ARI post-integration",
          ylab = "Celltype F1-classification score post-integration",
          palette = "jco",
          add = "reg.line", 
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE) +
  facet_wrap(.~`Method`, scales = "fixed") +
  stat_cor(
    aes(label = ..r.label..),
    method = "spearman", 
    label.x = 0.97, 
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
  "outs/control/figures/05_pbmc_ds_ablate_f1_score_batch_ari_corr.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

##### Analysis of PBMC 2 batch balanced data - hierarchical agglomeration ##### 

# Remove all previous objects and garbage collect
rm(list = ls())
gc()

# Change to results dir for control data 
setwd("results/control/")

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  clus_files,
  value = TRUE
)
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
clus_concord_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  clus_concord_files,
  value = TRUE
)
clus_concord_loaded <- lapply(clus_concord_files, fread)
clus_concord_concat <- Reduce(rbind, clus_concord_loaded)

# Load in and concatenate dge concordance summaries
setwd("../dge_concord_stats/")
dge_files <- list.files()
dge_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  dge_files,
  value = TRUE
)
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)

# Load in and concatenate knn classification summaries
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

# Change to PBMC data dir  
setwd("resources/h5ad_files/int_datasets/pbmc_2_batch_hierarchical_balanced/")

# Convert balanced pbmc data files to h5seurat format
Convert(
  "tran_exp5_pbmc_batch1_balanced_hierarchical.h5ad", 
  dest = "h5seurat", 
  overwrite = TRUE
)
Convert(
  "tran_exp5_pbmc_batch2_balanced_hierarchical.h5ad",
  dest = "h5seurat",
  overwrite = TRUE
)

# Load h5ad files for both balanced pbmc datasets 
pbmc_1 <- SeuratDisk::LoadH5Seurat(
  "tran_exp5_pbmc_batch1_balanced_hierarchical.h5seurat"
)
pbmc_2 <- SeuratDisk::LoadH5Seurat(
  "tran_exp5_pbmc_batch2_balanced_hierarchical.h5seurat"
)

# Change to top level dir 
setwd("../../../../")

### Fig 2E) - plotting of pbmc balanced hierarchical dataset and examples of 
### downsampling and ablation on NK/T Cells 

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
  "outs/control/figures/05_pbmc_balanced_hierarchical_combined_celltypes.pdf",
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
  "outs/control/figures/05_pbmc_balanced_hierarchical_combined_batch.pdf",
  width = 6,
  height = 6
)

# Remove (ablate) NK/T-Cells from batch 2 and replot figures 
pbmc_combined_nkt_ablate <- pbmc_combined[, 
  -(which(pbmc_combined@meta.data$celltype %in% ("NK/T cell") &
            pbmc_combined@meta.data$batch %in% ("Batch 2")))
]
DimPlot(pbmc_combined_nkt_ablate, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_hierarchical_combined_celltypes_nkt_ablate.pdf"
  ),
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined_nkt_ablate, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_hierarchical_combined_batch_nkt_ablate.pdf"
  ),
  width = 6,
  height = 6
)

# Downsample NK/T Cells from batch 2 and replot figures
nkt_batch_2_indices <- (
  which(pbmc_combined@meta.data$celltype %in% ("NK/T cell") &
          pbmc_combined@meta.data$batch %in% ("Batch 2"))
)
nkt_batch_2_indices_sample_90 <- sample(
  x = nkt_batch_2_indices,
  size = round(length(nkt_batch_2_indices)*0.9, 0)
)
pbmc_combined_nkt_ds <- pbmc_combined[, -(nkt_batch_2_indices_sample_90)]

DimPlot(pbmc_combined_nkt_ds, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_hierarchical_combined_celltypes_nkt_10_pct_ds.pdf"
  ),
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined_nkt_ds, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_combined_batch_nkt_10_pct_ds.pdf"
  ),
  width = 6,
  height = 6
)

### Fig 2E) - plotting of KNN classification results with the hierarchical
### celltype results 

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

# Create function to format facet labels (downsampled celltypes)
ds_celltype_labelled <- function(variable,value){
  return(paste0("Cell-type affected = ", value))
}

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 1
  ) +
  ylim(0, 1) + 
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected celltype F1-classification score post-integration"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/05_pbmc_hierarchical_ds_ablate_allmethod_knn_f1_score.pdf",
  width = 8,
  height = 14,
  device = cairo_pdf
)

#############################################################################
### REDO ALL PLOTS ABOVE BUT WITHOUT LIGER ###
#############################################################################

gc()
rm(list =ls())
`%ni%` <- Negate(`%in%`)
setwd("results/control/")

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
dge_concat <- dge_concat[
  dge_concat$`Method 1` != "liger"
]
dge_concat <- dge_concat[
  dge_concat$`Method 2` != "liger"
]

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

# Change to PBMC data dir  
setwd("resources/h5ad_files/int_datasets/pbmc_2_batch_base_balanced/")

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

### Fig 2A) - plotting of pbmc balanced dataset and examples of downsampling
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

### Fig 2B) - summary of ablation and downsampling effects on batch and 
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
  name = "Scaled median \ncell-type ARI", 
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
  name = "Affected cell-type",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
celltype_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/05_celltype_ari_ds_effects_heatmap_no_liger.pdf",
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

# Perform the exact same analysis/heatmap as above, but for the balanced 
# celltype ARI 

median_balanced_celltype_ari_results <- imba_clus_merged %>% 
  group_by(Method, type, `Downsampled celltypes`) %>% 
  summarize(
    `Median balanced celltype ARI` = median(`Celltype ARI Balanced`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Melt and format for ComplexHeatMap plotting 
median_balanced_celltype_ari_results_vals_long <- reshape2::dcast(
  median_balanced_celltype_ari_results,
  formula = type + `Downsampled celltypes` ~ `Method`,
  value.var = "Median balanced celltype ARI"
)
median_balanced_celltype_ari_results_vals_long$type <- factor(
  median_balanced_celltype_ari_results_vals_long$type,
  levels = c("Control", "Downsampled", "Ablated")
)
median_balanced_celltype_ari_results_vals_long <- 
  median_balanced_celltype_ari_results_vals_long[
    order(
      median_balanced_celltype_ari_results_vals_long$type,
      median_balanced_celltype_ari_results_vals_long$`Downsampled celltypes`
    ),
  ]
rownames(median_balanced_celltype_ari_results_vals_long) <- c(
    paste0("R_", seq(1, nrow(median_balanced_celltype_ari_results_vals_long)))
  )
colnames(median_balanced_celltype_ari_results_vals_long)[1] <- c(
    "Type"
)
median_balanced_celltype_ari_long_vals_only <- 
  median_balanced_celltype_ari_results_vals_long[
    ,-c(1,2)
  ]
median_balanced_celltype_ari_long_vals_only_scaled <- scale(
  median_balanced_celltype_ari_long_vals_only,
  center = TRUE,
  scale = TRUE
)
median_balanced_celltype_ari_long_type <- 
  median_balanced_celltype_ari_results_vals_long[
    ,1, drop = FALSE
  ]
median_balanced_celltype_ari_long_celltype <- 
  median_balanced_celltype_ari_results_vals_long[
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
  as.matrix(median_balanced_celltype_ari_long_vals_only_scaled), 
  name = "Scaled median \ncell-type balanced ARI", 
  width = unit(5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)
ht2 = Heatmap(
  as.matrix(median_balanced_celltype_ari_long_type), 
  name = "Type",
  col = col_type,
  width = unit(0.5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
ht3 = Heatmap(
  as.matrix(median_balanced_celltype_ari_long_celltype), 
  name = "Affected cell-type",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
celltype_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/05_celltype_balanced_ari_ds_effects_heatmap_no_liger.pdf",
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
  name = "Affected cell-type",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
batch_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/05_batch_ari_ds_effects_heatmap_no_liger.pdf", 
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

### Fig 2C) - results of celltype downsampling and ablation on  
### KNN classification scores 

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

# Format celltype names
imba_knn_merged_celltype$Celltype <- plyr::mapvalues(
  imba_knn_merged_celltype$Celltype,
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

# Create function to format facet labels (downsampled celltypes)
ds_celltype_labelled <- function(variable,value){
  return(paste0("Cell-type affected = ", value))
}

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 2
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected cell-type F1-classification score post-integration"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/05_pbmc_ds_ablate_allmethod_knn_f1_score_no_liger.pdf",
  width = 12,
  height = 14,
  device = cairo_pdf
)

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  ylim(0, 1) +
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 2
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected cell-type F1-classification score post-integration"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/05_pbmc_ds_ablate_allmethod_knn_f1_score_no_liger_0_1_y.pdf",
  width = 12,
  height = 14,
  device = cairo_pdf
)

# Get and save the standard deviation of the median values of F1, based on
# method utilized, celltype downsampled, and subset utilized
imba_knn_merged_celltype_medians  <- imba_knn_merged_celltype %>% 
  group_by(Method, Celltype, type) %>% 
  summarize(
    `Median F1-score per subset` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

imba_knn_merged_celltype_medians_stdev  <- imba_knn_merged_celltype_medians %>% 
  group_by(Celltype) %>% 
  summarize(
    `Stdev Method Median F1-score per subset` = sd(
      `Median F1-score per subset`, 
      na.rm = FALSE
    ),
    .groups = "keep"
  ) %>%
  as.data.frame()

colnames(imba_knn_merged_celltype_medians_stdev) <- c(
  "Cell-type", 
  "Standard deviation of medians for F1-score, across methods, replicates, and experiment types"
)
fwrite(
  imba_knn_merged_celltype_medians_stdev,
  "outs/control/results/05_baseline_pbmc_f1_score_stdevs_per_celltype.tsv",
  sep = "\t", 
  row.names = FALSE,
  col.names = TRUE
)

### Fig 2D) - correlation of batch and celltype ARI values with KNN 
### classification scores to show discordance of these results 

# Merge imbalance, clustering, and knn classification results together
imba_knn_cluster_merged <- merge(
  imba_knn_merged,
  clus_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Method"
  )
)
imba_knn_cluster_merged <- distinct(imba_knn_cluster_merged)

# Format celltype names
imba_knn_cluster_merged$Celltype <- plyr::mapvalues(
  imba_knn_cluster_merged$Celltype,
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

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_cluster_merged$type <- ifelse(
  imba_knn_cluster_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_cluster_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Get the median celltype classification score for each of these,
# grouped by method, type and celltype ARI score
median_imba_knn_cluster_results_cari <- imba_knn_cluster_merged %>% 
  group_by(Method, type, `Celltype ARI Imbalanced`) %>% 
  summarize(
    `Median F1-score` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Plot the result for correlation with the Celltype ARI values 
# THIS IS THE MEDIAN SCORE - must be indicated in figure legend
ggscatter(median_imba_knn_cluster_results_cari, 
          x = "Celltype ARI Imbalanced", 
          y = "Median F1-score", 
          size = 0.4,
          combine = TRUE,
          xlab = "Cell-type ARI post-integration",
          ylab = "Cell-type F1-classification score post-integration",
          palette = "jco",
          add = "reg.line", 
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE) +
  facet_wrap(.~`Method`, scales = "fixed") +
  stat_cor(
    aes(label = ..r.label..),
    method = "spearman", 
    label.x = 0.6, 
    label.y = 0.4, 
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
  "outs/control/figures/05_pbmc_ds_ablate_f1_score_celltype_ari_corr_no_liger.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Get the median celltype classification score for each of these,
# grouped by method, type and Batch ARI score
median_imba_knn_cluster_results_bari <- imba_knn_cluster_merged %>% 
  group_by(Method, type, `Batch ARI`) %>% 
  summarize(
    `Median F1-score` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

# Plot the result for correlation with the Batch ARI values 
# THIS IS THE MEDIAN SCORE - must be indicated in figure legend
ggscatter(median_imba_knn_cluster_results_bari, 
          x = "Batch ARI", 
          y = "Median F1-score", 
          size = 0.4,
          combine = TRUE,
          xlab = "Batch ARI post-integration",
          ylab = "Cell-type F1-classification score post-integration",
          palette = "jco",
          add = "reg.line", 
          add.params = list(color = "blue", fill = "lightgray"),
          conf.int = TRUE) +
  facet_wrap(.~`Method`, scales = "fixed") +
  stat_cor(
    aes(label = ..r.label..),
    method = "spearman", 
    label.x = 0.97, 
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
  "outs/control/figures/05_pbmc_ds_ablate_f1_score_batch_ari_corr_no_liger.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

##### Analysis of PBMC 2 batch balanced data - hierarchical agglomeration ##### 

# Remove all previous objects and garbage collect
rm(list = ls())
gc()

# Change to results dir for control data 
setwd("results/control/")

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_files <- grep(
  "pbmc_2_batch_hierarchical_balanced",
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
  "pbmc_2_batch_hierarchical_balanced",
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
  "pbmc_2_batch_hierarchical_balanced",
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
  "pbmc_2_batch_hierarchical_balanced",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
knn_concat <- knn_concat[knn_concat$Method != "liger"]

# Change to top level dir 
setwd("../../..")

# Change to PBMC data dir  
setwd("resources/h5ad_files/int_datasets/pbmc_2_batch_hierarchical_balanced/")

# Convert balanced pbmc data files to h5seurat format
Convert(
  "tran_exp5_pbmc_batch1_balanced_hierarchical.h5ad", 
  dest = "h5seurat", 
  overwrite = TRUE
)
Convert(
  "tran_exp5_pbmc_batch2_balanced_hierarchical.h5ad",
  dest = "h5seurat",
  overwrite = TRUE
)

# Load h5ad files for both balanced pbmc datasets 
pbmc_1 <- SeuratDisk::LoadH5Seurat(
  "tran_exp5_pbmc_batch1_balanced_hierarchical.h5seurat"
)
pbmc_2 <- SeuratDisk::LoadH5Seurat(
  "tran_exp5_pbmc_batch2_balanced_hierarchical.h5seurat"
)

# Change to top level dir 
setwd("../../../../")

### Fig 2E) - plotting of pbmc balanced hierarchical dataset and examples of 
### downsampling and ablation on NK/T Cells 

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
  "outs/control/figures/05_pbmc_balanced_hierarchical_combined_celltypes.pdf",
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
  "outs/control/figures/05_pbmc_balanced_hierarchical_combined_batch.pdf",
  width = 6,
  height = 6
)

# Remove (ablate) NK/T-Cells from batch 2 and replot figures 
pbmc_combined_nkt_ablate <- pbmc_combined[, 
                                          -(which(pbmc_combined@meta.data$celltype %in% ("NK/T cell") &
                                                    pbmc_combined@meta.data$batch %in% ("Batch 2")))
]
DimPlot(pbmc_combined_nkt_ablate, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_hierarchical_combined_celltypes_nkt_ablate.pdf"
  ),
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined_nkt_ablate, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_hierarchical_combined_batch_nkt_ablate.pdf"
  ),
  width = 6,
  height = 6
)

# Downsample NK/T Cells from batch 2 and replot figures
nkt_batch_2_indices <- (
  which(pbmc_combined@meta.data$celltype %in% ("NK/T cell") &
          pbmc_combined@meta.data$batch %in% ("Batch 2"))
)
nkt_batch_2_indices_sample_90 <- sample(
  x = nkt_batch_2_indices,
  size = round(length(nkt_batch_2_indices)*0.9, 0)
)
pbmc_combined_nkt_ds <- pbmc_combined[, -(nkt_batch_2_indices_sample_90)]

DimPlot(pbmc_combined_nkt_ds, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_hierarchical_combined_celltypes_nkt_10_pct_ds.pdf"
  ),
  width = 6, 
  height = 6
)
DimPlot(pbmc_combined_nkt_ds, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/control/figures/",
    "05_pbmc_balanced_combined_batch_nkt_10_pct_ds.pdf"
  ),
  width = 6,
  height = 6
)

### Fig 2E) - plotting of KNN classification results with the hierarchical
### celltype results 

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

# Create function to format facet labels (downsampled celltypes)
ds_celltype_labelled <- function(variable,value){
  return(paste0("Cell-type affected = ", value))
}

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  ylim(0.5, 1) + 
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 1
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected cell-type F1-classification score post-integration"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/05_pbmc_hierarchical_ds_ablate_allmethod_knn_f1_score_no_liger_ylim_0.5_1.pdf",
  width = 8,
  height = 14,
  device = cairo_pdf
)

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  ylim(0, 1) +
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 1
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected cell-type F1-classification score post-integration"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/05_pbmc_hierarchical_ds_ablate_allmethod_knn_f1_score_no_liger_y_0_1.pdf",
  width = 8,
  height = 14,
  device = cairo_pdf
)

