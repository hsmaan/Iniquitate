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
library(networkD3)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Kevin's palette for plotting many catagoricals 
kev_palette <- c(
  "dodgerblue2", "#E31A1C",
  "green4",
  "#6A3D9A", 
  "#FF7F00", 
  "black", "gold1",
  "skyblue2", "#FB9A99", 
  "palegreen2",
  "#CAB2D6", 
  "#FDBF6F", 
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# Change to results dir for control data w/ cider implemented 
setwd("../../../results/control_w_cider/")

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
clus_concat <- clus_concat[clus_concat$Method != "liger"]
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
clus_concord_concat <- clus_concord_concat[
  clus_concord_concat$`Method 1` != "liger"
]
clus_concord_concat <- clus_concord_concat[
  clus_concord_concat$`Method 2` != "liger"
]
gc()

# Load in and concatenate cluster classification results 
setwd("../clustering_classification/")
clus_class_files <- list.files()
clus_class_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_class_files,
  value = TRUE
)
clus_class_loaded <- lapply(clus_class_files, fread)
clus_class_concat <- Reduce(rbind, clus_class_loaded)
clus_class_concat <- clus_class_concat[
  clus_class_concat$Method != "liger"
]
gc()

# Load in and concatenate KNN-based cluster classification results
setwd("../knn_clustering_classification/")
knn_clus_class_files <- list.files()
knn_clus_class_files <- grep(
  "pbmc_2_batch_base_balanced",
  knn_clus_class_files,
  value = TRUE
)
knn_clus_class_loaded <- lapply(knn_clus_class_files, fread)
knn_clus_class_concat <- Reduce(rbind, knn_clus_class_loaded)
knn_clus_class_concat <- knn_clus_class_concat[
  knn_clus_class_concat$Method != "liger"
]

# Change to top level dir 
setwd("../../../")

# Create the cider output directories (figures/results) if it doesn't exist
if (!dir.exists("outs/control_w_cider/figures")) {
  dir.create("outs/control_w_cider/figures", recursive = TRUE)
}
if (!dir.exists("outs/control_w_cider/results")) {
  dir.create("outs/control_w_cider/results", recursive = TRUE)
}

# Begin by plotting the ARI scores for each method 

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
  "outs/control_w_cider/figures/17_celltype_ari_ds_effects_heatmap_no_liger.pdf",
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
  "outs/control_w_cider/figures/17_celltype_balanced_ari_ds_effects_heatmap_no_liger.pdf",
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
  "outs/control_w_cider/figures/17_batch_ari_ds_effects_heatmap_no_liger.pdf", 
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

# Do the same analysis above for imbalanced cell-type ARI, but this time don't
# scale the ARI values - keep the absolute differences between methods 
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

imba_clus_merged$type <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

median_celltype_ari_results <- imba_clus_merged %>% 
  group_by(Method, type, `Downsampled celltypes`) %>% 
  summarize(
    `Median celltype ARI` = median(`Celltype ARI Imbalanced`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

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

median_celltype_ari_long_type <- median_celltype_ari_results_vals_long[
  ,1, drop = FALSE
]
median_celltype_ari_long_celltype <- median_celltype_ari_results_vals_long[
  ,2, drop = FALSE
]

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
  as.matrix(median_celltype_ari_long_vals_only), 
  name = "Median \ncell-type ARI", 
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
  "outs/control_w_cider/figures/17_celltype_ari_unscaled_ds_effects_heatmap_no_liger.pdf",
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

# Repeat the same analysis above, but now with unscaled median batch ARI
median_batch_ari_results <- imba_clus_merged %>% 
  group_by(Method, type, `Downsampled celltypes`) %>% 
  summarize(
    `Median batch ARI` = median(`Batch ARI`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

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
median_batch_ari_long_type <- median_batch_ari_results_vals_long[
  ,1, drop = FALSE
]
median_batch_ari_long_celltype <- median_batch_ari_results_vals_long[
  ,2, drop = FALSE
]

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
  as.matrix(median_batch_ari_long_vals_only), 
  name = "Median \nbatch ARI", 
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
  "outs/control_w_cider/figures/17_batch_ari_unscaled_ds_effects_heatmap_no_liger.pdf", 
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

# Repeat the same analysis above, but now with unscaled median balanced celltype ARI 
median_balanced_ari_results <- imba_clus_merged %>% 
  group_by(Method, type, `Downsampled celltypes`) %>% 
  summarize(
    `Median cell-type balanced ARI` = median(`Celltype ARI Balanced`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

median_balanced_ari_results_vals_long <- reshape2::dcast(
  median_balanced_ari_results,
  formula = type + `Downsampled celltypes` ~ `Method`,
  value.var = "Median cell-type balanced ARI"
)
median_balanced_ari_results_vals_long$type <- factor(
  median_balanced_ari_results_vals_long$type,
  levels = c("Control", "Downsampled", "Ablated")
)
median_balanced_ari_results_vals_long <- median_balanced_ari_results_vals_long[
  order(
    median_balanced_ari_results_vals_long$type,
    median_balanced_ari_results_vals_long$`Downsampled celltypes`
  ),
]
rownames(median_balanced_ari_results_vals_long) <- c(
  paste0("R_", seq(1, nrow(median_balanced_ari_results_vals_long)))
)
colnames(median_balanced_ari_results_vals_long)[1] <- c(
  "Type"
)
median_balanced_ari_long_vals_only <- median_balanced_ari_results_vals_long[
  ,-c(1,2)
]
median_balanced_ari_long_type <- median_balanced_ari_results_vals_long[
  ,1, drop = FALSE
]
median_balanced_ari_long_celltype <- median_balanced_ari_results_vals_long[
  ,2, drop = FALSE
]

ht1 = Heatmap(
  as.matrix(median_balanced_ari_long_vals_only), 
  name = "Median cell-type unscaled \nbalanced ARI", 
  width = unit(5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE
)
ht2 = Heatmap(
  as.matrix(median_balanced_ari_long_type), 
  name = "Type",
  col = col_type,
  width = unit(0.5, "cm"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
ht3 = Heatmap(
  as.matrix(median_balanced_ari_long_celltype), 
  name = "Affected cell-type",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)

celltype_balanced_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control_w_cider/figures/17_celltype_balanced_ari_unscaled_ds_effects_heatmap_no_liger.pdf", 
  width = 8, 
  height = 6
)
draw(
  celltype_balanced_ari_hm,
  column_title = "Integration method",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

# Results of celltype downsampling and ablation on  
# cluster-based classification accuracy 

# Merge imbalance and cluster classification results together
imba_clus_merged <- merge(
  imba_concat,
  clus_class_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_clus_merged <- distinct(imba_clus_merged)

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_clus_merged_celltype <- imba_clus_merged[
  imba_clus_merged$Celltype == imba_clus_merged$`Downsampled celltypes` |
    imba_clus_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_clus_merged_celltype$type <- ifelse(
  imba_clus_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format celltype names
imba_clus_merged_celltype$Celltype <- plyr::mapvalues(
  imba_clus_merged_celltype$Celltype,
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

ggplot(data = imba_clus_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
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
  "outs/control_w_cider/figures/17_pbmc_ds_ablate_allmethod_clus_f1_score_no_liger.pdf",
  width = 16,
  height = 14,
  device = cairo_pdf
)

ggplot(data = imba_clus_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
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
  "outs/control_w_cider/figures/17_pbmc_ds_ablate_allmethod_clus_f1_score_no_liger_0_1_y.pdf",
  width = 16,
  height = 14,
  device = cairo_pdf
)

# Get and save the standard deviation of the median values of F1, based on
# method utilized, celltype downsampled, and subset utilized
imba_clus_merged_celltype_medians  <- imba_clus_merged_celltype %>% 
  group_by(Method, Celltype, type) %>% 
  summarize(
    `Median F1-score per subset` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

imba_clus_merged_celltype_medians_stdev  <- imba_clus_merged_celltype_medians %>% 
  group_by(Celltype) %>% 
  summarize(
    `Stdev Method Median F1-score per subset` = sd(
      `Median F1-score per subset`, 
      na.rm = FALSE
    ),
    .groups = "keep"
  ) %>%
  as.data.frame()

colnames(imba_clus_merged_celltype_medians_stdev) <- c(
  "Cell-type", 
  "Standard deviation of medians for F1-score, across methods, replicates, and experiment types"
)
fwrite(
  imba_clus_merged_celltype_medians_stdev,
  "outs/control_w_cider/results/17_baseline_pbmc_clus_class_f1_score_stdevs_per_celltype.tsv",
  sep = "\t", 
  row.names = FALSE,
  col.names = TRUE
)

# Results of celltype downsampling and ablation on  
# KNN-cluster-based classification accuracy
# The difference between this analysis is that the previous used logistic
# regression on the cluster labels while this one uses KNN on the cluster labels 

# Merge imbalance and knn cluster classification results together
imba_knn_clus_merged <- merge(
  imba_concat,
  knn_clus_class_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_knn_clus_merged <- distinct(imba_knn_clus_merged)

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_clus_merged_celltype <- imba_knn_clus_merged[
  imba_knn_clus_merged$Celltype == imba_knn_clus_merged$`Downsampled celltypes` |
    imba_knn_clus_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_clus_merged_celltype$type <- ifelse(
  imba_knn_clus_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_clus_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format celltype names
imba_knn_clus_merged_celltype$Celltype <- plyr::mapvalues(
  imba_knn_clus_merged_celltype$Celltype,
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

ggplot(data = imba_knn_clus_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
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
  "outs/control_w_cider/figures/17_pbmc_ds_ablate_allmethod_knn_clus_f1_score_no_liger.pdf",
  width = 16,
  height = 14,
  device = cairo_pdf
)

ggplot(data = imba_knn_clus_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
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
  "outs/control_w_cider/figures/17_pbmc_ds_ablate_allmethod_knn_clus_f1_score_no_liger_0_1_y.pdf",
  width = 16,
  height = 14,
  device = cairo_pdf
)

# Get and save the standard deviation of the median values of F1, based on
# method utilized, celltype downsampled, and subset utilized
imba_knn_clus_merged_celltype_medians  <- imba_knn_clus_merged_celltype %>% 
  group_by(Method, Celltype, type) %>% 
  summarize(
    `Median F1-score per subset` = median(`F1-score`, na.rm = FALSE),
    .groups = "keep"
  ) %>%
  as.data.frame()

imba_knn_clus_merged_celltype_medians_stdev  <- imba_knn_clus_merged_celltype_medians %>% 
  group_by(Celltype) %>% 
  summarize(
    `Stdev Method Median F1-score per subset` = sd(
      `Median F1-score per subset`, 
      na.rm = FALSE
    ),
    .groups = "keep"
  ) %>%
  as.data.frame()

colnames(imba_knn_clus_merged_celltype_medians_stdev) <- c(
  "Cell-type", 
  "Standard deviation of medians for F1-score, across methods, replicates, and experiment types"
)
fwrite(
  imba_knn_clus_merged_celltype_medians_stdev,
  "outs/control_w_cider/results/17_baseline_pbmc_knn_clus_f1_score_stdevs_per_celltype.tsv",
  sep = "\t", 
  row.names = FALSE,
  col.names = TRUE
)