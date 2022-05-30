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

# Change to results dir for PDAC comp data 
setwd("../../../results/pdac_comp/")

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)
gc()

# Load in and concatenate celltype summary files
setwd("../celltype_imbalance_summaries")
cimba_files <- list.files()
cimba_loaded <- lapply(cimba_files, fread)
cimba_concat <- Reduce(rbind, cimba_loaded)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)
gc()

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
clus_concord_loaded <- lapply(clus_concord_files, fread)
clus_concord_concat <- Reduce(rbind, clus_concord_loaded)
gc()

# Load in and concatenate dge concordance summaries
setwd("../dge_concord_stats/")
dge_files <- list.files()
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)
gc()

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("../dge_ranking_stats_marker_sub")
dge_rank_files <- list.files()
dge_rank_loaded <- lapply(dge_rank_files, fread)
dge_rank_concat <- Reduce(rbind, dge_rank_loaded)
gc()

# Load in markers and corresponding celltypes 
setwd("../marker_results/")
base_marker_genes <- fread(
  "peng_pdac_tumor_annot_8_batch_preintegration_marker_selection.tsv"
) 

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
gc()

# Change to top level dir 
setwd("../../..")

# Make pdac comp output dir if it doesn't exist
if (!dir.exists("outs/pdac_comp/figures")) {
  dir.create("outs/pdac_comp/figures", recursive = TRUE)
}
if (!dir.exists("outs/pdac_comp/results")) {
  dir.create("outs/pdac_comp/results", recursive = TRUE)
}

### Fig 6A - Differences in number of clusters before and after downsampling

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

# Create plot of effects on cluster number, facetted by method, filled by
# type (control, ablation, downsampling), and plotted based on downsampled
# celltype/compartment
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
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~Method, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Compartment downsampled",
    y = "Number of Leiden clusters post-integration"
  ) +
  coord_flip() +
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
  "outs/pdac_comp/figures/12_pdac_compartment_ds_ablate_methods_nclusters.pdf",
  width = 12,
  height = 7,
  device = cairo_pdf
)

### Fig 6B) - summary of ablation and downsampling effects on batch and 
### celltype/compartment ARI values (base metrics), dependant on method 

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

# Get median celltype/compartment ARI based on each method and whether or not
# it's a control, downsampling, or ablation, and by celltype/compartment
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
  "Epithelial normal" = dark_2_cols[1], 
  "Epithelial tumor" = dark_2_cols[2], 
  "Microenvironment" = dark_2_cols[3],
  "None" = "Black"
)

ht1 = Heatmap(
  as.matrix(median_celltype_ari_long_vals_only_scaled), 
  name = "Scaled median \ncompartment ARI", 
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
  name = "Affected compartment",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
celltype_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/pdac_comp/figures/12_celltype_ari_ds_effects_heatmap.pdf",
  width = 8, 
  height = 4
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
# it's a control, downsampling, or ablation, and by celltype/compartment
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
  "Epithelial normal" = dark_2_cols[1], 
  "Epithelial tumor" = dark_2_cols[2], 
  "Microenvironment" = dark_2_cols[3],
  "None" = "Black"
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
  name = "Affected compartment",
  col = col_celltype,
  width = unit(0.5, "cm"),  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE
)
batch_ari_hm <- ht1 + ht2 + ht3
CairoPDF(
  "outs/pdac_comp/figures/12_batch_ari_ds_effects_heatmap.pdf", 
  width = 8, 
  height = 4
)
draw(
  batch_ari_hm,
  column_title = "Integration method",
  column_title_side = "bottom",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)
dev.off()

### Fig 6C - summary of downsampling results on KNN classification scores
### of given subsets 

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

# Subset for only cases where the celltype/compartment downsampled is equal to 
# the celltype being classified
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

# Create function to format facet labels (downsampled celltypes/compartments)
ds_celltype_labelled <- function(variable,value){
  return(paste0("Compartment affected = ", value))
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
    ncol = 3
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected compartment F1-classification score post-integration"
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
  "outs/pdac_comp/figures/12_pdac_compartment_ds_ablate_methods_knn_f1_score.pdf",
  width = 16,
  height = 7,
  device = cairo_pdf
)

### Supplementary figure - plot the total number of cells across all batches
### and celltypes for the control data

# Indicate which cases are control and which are ds or ablate 
cimba_concat$type <- ifelse(
  cimba_concat$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    cimba_concat$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Pick just one control case - they all have the same values
cimba_concat_selected <- cimba_concat[
  (cimba_concat$Replicate == 0) &
    (cimba_concat$type == "Control") &
      (cimba_concat$`Proportion downsampled` == 0)
]

# Format the dataframe
cimba_concat_selected_formatted <- reshape2::melt(
  cimba_concat_selected,
  id.vars = "celltype",
  measure.vars = c(
    "celltype_count_batch_0",
    "celltype_count_batch_1",
    "celltype_count_batch_2",
    "celltype_count_batch_3",
    "celltype_count_batch_4",
    "celltype_count_batch_5",
    "celltype_count_batch_6",
    "celltype_count_batch_7"
  )
)
colnames(cimba_concat_selected_formatted) <- c("celltype", "batch", "n_cells")
cimba_concat_selected_formatted$batch <- plyr::mapvalues(
  cimba_concat_selected_formatted$batch,
  from = c(
    "celltype_count_batch_0",
    "celltype_count_batch_1",
    "celltype_count_batch_2",
    "celltype_count_batch_3",
    "celltype_count_batch_4",
    "celltype_count_batch_5",
    "celltype_count_batch_6",
    "celltype_count_batch_7"
  ),
  to = c(
    "Batch 1",
    "Batch 2",
    "Batch 3",
    "Batch 4",
    "Batch 5",
    "Batch 6",
    "Batch 7",
    "Batch 8"
  )
)

# Plot the results across batches 
ggplot(
  data = cimba_concat_selected_formatted, 
  aes(
    x = batch,
    y = n_cells
  )) +
  geom_bar(stat = "identity", aes(fill = celltype), position = "dodge2") +
  scale_fill_manual( 
    breaks = c("Epithelial normal", "Epithelial tumor", "Microenvironment"),
    values = c(dark_2_cols[1], dark_2_cols[2], dark_2_cols[3])
  ) +
  labs(
    fill = "Compartment",
    x = "",
    y = "Number of cells"
  ) +
  coord_flip() +
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
  "outs/pdac_comp/figures/12_pdac_compartment_batches_celltype_n_cells.pdf",
  width = 7,
  height = 5,
  device = cairo_pdf
)
