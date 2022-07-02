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

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("../dge_ranking_stats_marker_sub")
dge_rank_files <- list.files()
dge_rank_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_rank_files,
  value = TRUE
)
dge_rank_loaded <- lapply(dge_rank_files, fread)
dge_rank_concat <- Reduce(rbind, dge_rank_loaded)
gc()

# Load in markers and corresponding celltypes 
setwd("../marker_results/")
base_marker_genes <- fread(
  "pbmc_2_batch_base_balanced_preintegration_marker_selection.tsv"
) 

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

# Load in and concatenate raw annotation results 
setwd("../annotation_results/")
anno_files <- list.files()
anno_files <- grep(
  "pbmc_2_batch_base_balanced",
  anno_files,
  value = TRUE
)
anno_loaded <- lapply(anno_files, fread)
anno_concat <- Reduce(rbind, anno_loaded)
gc()

# Load in and concatenate annotation summary results 
setwd("../annotation_scores")
anno_scores_files <- list.files()
anno_scores_files <- grep(
  "pbmc_2_batch_base_balanced",
  anno_scores_files,
  value = TRUE
)
anno_scores_loaded <- lapply(anno_scores_files, fread)
anno_scores_concat <- Reduce(rbind, anno_scores_loaded)
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
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
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
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/07_pbmc_ds_ablate_allmethod_cluster_number.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)  

### Supplmentary figure - Correlation between the number of clusters per method 
### and the adjusted rand index for celltype and batch ARI, on a per-method 
### level
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

### Figures 3B-3D - concordance of DGE for marker genes before and after 
### downsampling 

# Merge together imbalanced and dge ranking datasets 
imba_dge_merged <- merge(
  imba_concat,
  dge_rank_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
) 
imba_dge_merged <- distinct(imba_dge_merged)

# Indicate which samples are controls and which are real runs
imba_dge_merged$type <- ifelse(
  imba_dge_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_dge_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

### Supplementary figures 

# Plot the max rank of the marker gene across all methods (method agnostic) -
# ordered by median variability of the rank - Supplementary figure
gene_rank_variance <- imba_dge_merged %>% 
  group_by(Gene) %>%
  summarize(var = sd(`Max rank`))
gene_rank_variance_sorted <- gene_rank_variance[
  order(gene_rank_variance$var, decreasing = TRUE),
]
ggplot(data = imba_dge_merged, aes(
  x = `Gene`,
  y = `Max rank`
)
) +
  geom_boxplot(
    aes(fill = type),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  facet_grid(
    .~factor(type, levels = c("Control", "Downsampled", "Ablated")), 
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Maximum rank in differential expression analysis across clusters"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(gene_rank_variance_sorted$Gene)
  ) + 
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "None")
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "marker_gene_max_ranks.pdf"
  ),
  width = 14,
  height = 12,
  device = cairo_pdf
)

# Plot the max rank of the marker gene across all methods, taking method into
# account - ordered by median variability of the rank. Omit gene names  
ggplot(data = imba_dge_merged, aes(
  x = `Gene`,
  y = `Max rank`
)
) +
  geom_boxplot(
    aes(fill = type),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  facet_grid(
    Method ~ .,
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  scale_color_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Maximum rank in differential expression analysis across clusters"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(gene_rank_variance_sorted$Gene)
  ) + 
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) + 
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "None")
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "marker_gene_max_ranks_per_method.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

# For each method, subset, plot and save the variability results
# Supplementary figure
methods <- sort(unique(imba_dge_merged$Method))
dge_method_plt <- function(imba_dge_df, method) {
  imba_dge_df_method <- imba_dge_df[which(imba_dge_df$Method %in% method)]
  ggplot(
    data = imba_dge_df_method, 
    aes(
      x = `Gene`,
      y = `Max rank`
    )
  ) +
    geom_boxplot(
      fill = "dodgerblue2",
      notch = FALSE,
      alpha = 0.8 
    ) + 
    facet_wrap(
      .~factor(type, levels = c("Control", "Downsampled", "Ablated")), 
      scales = "fixed"
    ) +
    scale_fill_manual( 
      breaks = c("Control", "Downsampled", "Ablated"),
      values = c("forestgreen", "darkorchid3", "firebrick2")
    ) +
    labs(
      title = method,
      fill = "Type",
      x = "Marker gene",
      y = "Maximum rank in differential expression analysis across clusters"
    ) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(factor(imba_dge_df$Gene)))) + 
    theme_few() +
    theme(plot.title = element_text(size = 24)) +
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16)) +
    theme(strip.text.x = element_text(size = 16)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(legend.title = element_text(size = 16)) +
    theme(legend.text = element_text(size = 14))
}
lapply(methods, function(x) {
  dge_method_plt(imba_dge_merged, x)
  ggsave(
    paste0(
      "outs/control/figures/07_pbmc_ds_ablate_",
      x,
      "_dge_rankings.pdf"
    ),
    width = 14,
    height = 10,
    device = cairo_pdf
  )
})

### Fig 3B - Heatmap of marker gene perturbation scores across methods,
### for control, downsampled and ablated cases 

# Get the standard deviation of each marker gene, with method and status 
# factored in 
gene_rank_var_method_type <- imba_dge_merged %>% 
  group_by(Gene, Method, type) %>%
  summarize(var = sd(`Max rank`))

# Create three separate matrices for each of the three types (control, 
# downsampled, ablated)
gene_rank_var_method_type_ctrl <- gene_rank_var_method_type[
  which(gene_rank_var_method_type$type %in% "Control"), 
]
gene_rank_var_method_type_ctrl_long <- reshape2::dcast(
  gene_rank_var_method_type_ctrl,
  Gene~Method,
  value.var = "var"
)
gene_rank_var_method_type_ctrl_long_mat <- as.matrix(
  gene_rank_var_method_type_ctrl_long[, -1]
)
rownames(gene_rank_var_method_type_ctrl_long_mat) <- 
  gene_rank_var_method_type_ctrl_long[, 1]

gene_rank_var_method_type_ds <- gene_rank_var_method_type[
  which(gene_rank_var_method_type$type %in% "Downsampled"), 
]
gene_rank_var_method_type_ds_long <- reshape2::dcast(
  gene_rank_var_method_type_ds,
  Gene~Method,
  value.var = "var"
)
gene_rank_var_method_type_ds_long_mat <- as.matrix(
  gene_rank_var_method_type_ds_long[, -1]
)
rownames(gene_rank_var_method_type_ds_long_mat) <- 
  gene_rank_var_method_type_ds_long[, 1]

gene_rank_var_method_type_abla <- gene_rank_var_method_type[
  which(gene_rank_var_method_type$type %in% "Ablated"), 
]
gene_rank_var_method_type_abla_long <- reshape2::dcast(
  gene_rank_var_method_type_abla,
  Gene~Method,
  value.var = "var"
)
gene_rank_var_method_type_abla_long_mat <- as.matrix(
  gene_rank_var_method_type_abla_long[, -1]
)
rownames(gene_rank_var_method_type_abla_long_mat) <- 
  gene_rank_var_method_type_abla_long[, 1]

# Create three separate heatmaps for each of the types
# (use max 10 here as this indicates significant enough changes)
col_pert = circlize::colorRamp2(
  c(
    0,
    10
  ),
  c("white", "firebrick1")
)
ht1 <- Heatmap(
  gene_rank_var_method_type_ctrl_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Control",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  col = col_pert
)
ht2 <- Heatmap(
  gene_rank_var_method_type_ds_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Downsampled",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)
ht3 <- Heatmap(
  gene_rank_var_method_type_abla_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Ablated",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)

# Create a heatmap indicating which celltype is associated with each
# marker gene
base_marker_gene_dup_added <- base_marker_genes %>%
  group_by(`Top 10 marker genes (union across batches)`) %>%
  summarize(Celltype = paste0(unique(Celltype), collapse = ', ')) %>%
  as.data.frame
base_marker_gene_dup_added_mat <- as.matrix(
  base_marker_gene_dup_added[, -1]
)
rownames(base_marker_gene_dup_added_mat) <- base_marker_gene_dup_added[, 1] 

# Color palette 
unique_dup_celltypes_sorted <- sort(unique(base_marker_gene_dup_added$Celltype))
unique_kev_palette_vals <- kev_palette[1:length(unique_dup_celltypes_sorted)]
col_celltype = unique_kev_palette_vals
names(col_celltype) <- unique_dup_celltypes_sorted

# Celltype heatmap 
ht4 <- Heatmap(
  base_marker_gene_dup_added_mat, 
  name = "Celltype associated \nwith marker gene", 
  width = unit(1, "cm"),
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = TRUE,
  col = col_celltype
)

# Plot and save all of the heatmaps together 
marker_pert_hm <- ht1 + ht2 + ht3 + ht4
CairoPDF(
  "outs/control/figures/07_marker_gene_pert_pbmc_control_heatmap.pdf",
  width = 14, 
  height = 6
)
draw(
  marker_pert_hm,
  column_title = "Integration method",
  column_title_side = "bottom",
  row_title_side = "right",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold")
)
dev.off()

# Plot one version without the celltype labels 
ht1 <- Heatmap(
  gene_rank_var_method_type_ctrl_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Control",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  col = col_pert
)
ht2 <- Heatmap(
  gene_rank_var_method_type_ds_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Downsampled",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)
ht3 <- Heatmap(
  gene_rank_var_method_type_abla_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Ablated",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)
marker_per_hm_no_celltype <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/07_marker_gene_pert_pbmc_control_heatmap_no_ctype.pdf",
  width = 14, 
  height = 6
)
draw(
  marker_per_hm_no_celltype,
  column_title = "Integration method",
  column_title_side = "bottom",
  row_title_side = "right",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
)
dev.off()

### Fig 3C - Concordance of DGE genes - top 10 most variable marker
### genes across methods 

# Pick 10 exemplary genes that show the highest variability - take top 10 with 
# the exclusion of the mitochondrial gene 
top_10_variable_dge <- gene_rank_variance_sorted$Gene[1:11]
top_10_variable_dge <- top_10_variable_dge[-6]

# Subset data for top 10 with the exclusion of mitochondrial gene and 
# plot results for ranking change with ablation and downsampling 
imba_dge_merged_top_10_var_genes <- imba_dge_merged[
  imba_dge_merged$Gene %in% top_10_variable_dge
]
ggplot(imba_dge_merged_top_10_var_genes, aes(x = Gene, y = `Max rank`)) +
  geom_boxplot(
    aes(
      fill = factor(type, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(.~Method, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Maximum rank in differential expression analysis across clusters"
  ) +
  coord_flip() +
  theme_few() +
  scale_x_discrete(limits = rev(top_10_variable_dge)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_dge_rankings_top_10_var_genes.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Get stdev in rank of each gene, subset by type and method 
gene_rank_variance_grouped <- imba_dge_merged %>% 
  group_by(Gene, Method, type) %>%
  summarize(`Max rank stdev` = sd(`Max rank`)) %>%
  as.data.table

# Get the plots for each method, and the standard deviations of max
# rank of the genes
methods <- sort(unique(imba_dge_merged$Method))
dge_stdev_method_plt <- function(imba_dge_stdev_df, method) {
  imba_dge_stdev_df_sub <- imba_dge_stdev_df[
    which(imba_dge_stdev_df$Method %in% method)
  ]
  ggplot(
    data = imba_dge_stdev_df_sub, 
    aes(
      x = `Gene`,
      y = `Max rank stdev`
    ) 
  ) +
    geom_bar(
      aes(
        fill = factor(
        imba_dge_stdev_df_sub$type, 
        levels = c("Control", "Downsampled", "Ablated")
        )
      ),
      stat = "identity",
      alpha = 0.8 
    ) + 
    facet_wrap(
      .~factor(
          imba_dge_stdev_df_sub$type, 
          levels = c("Control", "Downsampled", "Ablated")
      ), 
      scales = "fixed"
    ) +
    scale_fill_manual( 
      breaks = c("Control", "Downsampled", "Ablated"),
      values = c("forestgreen", "darkorchid3", "firebrick2")
    ) +
    labs(
      fill = "Type",
      x = "Marker gene",
      y = "Standard deviation of maximum rank in differential expression"
    ) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(factor(imba_dge_stdev_df_sub$Gene)))) + 
    theme_few() +
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16)) +
    theme(strip.text.x = element_text(size = 16)) +
    theme(plot.title = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(legend.title = element_text(size = 16)) +
    theme(legend.text = element_text(size = 14))
}
lapply(methods, function(x) {
  dge_stdev_method_plt(gene_rank_variance_grouped, x)
  ggsave(
    paste0(
      "outs/control/figures/07_pbmc_ds_ablate_",
      x,
      "_dge_rankings_method_type_stdev.pdf"
    ),
    width = 14,
    height = 10,
    device = cairo_pdf
  )
})

# Plot the top 10 most variable genes in terms of rank (across all methods)
# based on type and method 
gene_rank_variance_grouped_top_10_sub <- gene_rank_variance_grouped[
  which(gene_rank_variance_grouped$Gene %in% top_10_variable_dge)
]
ggplot(
  data = gene_rank_variance_grouped_top_10_sub, 
  aes(
    x = `Gene`,
    y = `Max rank stdev`
  ) 
) +
  geom_bar(
    aes(
      fill = factor(
        gene_rank_variance_grouped_top_10_sub$type, 
        levels = c("Control", "Downsampled", "Ablated")
      )
    ),
    stat = "identity",
    alpha = 0.8,
    position = "dodge2"
  ) + 
  facet_wrap(
    .~Method,
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Standard deviation of maximum rank in differential expression"
  ) +
  coord_flip() +
  theme_few() +
  scale_x_discrete(limits = rev(top_10_variable_dge)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "_dge_rankings_method_type_stdev_top_10.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Get stdev in rank of each gene, subset by type only (no method) 
gene_rank_variance_grouped_nomethod <- imba_dge_merged %>% 
  group_by(Gene, type) %>%
  summarize(`Max rank stdev` = sd(`Max rank`)) %>%
  as.data.table

# Plot the values of standard deviations of each gene across all methods 
ggplot(
  data = gene_rank_variance_grouped_nomethod, 
  aes(
    x = `Gene`,
    y = `Max rank stdev`
  ) 
) +
  geom_bar(
    aes(
      fill = factor(
        gene_rank_variance_grouped_nomethod$type, 
        levels = c("Control", "Downsampled", "Ablated")
      )
    ),
    stat = "identity",
    alpha = 0.8 
  ) + 
  facet_wrap(
    .~factor(
      gene_rank_variance_grouped_nomethod$type, 
      levels = c("Control", "Downsampled", "Ablated")
    ), 
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Standard deviation of maximum rank in differential expression"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(levels(factor(gene_rank_variance_grouped_nomethod$Gene)))
  ) + 
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "_dge_rankings_type_stdev.pdf"
  ),
  width = 14,
  height = 10,
  device = cairo_pdf
)

# Plot the correspondence of the top 10 most highly variable marker genes with
# their respective celltypes
top_10_variable_dge_markersub <- base_marker_genes[
  which(
    base_marker_genes$`Top 10 marker genes (union across batches)` %in% 
      top_10_variable_dge
    )
]

# Order the results 
top_10_variable_dge_markersub <- top_10_variable_dge_markersub[
  match( 
    top_10_variable_dge,
    top_10_variable_dge_markersub$`Top 10 marker genes (union across batches)`
  )
]
colnames(top_10_variable_dge_markersub) <- c(
  "Associated celltype",
  "Gene"
)

# Add the standard deviation values across subsets overall 
top_10_variable_dge_markersub <- merge(
  top_10_variable_dge_markersub,
  gene_rank_variance_sorted,
  join = "left",
  by = c(
    "Gene"
  ),
  all.y = FALSE
)

# Order based on the standard deviation
top_10_variable_dge_markersub <- top_10_variable_dge_markersub[
  order(top_10_variable_dge_markersub$var, decreasing = )
]

# Plot the SankeyNetwork diagram, save and then convert saved html to pdf
nodes <- data.frame(
  name = c(
    as.character(top_10_variable_dge_markersub$Gene), 
     as.character(top_10_variable_dge_markersub$`Associated celltype`) %>% 
      unique()
  )
)
top_10_variable_dge_markersub$IDsource <- match(
  top_10_variable_dge_markersub$Gene, nodes$name)-1 
top_10_variable_dge_markersub$IDtarget <- match(
  top_10_variable_dge_markersub$`Associated celltype`, nodes$name)-1

top_10_var_marker_genes_sankey <- sankeyNetwork(
  Links = top_10_variable_dge_markersub,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "var",
  NodeID = "name",
  sinksRight = FALSE
)
saveNetwork(
  top_10_var_marker_genes_sankey,
  file = paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "_top_var_dges_with_assoc_celltypes.html"
  ),
  selfcontained = TRUE
)

### Fig 3D) - Correlation of marker gene instability and the celltype
### downsampled - determine if this can impact the instability or 
### if there is a correlation

# Get stdev in rank of each gene, subset by type, method, downsampled celltype
# and marker(s) associated with celltype
gene_rank_variance_grouped_celltype_specific <- imba_dge_merged %>% 
  group_by(Gene, Method, type, `Downsampled celltypes`) %>%
  summarize(`Max rank stdev` = sd(`Max rank`)) %>%
  as.data.table

# Merge the DGE stdev summary stats with the marker data indicating
# which celltype each marker corresponds to - ALLOW A CARTESIAN PRODUCT
# HERE AS EACH MARKER MAY BE SPECIFIC TO MORE THAN ONE CELLTYPE. The cartesian
# product should be valid, as we are considering correlations between celltype-
# specific markers and change in DGE status
base_marker_genes_copy <- base_marker_genes
colnames(base_marker_genes_copy) <- c(
  "Associated celltype",
  "Gene"
)
gene_rank_variance_grouped_celltype_specific_marker <- merge(
  gene_rank_variance_grouped_celltype_specific,
  base_marker_genes_copy,
  by = c(
    "Gene"
  ),
  allow.cartesian = TRUE
)
gene_rank_variance_grouped_celltype_specific_marker <- distinct(
  gene_rank_variance_grouped_celltype_specific_marker
)

# Format celltype names for plotting 
gene_rank_variance_grouped_celltype_specific_marker$`Downsampled celltypes` <-
  plyr::mapvalues(
    gene_rank_variance_grouped_celltype_specific_marker$`Downsampled celltypes`,
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

gene_rank_variance_grouped_celltype_specific_marker$`Associated celltype` <-
  plyr::mapvalues(
    gene_rank_variance_grouped_celltype_specific_marker$`Associated celltype`,
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

# For first plot, collapse/summarize further by averaging across all of the 
# genes (mean - not median)
gene_rank_variance_grouped_celltype_specific_marker_mean <-
  gene_rank_variance_grouped_celltype_specific_marker %>%
  group_by(Method, type, `Downsampled celltypes`, `Associated celltype`) %>%
  summarize(`Mean max rank stdev` = mean(`Max rank stdev`)) %>%
  as.data.table

# Remove 'None' from this - only considering cases where the downsampled 
# celltype is equivalent to the associated celltype of the given marker 
gene_rank_variance_grouped_celltype_specific_marker_mean <-
  gene_rank_variance_grouped_celltype_specific_marker_mean[
    gene_rank_variance_grouped_celltype_specific_marker_mean$
      `Downsampled celltypes` %ni% "None"
  ]

# Plot heatmaps specific to the overall results of each method - first for 
# the downsampled results 
# MARKER GENE PERTURBATION SCORE HAS TO BE DEFINED IN RESULTS/METHODS
gene_rank_variance_grouped_celltype_specific_marker_mean_ds <- 
  gene_rank_variance_grouped_celltype_specific_marker_mean[
    gene_rank_variance_grouped_celltype_specific_marker_mean$type %in% 
      "Downsampled"
  ]

# Clip values to 10, as this is considered a fair "maximum acceptable" 
# perturbation
gene_rank_variance_grouped_celltype_specific_marker_mean_ds$
  `Mean max rank stdev clipped` <- pmin(
    10,
    gene_rank_variance_grouped_celltype_specific_marker_mean_ds$`Mean max rank stdev`
  )
  
ggplot(
  data = gene_rank_variance_grouped_celltype_specific_marker_mean_ds,
  aes(
    x = `Downsampled celltypes`,
    y = `Associated celltype`
  ) 
) + 
  geom_tile(
    aes(fill = `Mean max rank stdev clipped`)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkorchid3"
  ) +
  facet_wrap(.~Method, scales = "free") +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(
    size = 12, 
    angle = 90, 
    vjust = 1, 
    hjust = 1)
  ) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  labs(
    fill = "Average marker gene \nperturbation score",
    x = "Downsampled celltype",
    y = "Celltype associated with marker genes"
  )
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_only_",
    "_dge_rankings_celltype_marker_celltype_ds_compare.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

# Plot heatmaps specific to the overall results of each method - now for 
# the ablated results 
# MARKER GENE PERTURBATION SCORE HAS TO BE DEFINED IN RESULTS/METHODS
gene_rank_variance_grouped_celltype_specific_marker_mean_ablated <- 
  gene_rank_variance_grouped_celltype_specific_marker_mean[
    gene_rank_variance_grouped_celltype_specific_marker_mean$type %in% 
      "Ablated"
  ]

# Clip values to 10, as this is considered a fair "maximum acceptable" 
# perturbation
gene_rank_variance_grouped_celltype_specific_marker_mean_ablated$
  `Mean max rank stdev clipped` <- pmin(
    10,
    gene_rank_variance_grouped_celltype_specific_marker_mean_ablated$`Mean max rank stdev`
  )

ggplot(
  data = gene_rank_variance_grouped_celltype_specific_marker_mean_ablated,
  aes(
    x = `Downsampled celltypes`,
    y = `Associated celltype`
  ) 
) + 
  geom_tile(
    aes(fill = `Mean max rank stdev clipped`)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "firebrick2"
  ) +
  facet_wrap(.~Method, scales = "free") +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(
    size = 12, 
    angle = 90, 
    vjust = 1, 
    hjust = 1)
  ) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  labs(
    fill = "Average marker gene \nperturbation score",
    x = "Ablated celltype",
    y = "Celltype associated with marker genes"
  )
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ablated_only_",
    "_dge_rankings_celltype_marker_celltype_ablated_compare.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

### Fig 3E) - Impact of imbalance on annotation of data using an external
### reference

# Combine annotation score and imbalance dataframes
imba_anno_merged <- merge(
  imba_concat,
  anno_scores_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
) 
imba_anno_merged <- distinct(imba_anno_merged)

# Indicate which samples are controls and which are real runs
imba_anno_merged$type <- ifelse(
  imba_anno_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_anno_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format celltype names for plotting 
imba_anno_merged$`Downsampled celltypes` <- plyr::mapvalues(
  imba_anno_merged$`Downsampled celltypes`,
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

# Reformat to get the celltype-specific F1-scores for both L1 and L2 
# annotations 
imba_anno_merged_score_format <- imba_anno_merged[
  ,
  c(
    "Downsampled celltypes",
    "B cell", 
    "CD4 T cell", 
    "CD8 T cell", 
    "Monocyte_CD14", 
    "Monocyte_FCGR3A",
    "NK cell",
    "Overall balanced accuracy",
    "Overall F1-score",
    "Subset",
    "Score type",
    "type"
  ),
  with = FALSE
] 
imba_anno_merged_score_format <- melt(
  imba_anno_merged_score_format,
  id.vars = c(
    "Downsampled celltypes",
    "Overall balanced accuracy",
    "Overall F1-score",
    "Subset",
    "Score type",
    "type"
  ),
  value.name = "Score",
  variable.name = "Celltype scored"
) 
imba_anno_merged_score_format <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$`Score type` %in% "f1-score", 
]

# Format celltype names for plotting 
imba_anno_merged_score_format$`Celltype scored` <- plyr::mapvalues(
  imba_anno_merged_score_format$`Celltype scored`,
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
  
# Plot the celltype specific downsampling results, for both L1, and L2 
# annotations - F1 scores as a function of type, and celltype downsampled, 
# for each annotation
imba_anno_merged_score_format_l1 <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L1", 
]
imba_anno_merged_score_format_l2 <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L2", 
]

ggplot(data = imba_anno_merged_score_format_l1, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L1 annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_f1_scores_celltype_ds.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

ggplot(data = imba_anno_merged_score_format_l2, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L2 annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_f1_scores_celltype_ds.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

## Focus on T-cells here, plot the variability of T cell annotation
## accuracy within downsampled, ablated, and control results

# Get the variance of the F1-score, based on celltype scored and 
# downsampling status
imba_anno_merged_score_format_l1_sd <- imba_anno_merged_score_format_l1 %>% 
  group_by(`Celltype scored`, type) %>%
  summarize(F1_sd = sd(`Score`)) 
imba_anno_merged_score_format_l2_sd <- imba_anno_merged_score_format_l2 %>% 
  group_by(`Celltype scored`, type) %>%
  summarize(F1_sd = sd(`Score`)) 

# Plot these results, emphasizing high variance of F1-scores for T cells 
ggplot(data = imba_anno_merged_score_format_l1_sd, aes(
  x = `Celltype scored`,
  y = F1_sd
)) +
  geom_bar(
    aes(fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated"))),
    alpha = 0.8,
    stat = "identity",
    position = "dodge",
    width = 3
  ) + 
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype",
    y = "Standard deviation of L1 F1-scores"
  ) +
  theme_few() +
  coord_flip() +
  scale_x_discrete(
    limits = rev(imba_anno_merged_score_format_l1_sd$`Celltype scored`)
  ) + 
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_f1_scores_celltype_sdev.pdf"
  ),
  width = 7,
  height = 5,
  device = cairo_pdf
)  

ggplot(data = imba_anno_merged_score_format_l2_sd, aes(
  x = `Celltype scored`,
  y = F1_sd
)) +
  geom_bar(
    aes(fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated"))),
    alpha = 0.8,
    stat = "identity",
    position = "dodge",
    width = 3
  ) + 
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype",
    y = "Standard deviation of L2 F1-scores"
  ) +
  theme_few() +
  coord_flip() +
  scale_x_discrete(
    limits = rev(imba_anno_merged_score_format_l2_sd$`Celltype scored`)
  ) + 
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_f1_scores_celltype_sdev.pdf"
  ),
  width = 7,
  height = 5,
  device = cairo_pdf
)  

# Determine the celltypes that are being predicted for both CD4 and CD8 
# T cells, L1 and L2, for control, downsampled and ablated experiments 

# Merge all of the annotation results and imbalance results
imba_anno_merged_all <- merge(
  imba_concat,
  anno_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
)

# Indicate which samples are controls and which are real runs
imba_anno_merged_all$type <- ifelse(
  imba_anno_merged_all$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_anno_merged_all$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format celltype names for plotting
imba_anno_merged_all$`Downsampled celltypes` <- plyr::mapvalues(
  imba_anno_merged_all$`Downsampled celltypes`,
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
imba_anno_merged_all$`Real celltype` <- plyr::mapvalues(
  imba_anno_merged_all$`Real celltype`,
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

# Subset for L1 and L2 cases for only the T cells 
imba_anno_merged_all_tcell <- imba_anno_merged_all[
  imba_anno_merged_all$`Real celltype` %in% c("CD4+ T cell", "CD8+ T cell") 
]

# Plot the results for the celltypes that are most predicted for the T-cell
# subsets - for L1 and L2 
ggplot(data = imba_anno_merged_all_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L1`
)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Real celltype`) +
  labs(
    fill = "Predicted L1 \ncelltype",
    x = "Type",
    y = "Proportion of L1 predictions"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_t_cell_preds_barplot.pdf"
  ),
  width = 10,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L2`
)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Real celltype`) +
  labs(
    fill = "Predicted L2 \ncelltype",
    x = "Type",
    y = "Proportion of L2 predictions"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_t_cell_preds_barplot.pdf"
  ),
  width = 10,
  height = 7,
  device = cairo_pdf
) 

# Plot the same results as above, separately this time for CD4 and CD8 T cells,
# and L1 and L2 predictions, but conditioned on the celltype that has been
# downsampled

# Subsample data for both celltypes
imba_anno_merged_all_cd4_tcell <- imba_anno_merged_all_tcell[
  imba_anno_merged_all_tcell$`Real celltype` %in% "CD4+ T cell"
]
imba_anno_merged_all_cd8_tcell <- imba_anno_merged_all_tcell[
  imba_anno_merged_all_tcell$`Real celltype` %in% "CD8+ T cell"
]

ggplot(data = imba_anno_merged_all_cd4_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L1`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L1 \ncelltype",
    x = "",
    y = "Proportion of L1 predictions for CD4+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_cd4t_cell_preds",
    "_celltype_ds_specific_barplots.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_cd4_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L2`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L2 \ncelltype",
    x = "",
    y = "Proportion of L2 predictions for CD4+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_cd4t_cell_preds",
    "_celltype_ds_specific_barplots.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_cd8_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L1`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L1 \ncelltype",
    x = "",
    y = "Proportion of L1 predictions for CD8+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_cd8t_cell_preds",
    "_celltype_ds_specific_barplots.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_cd8_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L2`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L2 \ncelltype",
    x = "",
    y = "Proportion of L2 predictions for CD8+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_cd8t_cell_preds",
    "_celltype_ds_specific_barplots.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 


# Plot the celltype specific downsampling results, for the L1 and L2 baseline
# annotations (using just HVG KNN) - F1 scores as a function of type, and 
# celltype downsampled, for each annotation
imba_anno_merged_score_format_l1_baseline <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L1 baseline", 
]
imba_anno_merged_score_format_l2_baseline <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L2 baseline", 
]

ggplot(data = imba_anno_merged_score_format_l1_baseline, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_point(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L1 baseline annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_baseline_annotation_f1_scores_celltype_ds.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

ggplot(data = imba_anno_merged_score_format_l2_baseline, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_point(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L2 baseline annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_baseline_annotation_f1_scores_celltype_ds.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

#############################################################################
### REDO ALL FIGURES ABOVE WITHOUT LIGER ###
#############################################################################

rm(list = ls())
gc()

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
dge_concat <- dge_concat[dge_concat$`Method 1` != "liger"]
dge_concat <- dge_concat[dge_concat$`Method 2` != "liger"]
gc()

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("../dge_ranking_stats_marker_sub")
dge_rank_files <- list.files()
dge_rank_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_rank_files,
  value = TRUE
)
dge_rank_loaded <- lapply(dge_rank_files, fread)
dge_rank_concat <- Reduce(rbind, dge_rank_loaded)
dge_rank_concat <- dge_rank_concat[dge_rank_concat$Method != "liger"]
gc()

# Load in markers and corresponding celltypes 
setwd("../marker_results/")
base_marker_genes <- fread(
  "pbmc_2_batch_base_balanced_preintegration_marker_selection.tsv"
) 

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
gc()

# Load in and concatenate raw annotation results 
setwd("../annotation_results/")
anno_files <- list.files()
anno_files <- grep(
  "pbmc_2_batch_base_balanced",
  anno_files,
  value = TRUE
)
anno_loaded <- lapply(anno_files, fread)
anno_concat <- Reduce(rbind, anno_loaded)
gc()

# Load in and concatenate annotation summary results 
setwd("../annotation_scores")
anno_scores_files <- list.files()
anno_scores_files <- grep(
  "pbmc_2_batch_base_balanced",
  anno_scores_files,
  value = TRUE
)
anno_scores_loaded <- lapply(anno_scores_files, fread)
anno_scores_concat <- Reduce(rbind, anno_scores_loaded)
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
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
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
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16))
ggsave(
  "outs/control/figures/07_pbmc_ds_ablate_allmethod_cluster_number_no_liger.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)  

### Supplmentary figure - Correlation between the number of clusters per method 
### and the adjusted rand index for celltype and batch ARI, on a per-method 
### level
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
  "outs/control/figures/07_pbmc_ds_ablate_clus_num_celltype_ari_corr_no_liger.pdf",
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
  "outs/control/figures/07_pbmc_ds_ablate_clus_num_batch_ari_corr_no_liger.pdf",
  width = 14,
  height = 8,
  device = cairo_pdf
)

### Figures 3B-3D - concordance of DGE for marker genes before and after 
### downsampling 

# Merge together imbalanced and dge ranking datasets 
imba_dge_merged <- merge(
  imba_concat,
  dge_rank_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
) 
imba_dge_merged <- distinct(imba_dge_merged)

# Indicate which samples are controls and which are real runs
imba_dge_merged$type <- ifelse(
  imba_dge_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_dge_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

### Supplementary figures 

# Plot the max rank of the marker gene across all methods (method agnostic) -
# ordered by median variability of the rank - Supplementary figure
gene_rank_variance <- imba_dge_merged %>% 
  group_by(Gene) %>%
  summarize(var = sd(`Max rank`))
gene_rank_variance_sorted <- gene_rank_variance[
  order(gene_rank_variance$var, decreasing = TRUE),
]
ggplot(data = imba_dge_merged, aes(
  x = `Gene`,
  y = `Max rank`
)
) +
  geom_boxplot(
    aes(fill = type),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  facet_grid(
    .~factor(type, levels = c("Control", "Downsampled", "Ablated")), 
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Maximum rank in differential expression analysis across clusters"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(gene_rank_variance_sorted$Gene)
  ) + 
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "None")
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "marker_gene_max_ranks_no_liger.pdf"
  ),
  width = 14,
  height = 12,
  device = cairo_pdf
)

# Plot the max rank of the marker gene across all methods, taking method into
# account - ordered by median variability of the rank. Omit gene names  
ggplot(data = imba_dge_merged, aes(
  x = `Gene`,
  y = `Max rank`
)
) +
  geom_boxplot(
    aes(fill = type),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  facet_grid(
    Method ~ .,
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  scale_color_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Maximum rank in differential expression analysis across clusters"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(gene_rank_variance_sorted$Gene)
  ) + 
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) + 
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  theme(legend.position = "None")
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "marker_gene_max_ranks_per_method_no_liger.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

# For each method, subset, plot and save the variability results
# Supplementary figure
methods <- sort(unique(imba_dge_merged$Method))
dge_method_plt <- function(imba_dge_df, method) {
  imba_dge_df_method <- imba_dge_df[which(imba_dge_df$Method %in% method)]
  ggplot(
    data = imba_dge_df_method, 
    aes(
      x = `Gene`,
      y = `Max rank`
    )
  ) +
    geom_boxplot(
      fill = "dodgerblue2",
      notch = FALSE,
      alpha = 0.8 
    ) + 
    facet_wrap(
      .~factor(type, levels = c("Control", "Downsampled", "Ablated")), 
      scales = "fixed"
    ) +
    scale_fill_manual( 
      breaks = c("Control", "Downsampled", "Ablated"),
      values = c("forestgreen", "darkorchid3", "firebrick2")
    ) +
    labs(
      title = method,
      fill = "Type",
      x = "Marker gene",
      y = "Maximum rank in differential expression analysis across clusters"
    ) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(factor(imba_dge_df$Gene)))) + 
    theme_few() +
    theme(plot.title = element_text(size = 24)) +
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16)) +
    theme(strip.text.x = element_text(size = 16)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(legend.title = element_text(size = 16)) +
    theme(legend.text = element_text(size = 14))
}
lapply(methods, function(x) {
  dge_method_plt(imba_dge_merged, x)
  ggsave(
    paste0(
      "outs/control/figures/07_pbmc_ds_ablate_",
      x,
      "_dge_rankings_no_liger.pdf"
    ),
    width = 14,
    height = 10,
    device = cairo_pdf
  )
})

### Fig 3B - Heatmap of marker gene perturbation scores across methods,
### for control, downsampled and ablated cases 

# Get the standard deviation of each marker gene, with method and status 
# factored in 
gene_rank_var_method_type <- imba_dge_merged %>% 
  group_by(Gene, Method, type) %>%
  summarize(var = sd(`Max rank`))

# Create three separate matrices for each of the three types (control, 
# downsampled, ablated)
gene_rank_var_method_type_ctrl <- gene_rank_var_method_type[
  which(gene_rank_var_method_type$type %in% "Control"), 
]
gene_rank_var_method_type_ctrl_long <- reshape2::dcast(
  gene_rank_var_method_type_ctrl,
  Gene~Method,
  value.var = "var"
)
gene_rank_var_method_type_ctrl_long_mat <- as.matrix(
  gene_rank_var_method_type_ctrl_long[, -1]
)
rownames(gene_rank_var_method_type_ctrl_long_mat) <- 
  gene_rank_var_method_type_ctrl_long[, 1]

gene_rank_var_method_type_ds <- gene_rank_var_method_type[
  which(gene_rank_var_method_type$type %in% "Downsampled"), 
]
gene_rank_var_method_type_ds_long <- reshape2::dcast(
  gene_rank_var_method_type_ds,
  Gene~Method,
  value.var = "var"
)
gene_rank_var_method_type_ds_long_mat <- as.matrix(
  gene_rank_var_method_type_ds_long[, -1]
)
rownames(gene_rank_var_method_type_ds_long_mat) <- 
  gene_rank_var_method_type_ds_long[, 1]

gene_rank_var_method_type_abla <- gene_rank_var_method_type[
  which(gene_rank_var_method_type$type %in% "Ablated"), 
]
gene_rank_var_method_type_abla_long <- reshape2::dcast(
  gene_rank_var_method_type_abla,
  Gene~Method,
  value.var = "var"
)
gene_rank_var_method_type_abla_long_mat <- as.matrix(
  gene_rank_var_method_type_abla_long[, -1]
)
rownames(gene_rank_var_method_type_abla_long_mat) <- 
  gene_rank_var_method_type_abla_long[, 1]

# Create three separate heatmaps for each of the types
# (use max 10 here as this indicates significant enough changes)
col_pert = circlize::colorRamp2(
  c(
    0,
    10
  ),
  c("white", "firebrick1")
)
ht1 <- Heatmap(
  gene_rank_var_method_type_ctrl_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Control",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  col = col_pert
)
ht2 <- Heatmap(
  gene_rank_var_method_type_ds_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Downsampled",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)
ht3 <- Heatmap(
  gene_rank_var_method_type_abla_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Ablated",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)

# Create a heatmap indicating which celltype is associated with each
# marker gene
base_marker_gene_dup_added <- base_marker_genes %>%
  group_by(`Top 10 marker genes (union across batches)`) %>%
  summarize(Celltype = paste0(unique(Celltype), collapse = ', ')) %>%
  as.data.frame
base_marker_gene_dup_added_mat <- as.matrix(
  base_marker_gene_dup_added[, -1]
)
rownames(base_marker_gene_dup_added_mat) <- base_marker_gene_dup_added[, 1] 

# Color palette 
unique_dup_celltypes_sorted <- sort(unique(base_marker_gene_dup_added$Celltype))
unique_kev_palette_vals <- kev_palette[1:length(unique_dup_celltypes_sorted)]
col_celltype = unique_kev_palette_vals
names(col_celltype) <- unique_dup_celltypes_sorted

# Celltype heatmap 
ht4 <- Heatmap(
  base_marker_gene_dup_added_mat, 
  name = "Celltype associated \nwith marker gene", 
  width = unit(1, "cm"),
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = TRUE,
  col = col_celltype
)

# Plot and save all of the heatmaps together 
marker_pert_hm <- ht1 + ht2 + ht3 + ht4
CairoPDF(
  "outs/control/figures/07_marker_gene_pert_pbmc_control_heatmap_no_liger.pdf",
  width = 14, 
  height = 6
)
draw(
  marker_pert_hm,
  column_title = "Integration method",
  column_title_side = "bottom",
  row_title_side = "right",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold")
)
dev.off()

# Plot one version without the celltype labels 
ht1 <- Heatmap(
  gene_rank_var_method_type_ctrl_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Control",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  col = col_pert
)
ht2 <- Heatmap(
  gene_rank_var_method_type_ds_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Downsampled",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)
ht3 <- Heatmap(
  gene_rank_var_method_type_abla_long_mat, 
  name = "Marker gene \nperturbation score", 
  width = unit(3, "cm"),
  column_title = "Ablated",
  row_title = "Marker gene",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 10),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"), 
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_heatmap_legend = FALSE,
  col = col_pert
)
marker_per_hm_no_celltype <- ht1 + ht2 + ht3
CairoPDF(
  "outs/control/figures/07_marker_gene_pert_pbmc_control_heatmap_no_ctype_no_liger.pdf",
  width = 14, 
  height = 6
)
draw(
  marker_per_hm_no_celltype,
  column_title = "Integration method",
  column_title_side = "bottom",
  row_title_side = "right",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
)
dev.off()

### Fig 3C - Concordance of DGE genes - top 10 most variable marker
### genes across methods 

# Pick 10 exemplary genes that show the highest variability - take top 10 with 
# the exclusion of the mitochondrial gene 
top_10_variable_dge <- gene_rank_variance_sorted$Gene[1:11]
top_10_variable_dge <- top_10_variable_dge[-6]

# Subset data for top 10 with the exclusion of mitochondrial gene and 
# plot results for ranking change with ablation and downsampling 
imba_dge_merged_top_10_var_genes <- imba_dge_merged[
  imba_dge_merged$Gene %in% top_10_variable_dge
]
ggplot(imba_dge_merged_top_10_var_genes, aes(x = Gene, y = `Max rank`)) +
  geom_boxplot(
    aes(
      fill = factor(type, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(.~Method, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Maximum rank in differential expression analysis across clusters"
  ) +
  coord_flip() +
  theme_few() +
  scale_x_discrete(limits = rev(top_10_variable_dge)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_dge_rankings_top_10_var_genes_no_liger.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Get stdev in rank of each gene, subset by type and method 
gene_rank_variance_grouped <- imba_dge_merged %>% 
  group_by(Gene, Method, type) %>%
  summarize(`Max rank stdev` = sd(`Max rank`)) %>%
  as.data.table

# Get the plots for each method, and the standard deviations of max
# rank of the genes
methods <- sort(unique(imba_dge_merged$Method))
dge_stdev_method_plt <- function(imba_dge_stdev_df, method) {
  imba_dge_stdev_df_sub <- imba_dge_stdev_df[
    which(imba_dge_stdev_df$Method %in% method)
  ]
  ggplot(
    data = imba_dge_stdev_df_sub, 
    aes(
      x = `Gene`,
      y = `Max rank stdev`
    ) 
  ) +
    geom_bar(
      aes(
        fill = factor(
          imba_dge_stdev_df_sub$type, 
          levels = c("Control", "Downsampled", "Ablated")
        )
      ),
      stat = "identity",
      alpha = 0.8 
    ) + 
    facet_wrap(
      .~factor(
        imba_dge_stdev_df_sub$type, 
        levels = c("Control", "Downsampled", "Ablated")
      ), 
      scales = "fixed"
    ) +
    scale_fill_manual( 
      breaks = c("Control", "Downsampled", "Ablated"),
      values = c("forestgreen", "darkorchid3", "firebrick2")
    ) +
    labs(
      fill = "Type",
      x = "Marker gene",
      y = "Standard deviation of maximum rank in differential expression"
    ) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(factor(imba_dge_stdev_df_sub$Gene)))) + 
    theme_few() +
    theme(axis.title.x = element_text(size = 16)) +
    theme(axis.title.y = element_text(size = 16)) +
    theme(strip.text.x = element_text(size = 16)) +
    theme(plot.title = element_text(size = 14)) +
    theme(axis.text.x = element_text(size = 14)) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(legend.title = element_text(size = 16)) +
    theme(legend.text = element_text(size = 14))
}
lapply(methods, function(x) {
  dge_stdev_method_plt(gene_rank_variance_grouped, x)
  ggsave(
    paste0(
      "outs/control/figures/07_pbmc_ds_ablate_",
      x,
      "_dge_rankings_method_type_stdev_no_liger.pdf"
    ),
    width = 14,
    height = 10,
    device = cairo_pdf
  )
})

# Plot the top 10 most variable genes in terms of rank (across all methods)
# based on type and method 
gene_rank_variance_grouped_top_10_sub <- gene_rank_variance_grouped[
  which(gene_rank_variance_grouped$Gene %in% top_10_variable_dge)
]
ggplot(
  data = gene_rank_variance_grouped_top_10_sub, 
  aes(
    x = `Gene`,
    y = `Max rank stdev`
  ) 
) +
  geom_bar(
    aes(
      fill = factor(
        gene_rank_variance_grouped_top_10_sub$type, 
        levels = c("Control", "Downsampled", "Ablated")
      )
    ),
    stat = "identity",
    alpha = 0.8,
    position = "dodge2"
  ) + 
  facet_wrap(
    .~Method,
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Standard deviation of maximum rank in differential expression"
  ) +
  coord_flip() +
  theme_few() +
  scale_x_discrete(limits = rev(top_10_variable_dge)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "_dge_rankings_method_type_stdev_top_10_no_liger.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Get stdev in rank of each gene, subset by type only (no method) 
gene_rank_variance_grouped_nomethod <- imba_dge_merged %>% 
  group_by(Gene, type) %>%
  summarize(`Max rank stdev` = sd(`Max rank`)) %>%
  as.data.table

# Plot the values of standard deviations of each gene across all methods 
ggplot(
  data = gene_rank_variance_grouped_nomethod, 
  aes(
    x = `Gene`,
    y = `Max rank stdev`
  ) 
) +
  geom_bar(
    aes(
      fill = factor(
        gene_rank_variance_grouped_nomethod$type, 
        levels = c("Control", "Downsampled", "Ablated")
      )
    ),
    stat = "identity",
    alpha = 0.8 
  ) + 
  facet_wrap(
    .~factor(
      gene_rank_variance_grouped_nomethod$type, 
      levels = c("Control", "Downsampled", "Ablated")
    ), 
    scales = "fixed"
  ) +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Marker gene",
    y = "Standard deviation of maximum rank in differential expression"
  ) +
  coord_flip() +
  scale_x_discrete(
    limits = rev(levels(factor(gene_rank_variance_grouped_nomethod$Gene)))
  ) + 
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "_dge_rankings_type_stdev_no_liger.pdf"
  ),
  width = 14,
  height = 10,
  device = cairo_pdf
)

# Plot the correspondence of the top 10 most highly variable marker genes with
# their respective celltypes
top_10_variable_dge_markersub <- base_marker_genes[
  which(
    base_marker_genes$`Top 10 marker genes (union across batches)` %in% 
      top_10_variable_dge
  )
]

# Order the results 
top_10_variable_dge_markersub <- top_10_variable_dge_markersub[
  match( 
    top_10_variable_dge,
    top_10_variable_dge_markersub$`Top 10 marker genes (union across batches)`
  )
]
colnames(top_10_variable_dge_markersub) <- c(
  "Associated celltype",
  "Gene"
)

# Add the standard deviation values across subsets overall 
top_10_variable_dge_markersub <- merge(
  top_10_variable_dge_markersub,
  gene_rank_variance_sorted,
  join = "left",
  by = c(
    "Gene"
  ),
  all.y = FALSE
)

# Order based on the standard deviation
top_10_variable_dge_markersub <- top_10_variable_dge_markersub[
  order(top_10_variable_dge_markersub$var, decreasing = )
]

# Plot the SankeyNetwork diagram, save and then convert saved html to pdf
nodes <- data.frame(
  name = c(
    as.character(top_10_variable_dge_markersub$Gene), 
    as.character(top_10_variable_dge_markersub$`Associated celltype`) %>% 
      unique()
  )
)
top_10_variable_dge_markersub$IDsource <- match(
  top_10_variable_dge_markersub$Gene, nodes$name)-1 
top_10_variable_dge_markersub$IDtarget <- match(
  top_10_variable_dge_markersub$`Associated celltype`, nodes$name)-1

top_10_var_marker_genes_sankey <- sankeyNetwork(
  Links = top_10_variable_dge_markersub,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "var",
  NodeID = "name",
  sinksRight = FALSE
)
saveNetwork(
  top_10_var_marker_genes_sankey,
  file = paste0(
    "outs/control/figures/07_pbmc_ds_ablate_",
    "_top_var_dges_with_assoc_celltypes_no_liger.html"
  ),
  selfcontained = TRUE
)

### Fig 3D) - Correlation of marker gene instability and the celltype
### downsampled - determine if this can impact the instability or 
### if there is a correlation

# Get stdev in rank of each gene, subset by type, method, downsampled celltype
# and marker(s) associated with celltype
gene_rank_variance_grouped_celltype_specific <- imba_dge_merged %>% 
  group_by(Gene, Method, type, `Downsampled celltypes`) %>%
  summarize(`Max rank stdev` = sd(`Max rank`)) %>%
  as.data.table

# Merge the DGE stdev summary stats with the marker data indicating
# which celltype each marker corresponds to - ALLOW A CARTESIAN PRODUCT
# HERE AS EACH MARKER MAY BE SPECIFIC TO MORE THAN ONE CELLTYPE. The cartesian
# product should be valid, as we are considering correlations between celltype-
# specific markers and change in DGE status
base_marker_genes_copy <- base_marker_genes
colnames(base_marker_genes_copy) <- c(
  "Associated celltype",
  "Gene"
)
gene_rank_variance_grouped_celltype_specific_marker <- merge(
  gene_rank_variance_grouped_celltype_specific,
  base_marker_genes_copy,
  by = c(
    "Gene"
  ),
  allow.cartesian = TRUE
)
gene_rank_variance_grouped_celltype_specific_marker <- distinct(
  gene_rank_variance_grouped_celltype_specific_marker
)

# Format celltype names for plotting 
gene_rank_variance_grouped_celltype_specific_marker$`Downsampled celltypes` <-
  plyr::mapvalues(
    gene_rank_variance_grouped_celltype_specific_marker$`Downsampled celltypes`,
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

gene_rank_variance_grouped_celltype_specific_marker$`Associated celltype` <-
  plyr::mapvalues(
    gene_rank_variance_grouped_celltype_specific_marker$`Associated celltype`,
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

# For first plot, collapse/summarize further by averaging across all of the 
# genes (mean - not median)
gene_rank_variance_grouped_celltype_specific_marker_mean <-
  gene_rank_variance_grouped_celltype_specific_marker %>%
  group_by(Method, type, `Downsampled celltypes`, `Associated celltype`) %>%
  summarize(`Mean max rank stdev` = mean(`Max rank stdev`)) %>%
  as.data.table

# Remove 'None' from this - only considering cases where the downsampled 
# celltype is equivalent to the associated celltype of the given marker 
gene_rank_variance_grouped_celltype_specific_marker_mean <-
  gene_rank_variance_grouped_celltype_specific_marker_mean[
    gene_rank_variance_grouped_celltype_specific_marker_mean$
      `Downsampled celltypes` %ni% "None"
  ]

# Plot heatmaps specific to the overall results of each method - first for 
# the downsampled results 
# MARKER GENE PERTURBATION SCORE HAS TO BE DEFINED IN RESULTS/METHODS
gene_rank_variance_grouped_celltype_specific_marker_mean_ds <- 
  gene_rank_variance_grouped_celltype_specific_marker_mean[
    gene_rank_variance_grouped_celltype_specific_marker_mean$type %in% 
      "Downsampled"
  ]

# Clip values to 10, as this is considered a fair "maximum acceptable" 
# perturbation
gene_rank_variance_grouped_celltype_specific_marker_mean_ds$
  `Mean max rank stdev clipped` <- pmin(
    10,
    gene_rank_variance_grouped_celltype_specific_marker_mean_ds$`Mean max rank stdev`
  )

ggplot(
  data = gene_rank_variance_grouped_celltype_specific_marker_mean_ds,
  aes(
    x = `Downsampled celltypes`,
    y = `Associated celltype`
  ) 
) + 
  geom_tile(
    aes(fill = `Mean max rank stdev clipped`)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "darkorchid3"
  ) +
  facet_wrap(.~Method, scales = "free") +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(
    size = 12, 
    angle = 90, 
    vjust = 1, 
    hjust = 1)
  ) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  labs(
    fill = "Average marker gene \nperturbation score",
    x = "Downsampled celltype",
    y = "Celltype associated with marker genes"
  )
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ds_only_",
    "_dge_rankings_celltype_marker_celltype_ds_compare_no_liger.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

# Plot heatmaps specific to the overall results of each method - now for 
# the ablated results 
# MARKER GENE PERTURBATION SCORE HAS TO BE DEFINED IN RESULTS/METHODS
gene_rank_variance_grouped_celltype_specific_marker_mean_ablated <- 
  gene_rank_variance_grouped_celltype_specific_marker_mean[
    gene_rank_variance_grouped_celltype_specific_marker_mean$type %in% 
      "Ablated"
  ]

# Clip values to 10, as this is considered a fair "maximum acceptable" 
# perturbation
gene_rank_variance_grouped_celltype_specific_marker_mean_ablated$
  `Mean max rank stdev clipped` <- pmin(
    10,
    gene_rank_variance_grouped_celltype_specific_marker_mean_ablated$`Mean max rank stdev`
  )

ggplot(
  data = gene_rank_variance_grouped_celltype_specific_marker_mean_ablated,
  aes(
    x = `Downsampled celltypes`,
    y = `Associated celltype`
  ) 
) + 
  geom_tile(
    aes(fill = `Mean max rank stdev clipped`)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "firebrick2"
  ) +
  facet_wrap(.~Method, scales = "free") +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(
    size = 12, 
    angle = 90, 
    vjust = 1, 
    hjust = 1)
  ) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  labs(
    fill = "Average marker gene \nperturbation score",
    x = "Ablated celltype",
    y = "Celltype associated with marker genes"
  )
ggsave(
  paste0(
    "outs/control/figures/07_pbmc_ablated_only_",
    "_dge_rankings_celltype_marker_celltype_ablated_compare_no_liger.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

### Fig 3E) - Impact of imbalance on annotation of data using an external
### reference

# Combine annotation score and imbalance dataframes
imba_anno_merged <- merge(
  imba_concat,
  anno_scores_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
) 
imba_anno_merged <- distinct(imba_anno_merged)

# Indicate which samples are controls and which are real runs
imba_anno_merged$type <- ifelse(
  imba_anno_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_anno_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format celltype names for plotting 
imba_anno_merged$`Downsampled celltypes` <- plyr::mapvalues(
  imba_anno_merged$`Downsampled celltypes`,
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

# Reformat to get the celltype-specific F1-scores for both L1 and L2 
# annotations 
imba_anno_merged_score_format <- imba_anno_merged[
  ,
  c(
    "Downsampled celltypes",
    "B cell", 
    "CD4 T cell", 
    "CD8 T cell", 
    "Monocyte_CD14", 
    "Monocyte_FCGR3A",
    "NK cell",
    "Overall balanced accuracy",
    "Overall F1-score",
    "Subset",
    "Score type",
    "type"
  ),
  with = FALSE
] 
imba_anno_merged_score_format <- melt(
  imba_anno_merged_score_format,
  id.vars = c(
    "Downsampled celltypes",
    "Overall balanced accuracy",
    "Overall F1-score",
    "Subset",
    "Score type",
    "type"
  ),
  value.name = "Score",
  variable.name = "Celltype scored"
) 
imba_anno_merged_score_format <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$`Score type` %in% "f1-score", 
]

# Format celltype names for plotting 
imba_anno_merged_score_format$`Celltype scored` <- plyr::mapvalues(
  imba_anno_merged_score_format$`Celltype scored`,
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

# Plot the celltype specific downsampling results, for both L1, and L2 
# annotations - F1 scores as a function of type, and celltype downsampled, 
# for each annotation
imba_anno_merged_score_format_l1 <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L1", 
]
imba_anno_merged_score_format_l2 <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L2", 
]

ggplot(data = imba_anno_merged_score_format_l1, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L1 annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_f1_scores_celltype_ds_no_liger.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

ggplot(data = imba_anno_merged_score_format_l2, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_jitter(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L2 annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_f1_scores_celltype_ds_no_liger.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

## Focus on T-cells here, plot the variability of T cell annotation
## accuracy within downsampled, ablated, and control results

# Get the variance of the F1-score, based on celltype scored and 
# downsampling status
imba_anno_merged_score_format_l1_sd <- imba_anno_merged_score_format_l1 %>% 
  group_by(`Celltype scored`, type) %>%
  summarize(F1_sd = sd(`Score`)) 
imba_anno_merged_score_format_l2_sd <- imba_anno_merged_score_format_l2 %>% 
  group_by(`Celltype scored`, type) %>%
  summarize(F1_sd = sd(`Score`)) 

# Plot these results, emphasizing high variance of F1-scores for T cells 
ggplot(data = imba_anno_merged_score_format_l1_sd, aes(
  x = `Celltype scored`,
  y = F1_sd
)) +
  geom_bar(
    aes(fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated"))),
    alpha = 0.8,
    stat = "identity",
    position = "dodge",
    width = 3
  ) + 
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype",
    y = "Standard deviation of L1 F1-scores"
  ) +
  theme_few() +
  coord_flip() +
  scale_x_discrete(
    limits = rev(imba_anno_merged_score_format_l1_sd$`Celltype scored`)
  ) + 
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_f1_scores_celltype_sdev_no_liger.pdf"
  ),
  width = 7,
  height = 5,
  device = cairo_pdf
)  

ggplot(data = imba_anno_merged_score_format_l2_sd, aes(
  x = `Celltype scored`,
  y = F1_sd
)) +
  geom_bar(
    aes(fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated"))),
    alpha = 0.8,
    stat = "identity",
    position = "dodge",
    width = 3
  ) + 
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype",
    y = "Standard deviation of L2 F1-scores"
  ) +
  theme_few() +
  coord_flip() +
  scale_x_discrete(
    limits = rev(imba_anno_merged_score_format_l2_sd$`Celltype scored`)
  ) + 
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_f1_scores_celltype_sdev_no_liger.pdf"
  ),
  width = 7,
  height = 5,
  device = cairo_pdf
)  

# Determine the celltypes that are being predicted for both CD4 and CD8 
# T cells, L1 and L2, for control, downsampled and ablated experiments 

# Merge all of the annotation results and imbalance results
imba_anno_merged_all <- merge(
  imba_concat,
  anno_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
)

# Indicate which samples are controls and which are real runs
imba_anno_merged_all$type <- ifelse(
  imba_anno_merged_all$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_anno_merged_all$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format celltype names for plotting
imba_anno_merged_all$`Downsampled celltypes` <- plyr::mapvalues(
  imba_anno_merged_all$`Downsampled celltypes`,
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
imba_anno_merged_all$`Real celltype` <- plyr::mapvalues(
  imba_anno_merged_all$`Real celltype`,
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

# Subset for L1 and L2 cases for only the T cells 
imba_anno_merged_all_tcell <- imba_anno_merged_all[
  imba_anno_merged_all$`Real celltype` %in% c("CD4+ T cell", "CD8+ T cell") 
]

# Plot the results for the celltypes that are most predicted for the T-cell
# subsets - for L1 and L2 
ggplot(data = imba_anno_merged_all_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L1`
)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Real celltype`) +
  labs(
    fill = "Predicted L1 \ncelltype",
    x = "Type",
    y = "Proportion of L1 predictions"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_t_cell_preds_barplot_no_liger.pdf"
  ),
  width = 10,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L2`
)) +
  geom_bar(position = "fill") +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Real celltype`) +
  labs(
    fill = "Predicted L2 \ncelltype",
    x = "Type",
    y = "Proportion of L2 predictions"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_t_cell_preds_barplot_no_liger.pdf"
  ),
  width = 10,
  height = 7,
  device = cairo_pdf
) 

# Plot the same results as above, separately this time for CD4 and CD8 T cells,
# and L1 and L2 predictions, but conditioned on the celltype that has been
# downsampled

# Subsample data for both celltypes
imba_anno_merged_all_cd4_tcell <- imba_anno_merged_all_tcell[
  imba_anno_merged_all_tcell$`Real celltype` %in% "CD4+ T cell"
]
imba_anno_merged_all_cd8_tcell <- imba_anno_merged_all_tcell[
  imba_anno_merged_all_tcell$`Real celltype` %in% "CD8+ T cell"
]

ggplot(data = imba_anno_merged_all_cd4_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L1`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L1 \ncelltype",
    x = "",
    y = "Proportion of L1 predictions for CD4+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_cd4t_cell_preds",
    "_celltype_ds_specific_barplots_no_liger.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_cd4_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L2`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L2 \ncelltype",
    x = "",
    y = "Proportion of L2 predictions for CD4+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_cd4t_cell_preds",
    "_celltype_ds_specific_barplots_no_liger.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_cd8_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L1`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L1 \ncelltype",
    x = "",
    y = "Proportion of L1 predictions for CD8+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_annotation_cd8t_cell_preds",
    "_celltype_ds_specific_barplots_no_liger.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 

ggplot(data = imba_anno_merged_all_cd8_tcell, aes(
  x = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
  fill = `Predicted L2`
)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(.~`Downsampled celltypes`, scales = "free_x") +
  labs(
    fill = "Predicted L2 \ncelltype",
    x = "",
    y = "Proportion of L2 predictions for CD8+ T cells"
  ) +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) 
ggsave(
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_annotation_cd8t_cell_preds",
    "_celltype_ds_specific_barplots_no_liger.pdf"
  ),
  width = 12,
  height = 7,
  device = cairo_pdf
) 


# Plot the celltype specific downsampling results, for the L1 and L2 baseline
# annotations (using just HVG KNN) - F1 scores as a function of type, and 
# celltype downsampled, for each annotation
imba_anno_merged_score_format_l1_baseline <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L1 baseline", 
]
imba_anno_merged_score_format_l2_baseline <- imba_anno_merged_score_format[
  imba_anno_merged_score_format$Subset %in% "L2 baseline", 
]

ggplot(data = imba_anno_merged_score_format_l1_baseline, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_point(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L1 baseline annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l1_baseline_annotation_f1_scores_celltype_ds_no_liger.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

ggplot(data = imba_anno_merged_score_format_l2_baseline, aes(
  x = `Downsampled celltypes`,
  y = Score
)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) + 
  geom_point(color = "black", size = 0.4, alpha = 0.7) +
  facet_wrap(.~`Celltype scored`, scales = "fixed") +
  scale_fill_manual( 
    breaks = c("Control", "Downsampled", "Ablated"),
    values = c("forestgreen", "darkorchid3", "firebrick2")
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "L2 baseline annotation accuracy (F1-score)"
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
  paste0(
    "outs/control/figures/",
    "07_pbmc_ds_ablate_l2_baseline_annotation_f1_scores_celltype_ds_no_liger.pdf"
  ),
  width = 16,
  height = 7,
  device = cairo_pdf
)  

