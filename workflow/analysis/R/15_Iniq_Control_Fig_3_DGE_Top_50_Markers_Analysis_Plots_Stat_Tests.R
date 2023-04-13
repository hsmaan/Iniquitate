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
library(EMT)

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

# Change to results dir for control data with top 50 marker genes
setwd("../../../results/control_dge_top_50/")

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

# Change to top level dir 
setwd("../../..")

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
    "outs/control/figures/15_pbmc_ds_ablate_",
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
    "outs/control/figures/15_pbmc_ds_ablate_",
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
      "outs/control/figures/15_pbmc_ds_ablate_",
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
# (use max 50 here as this indicates significant enough changes)
col_pert = circlize::colorRamp2(
  c(
    0,
    50
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
base_marker_genes$Celltype <- plyr::mapvalues(
  base_marker_genes$Celltype,
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
base_marker_gene_dup_added <- base_marker_genes %>%
  group_by(`Top 50 marker genes (union across batches)`) %>%
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
  name = "Cell-type associated \nwith marker gene", 
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
  "outs/control/figures/15_marker_gene_pert_pbmc_control_heatmap_no_liger.pdf",
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
  "outs/control/figures/15_marker_gene_pert_pbmc_control_heatmap_no_ctype_no_liger.pdf",
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
    "outs/control/figures/15_pbmc_ds_ablate_dge_rankings_top_10_var_genes_no_liger.pdf"
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
      "outs/control/figures/15_pbmc_ds_ablate_",
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
    "outs/control/figures/15_pbmc_ds_ablate_",
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
    "outs/control/figures/15_pbmc_ds_ablate_",
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
    "outs/control/figures/15_pbmc_ds_ablate_",
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

# Clip values to 50, as this is considered a fair "maximum acceptable" 
# perturbation (since we're using 50 marker genes per celltype)
gene_rank_variance_grouped_celltype_specific_marker_mean_ds$
  `Mean max rank stdev clipped` <- pmin(
    50,
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
    x = "Downsampled cell-type",
    y = "Cell-type associated with marker genes"
  )
ggsave(
  paste0(
    "outs/control/figures/15_pbmc_ds_only_",
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

# Clip values to 50, as this is considered a fair "maximum acceptable" 
# perturbation (since we're using 50 marker genes per celltype)
gene_rank_variance_grouped_celltype_specific_marker_mean_ablated$
  `Mean max rank stdev clipped` <- pmin(
    50,
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
    x = "Ablated cell-type",
    y = "Cell-type associated with marker genes"
  )
ggsave(
  paste0(
    "outs/control/figures/15_pbmc_ablated_only_",
    "_dge_rankings_celltype_marker_celltype_ablated_compare_no_liger.pdf"
  ),
  width = 14,
  height = 9,
  device = cairo_pdf
)

# Format marker genes to indicate which celltypes they belong to (due to 
# duplicates)
base_marker_gene_dup_added <- base_marker_genes %>%
  group_by(`Top 50 marker genes (union across batches)`) %>%
  summarize(Celltype = paste0(unique(Celltype), collapse = ', ')) %>%
  as.data.frame
colnames(base_marker_gene_dup_added) <- c(
  "marker_gene", "celltype"
)

### Statistical tests II - analysis of marker gene DGE rank change 

# Merge together imbalanced and dge data 
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

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_dge_unique_celltypes <- unique(imba_dge_merged$`Downsampled celltypes`)
imba_dge_unique_celltypes <- imba_dge_unique_celltypes[
  imba_dge_unique_celltypes %ni% "None"
]
imba_dge_none_celltype_len <- length(
  which(
    imba_dge_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_dge_none_celltype_draw <- sample(
  imba_dge_unique_celltypes,
  imba_dge_none_celltype_len,
  replace = TRUE
)
imba_dge_merged$`Downsampled celltypes`[
  which(
    imba_dge_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_dge_none_celltype_draw

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

# Format column names for an ANOVA test 
colnames(imba_dge_merged) <- plyr::mapvalues(
  colnames(imba_dge_merged),
  from = c(
    "Gene",
    "Max rank",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "gene",
    "max_rank",
    "method",
    "downsampled_celltypes"
  )
)

## Perform an overall ANOVA test for concordance of max rank of marker gene vs 
## type, after factoring in marker gnee, method and celltype downsampled - format 
## and save the result
dge_model_fit <- lm(
  as.formula(
    paste0(
      "max_rank", 
      "~",
      "gene+",
      "method+",
      "downsampled_celltypes+",
      "type"
    )
  ),
  data = imba_dge_merged
)
dge_anova_result <- anova(dge_model_fit, test = "F")
dge_anova_result_dt <- as.data.table(dge_anova_result, keep.rownames = TRUE)
colnames(dge_anova_result_dt)[1] <- "Covariate"
dge_anova_result_dt$dataset_name <- "PBMC 2 batch base balanced"
dge_anova_result_dt$metric <- "Max DGE rank marker gene"
dge_anova_result_dt$last_covariate <- "type"
fwrite(
  dge_anova_result_dt,
  paste0(
    "outs/control/results/",
    "15_pbmc_base_dge_rank_aov_results_ctrl_gene_method_celltype_ds_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## Perform this ANOVA test independantly for each marker gene, now, so removing
## the conditioning on the marker gene, but keeping the conditioning on the 
## method and downsampled celltypes 
marker_gene_list <- as.list(unique(imba_dge_merged$gene))

# Create ANOVA function to perform anova for each marker gene
anova_marker_gene <- function(
    dataset, 
    dataset_name,
    marker_gene
){
  # Subset dataset for the given marker
  dataset_marker_sub <- dataset[
    which(dataset$gene %in% marker_gene)
  ]
  
  # Fit ANOVA model
  model_fit <- lm(
    as.formula(
      paste0(
        "max_rank", 
        "~",
        "method+",
        "downsampled_celltypes+",
        "type"
      )
    ),
    data = dataset_marker_sub
  )
  anova_result <- anova(model_fit, test = "F")
  
  # Format results and return
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$marker_gene <- marker_gene
  anova_result_dt$metric <- "Max DGE rank"
  anova_result_dt$last_covariate <- "type"
  return(anova_result_dt)
}

# Iterate over marker genes, get results and concatenate and save
marker_aov_results <- lapply(marker_gene_list, function(x) {
  aov_summary_df <- anova_marker_gene(
    dataset = imba_dge_merged,
    dataset_name = "PBMC 2 batch base balanced",
    marker_gene = x
  )
})

# Concatenate results, correct for FDR, and save 
marker_aov_results_concat <- Reduce(rbind, marker_aov_results)
marker_aov_results_concat$`FDR_q` <- p.adjust(
  marker_aov_results_concat$`Pr(>F)`
)
fwrite(
  marker_aov_results_concat,
  paste0(
    "outs/control/results/",
    "15_pbmc_base_dge_specific_rank_aov_results_ctrl_method_celltype_ds_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## Get summary of linear model fit for each marker gene, after accounting for 
## method - determine which celltypes downsampled are most significantly 
## correlated for a given marker gene in this setup

# Create function that returns summary of a linear model for a 
# marker gene 
lm_marker_gene <- function(
    dataset, 
    dataset_name,
    marker_gene
){
  # Subset dataset for the given marker
  dataset_marker_sub <- dataset[
    which(dataset$gene %in% marker_gene)
  ]
  
  # Fit linear regression model
  model_fit <- lm(
    as.formula(
      paste0(
        "max_rank", 
        "~",
        "method+",
        "downsampled_celltypes"
      )
    ),
    data = dataset_marker_sub
  )
  
  # Format summary results - append relevant information
  model_fit_summary <- summary(model_fit)$coefficients
  model_fit_summary_dt <- as.data.table(
    as.data.frame(
      model_fit_summary
    ),
    keep.rownames = TRUE
  )
  colnames(model_fit_summary_dt)[1] <- "coefficient"
  model_fit_summary_dt$marker_gene <- marker_gene
  model_fit_summary_dt$metric <- "Max DGE rank"
  model_fit_summary_dt$dataset_name <- "PBMC 2 batch base balanced"
  model_fit_summary_dt_merged <- merge(
    model_fit_summary_dt,
    base_marker_gene_dup_added,
    by = "marker_gene"
  )
  
  # Indicate if the top `celltype downsampled` coefficient corresponds to the 
  # celltype the marker gene is indicative of 
  model_fit_summary_dt_merged_coeff <- model_fit_summary_dt_merged[
    grep("downsampled_celltypes", model_fit_summary_dt_merged$coefficient)
  ]
  top_downsampled_celltype_idx <- which.max(
    abs(
      model_fit_summary_dt_merged_coeff$Estimate
    )
  )
  top_downsampled_celltype_uf <- model_fit_summary_dt_merged_coeff$coefficient[
    top_downsampled_celltype_idx
  ]
  top_downsampled_celltype_formatted <- str_split_fixed(
    top_downsampled_celltype_uf,
    "downsampled_celltypes",
    2
  )[,2]
  model_fit_summary_dt_merged$top_ds_celltype <- top_downsampled_celltype_formatted
  
  marker_associated_celltypes <- base_marker_gene_dup_added[
    base_marker_gene_dup_added$marker_gene == marker_gene, 
  ]$celltype
  model_fit_summary_dt_merged$marker_assoc_celltypes <- marker_associated_celltypes
  model_fit_summary_dt_merged$same_ds_lm_celltype <- ifelse(
    grepl(
      top_downsampled_celltype_formatted,
      marker_associated_celltypes,
      ignore.case = TRUE
    ),
    "Yes",
    "No"
  )
  
  # Return summary datatable
  return(model_fit_summary_dt_merged)
}

# Iterate over marker genes and get results 
marker_lm_results <- lapply(marker_gene_list, function(x) {
  aov_summary_df <- lm_marker_gene(
    dataset = imba_dge_merged,
    dataset_name = "PBMC 2 batch base balanced",
    marker_gene = x
  )
})

# Concatenate and FDR correct the results 
marker_lm_summaries_concat <- Reduce(rbind, marker_lm_results)
marker_lm_summaries_concat$q_val <- p.adjust(
  marker_lm_summaries_concat$`Pr(>|t|)`
)

# Save results for coefficients and top associated celltypes 
fwrite(
  marker_lm_summaries_concat,
  paste0(
    "outs/control/results/",
    "15_pbmc_base_dge_rank_coeffs_celltype_ds_method_importance_ranks_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## Get enrichment of top coefficient for celltype downsampled being the 
## celltype associated with the marker gene being perturbed - perform
## multinomial test for each downsampled celltype against all others

# Create function that takes a subset of marker associated celltypes 
# and returns a multinomial sample test for top associated downsampled
# celltype coefficient value
marker_celltype_multinomial <- function(
    celltype_name,
    dataset
) {
  # Get list of all celltypes (no-mixed) and their length
  celltypes_all <- unique(base_marker_genes$Celltype)
  celltypes_all_len <- length(celltypes_all)
  
  # Subset data for the given celltype 
  data_sub <- dataset[dataset$marker_assoc_celltypes %in% celltype_name,]
  
  # Get counts vector for all the different celltypes 
  celltype_counts_vector <- table(data_sub$top_ds_celltype)
  
  # Add any missing celltypes to vector
  missing_celltypes <- celltypes_all[
    celltypes_all %ni% names(celltype_counts_vector)
  ]
  for (celltype in missing_celltypes) {
    celltype_counts_vector[celltype] <- 0
  }
  
  # Perform multinomial test for observed celltypes in top downsampled
  # results
  prob <- rep(1/celltypes_all_len, celltypes_all_len)
  observed <- as.vector(celltype_counts_vector)
  p_val <- multinomial.test(
    observed, 
    prob, 
    ntrial = 1e10,
    MonteCarlo = TRUE
  )$p.value
  print(p_val)
  
  # Create results dataframe, sort by celltype, and return
  res_df <- as.data.frame(as.array(celltype_counts_vector))
  colnames(res_df) <- c(
    "Top_celltype_downsampled_coefficient",
    "Observed_value_for_top_celltype_downsampled"
  )
  res_df$Marker_associated_celltype <- celltype_name
  res_df$multinomial_p_val <- p_val
  res_df$dataset_name <- "PBMC 2 batch base balanced"
  res_df <- res_df[order(res_df$Top_celltype_downsampled_coefficient), ]
  print(res_df)
  return(res_df)
}

# Iterate over all celltypes (no mixing) and get multinomial test results
## THIS WILL NOT CONSIDER MARKER GENES THAT ARE MARKERS TO MORE THAN
## ONE CELLTYPE 
celltypes_list <- as.list(unique(base_marker_genes$Celltype))
multinomial_results <- lapply(celltypes_list, function(x) {
  results <- marker_celltype_multinomial(x, marker_lm_summaries_concat)
})

# Concatenate, fdr correct, and save the results
multinomial_results_concat <- Reduce(rbind, multinomial_results)
multinomial_results_concat$multinomial_q_val <- p.adjust(
  multinomial_results_concat$multinomial_p_val
)
fwrite(
  multinomial_results_concat,
  paste0(
    "outs/control/results/",
    "15_pbmc_base_ds_dge_marker_celltypes_top_ds_celltype_",
    "coeff_multinom_tests_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)