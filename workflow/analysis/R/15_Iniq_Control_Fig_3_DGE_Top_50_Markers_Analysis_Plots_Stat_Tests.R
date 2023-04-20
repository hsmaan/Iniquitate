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

# Create figures and results directories if they don't exist
if (!dir.exists("outs/control_dge_50/figures")) {
  dir.create("outs/control_dge_50/figures", recursive = TRUE)
}
if (!dir.exists("outs/control_dge_50/results")) {
  dir.create("outs/control_dge_50/results", recursive = TRUE)
}

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
    "outs/control_dge_50/figures/15_pbmc_ds_ablate_",
    "marker_gene_max_ranks_no_liger.pdf"
  ),
  width = 14,
  height = 14,
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
    "outs/control_dge_50/figures/15_pbmc_ds_ablate_",
    "marker_gene_max_ranks_per_method_no_liger.pdf"
  ),
  width = 14,
  height = 14,
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
      "outs/control_dge_50/figures/15_pbmc_ds_ablate_",
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
  "outs/control_dge_50/figures/15_marker_gene_pert_pbmc_control_heatmap_no_liger.pdf",
  width = 14, 
  height = 14
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
  "outs/control_dge_50/figures/15_marker_gene_pert_pbmc_control_heatmap_no_ctype_no_liger.pdf",
  width = 12, 
  height = 16
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
    "outs/control_dge_50/results/",
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
    "outs/control_dge_50/results/",
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
    "outs/control_dge_50/results/",
    "15_pbmc_base_dge_rank_coeffs_celltype_ds_method_importance_ranks_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)