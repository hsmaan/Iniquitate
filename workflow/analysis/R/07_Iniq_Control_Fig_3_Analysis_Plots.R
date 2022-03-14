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

### Fig 3C) - concordance of DGE for marker genes before and after 
### downsampling 

# Merge together imbalanced and dge ranking datasets 
imba_dge_merged <- merge(
  imba_concat,
  dge_rank_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
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

# For each method, subset, plot and save the variability results 
methods <- sort(unique(imba_dge_merged$Method))
dge_method_plt <- function(imba_dge_df, method) {
  imba_dge_df_method <- imba_dge_df[which(imba_dge_df$Method %in% method)]
  ggplot(data = imba_dge_df_method, aes(
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
      fill = "Type",
      x = "Gene",
      y = "Maximum rank in differential expression analysis across clusters"
    ) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(factor(imba_dge_df$Gene)))) + 
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

# Pick 10 exemplary genes that show the highest variability 
gene_rank_variance <- imba_dge_merged %>% 
  group_by(Gene) %>%
  summarize(var = sd(`Max rank`))
gene_rank_variance_sorted <- gene_rank_variance[
  order(gene_rank_variance$var, decreasing = TRUE),
]

# Take top 10 with the exclusion of the mitochondrial gene 
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
    x = "Gene",
    y = "Maximum rank in differential expression analysis across clusters"
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
      x = "Gene",
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
    x = "Gene",
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