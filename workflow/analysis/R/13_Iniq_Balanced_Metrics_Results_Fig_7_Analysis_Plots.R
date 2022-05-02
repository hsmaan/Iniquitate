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

# Change to outs dir for balanced metrics results
setwd("../../../outs/balanced_metrics/results/")


# Load all class and metrics results
bal_7A_cluster_df <- fread("04_7A_cluster_class_results.tsv", sep = "\t")
bal_7A_metrics_df <- fread("04_7A_metrics_results.tsv", sep = "\t")

bal_7B_cluster_df <- fread("04_7B_cluster_class_results.tsv", sep = "\t")
bal_7B_metrics_df <- fread("04_7B_metrics_results.tsv", sep = "\t")

bal_7C_cluster_df <- fread("04_7C_cluster_celltype_imbal_results.tsv", sep = "\t")
bal_7C_metrics_df <- fread("04_7C_metrics_results.tsv", sep = "\t")

bal_7D_cluster_df <- fread("04_7D_cluster_celltype_imbal_results.tsv", sep = "\t")
bal_7D_metrics_df <- fread("04_7D_metrics_results.tsv", sep = "\t")

bal_7E_cluster_df <- fread("04_7E_cluster_celltype_imbal_results.tsv", sep = "\t")
bal_7E_metrics_df <- fread("04_7E_metrics_results.tsv", sep = "\t")

# Change to top level dir
setwd("../../../")

# Make balanced metrics figures directory if it doesn't exist
if (!dir.exists("outs/balanced_metrics/figures")) {
  dir.create("outs/balanced_metrics/figures", recursive = TRUE)
}

### Fig 7A Analysis - First use case on simulated data 

# Plot the class coordinate results 
ggplot(data = bal_7A_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = class)) +
  labs(
    color = "Class",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_few() +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(
  "outs/balanced_metrics/figures/13_7A_trial_class_coordinates.pdf",
  width = 6,
  height = 6
)

# Plot the cluster coordinate results 
ggplot(data = bal_7A_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(kmeans))) +
  labs(
    color = "Cluster",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_few() +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(
  "outs/balanced_metrics/figures/13_7A_trial_cluster_coordinates.pdf",
  width = 6,
  height = 6
)

# Plot the comparisons of the balanced and imbalanced metrics for 7A
bal_7A_metrics_df$Metric_bare <- str_split_fixed(
  bal_7A_metrics_df$Metric,
  pattern = " ",
  n = 2
)[, 1]
ggplot(
  data = bal_7A_metrics_df, 
  aes(
    x = factor(
      bal_7A_metrics_df$Metric_bare,
      levels = rev(unique(bal_7A_metrics_df$Metric_bare))
    ), 
    y = Value
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = factor(
        bal_7A_metrics_df$Type, 
        levels = c("balanced", "imbalanced")
      )
    ), 
    position = position_dodge2()
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_few() +
  labs(
    fill = "Metric type",
    x = "Metric",
    y = "Value"
  ) +
  ylim(c(0, 1)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) +
  coord_flip() 
ggsave(
  "outs/balanced_metrics/figures/13_7A_metrics_barplot.pdf",
  width = 8,
  height = 8,
  device = cairo_pdf
)

### Fig 7B Analysis - Second use case on simulated data  

# Plot the class coordinate results 
ggplot(data = bal_7B_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = class)) +
  labs(
    color = "Class",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_few() +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(
  "outs/balanced_metrics/figures/13_7B_trial_class_coordinates.pdf",
  width = 6,
  height = 6
)

# Plot the cluster coordinate results 
ggplot(data = bal_7B_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(kmeans))) +
  labs(
    color = "Cluster",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_few() +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(
  "outs/balanced_metrics/figures/13_7B_trial_cluster_coordinates.pdf",
  width = 6,
  height = 6
)

# Plot the comparisons of the balanced and imbalanced metrics for 7B
bal_7B_metrics_df$Metric_bare <- str_split_fixed(
  bal_7B_metrics_df$Metric,
  pattern = " ",
  n = 2
)[, 1]
ggplot(
  data = bal_7B_metrics_df, 
  aes(
    x = factor(
      bal_7B_metrics_df$Metric_bare,
      levels = rev(unique(bal_7B_metrics_df$Metric_bare))
    ), 
    y = Value
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = factor(
        bal_7B_metrics_df$Type, 
        levels = c("balanced", "imbalanced")
      )
    ), 
    position = position_dodge2()
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_few() +
  labs(
    fill = "Metric type",
    x = "Metric",
    y = "Value"
  ) +
  ylim(c(0, 1)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) +
  coord_flip() 
ggsave(
  "outs/balanced_metrics/figures/13_7B_metrics_barplot.pdf",
  width = 8,
  height = 8,
  device = cairo_pdf
)

### Fig 7C Analysis - First use case on single-cell data 
palette_7C_celltype <- kev_palette[1:length(unique(bal_7C_cluster_df$celltype))]
palette_7C_cluster <- kev_palette[1:length(unique(bal_7C_cluster_df$cluster))]

# Plot the celltype coordinate results 
ggplot(data = bal_7C_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = celltype)) +
  labs(
    color = "Celltype",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = palette_7C_celltype) +
  theme_few() +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(
  "outs/balanced_metrics/figures/13_7C_trial_celltype_coordinates.pdf",
  width = 7,
  height = 7
)

# Plot the cluster coordinate results 
ggplot(data = bal_7C_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(cluster))) +
  labs(
    color = "Cluster",
    x = "x1",
    y = "x2"
  ) +
  scale_color_manual(values = palette_7C_cluster) +
  theme_few() +
  theme_classic() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) + 
  guides(colour = guide_legend(override.aes = list(size=3)))
ggsave(
  "outs/balanced_metrics/figures/13_7C_trial_cluster_coordinates.pdf",
  width = 7,
  height = 7
)

# Plot the comparisons of the balanced and imbalanced metrics for 7C
bal_7C_metrics_df$Metric_bare <- str_split_fixed(
  bal_7C_metrics_df$Metric,
  pattern = " ",
  n = 2
)[, 1]
ggplot(
  data = bal_7C_metrics_df, 
  aes(
    x = factor(
      bal_7C_metrics_df$Metric_bare,
      levels = rev(unique(bal_7C_metrics_df$Metric_bare))
    ), 
    y = Value
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = factor(
        bal_7C_metrics_df$Type, 
        levels = c("balanced", "imbalanced")
      )
    ), 
    position = position_dodge2()
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_few() +
  labs(
    fill = "Metric type",
    x = "Metric",
    y = "Value"
  ) +
  ylim(c(0, 1)) +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1) +
  coord_flip() 
ggsave(
  "outs/balanced_metrics/figures/13_7C_metrics_barplot.pdf",
  width = 8,
  height = 8,
  device = cairo_pdf
)


