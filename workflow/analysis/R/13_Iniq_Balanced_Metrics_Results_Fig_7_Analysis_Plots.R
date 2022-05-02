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

# Change to outs dir for balanced metrics results
setwd("../../../outs/balanced_metrics/results/")


# Load 7A class and metrics results
bal_7A_cluster_df <- fread("04_7A_cluster_class_results.tsv", sep = "\t")
bal_7A_metrics_df <- fread("04_7A_metrics_results.tsv", sep = "\t")

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

# 