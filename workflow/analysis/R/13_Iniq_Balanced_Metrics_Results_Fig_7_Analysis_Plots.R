library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(ggpubr)
library(ggbump)
library(dotwhisker)
library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)
library(cowplot)

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
p_7a_1 <- ggplot(data = bal_7A_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = class)) +
  labs(
    color = "Class",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Set1") +
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
p_7a_1
ggsave(
  "outs/balanced_metrics/figures/13_7A_trial_class_coordinates.pdf",
  width = 6,
  height = 6
)

# Plot the cluster coordinate results 
p_7a_2 <- ggplot(data = bal_7A_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(kmeans))) +
  labs(
    color = "Cluster",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
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
p_7a_2
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
p_7a_3 <- ggplot(
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
  theme_classic() +
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
p_7a_3
ggsave(
  "outs/balanced_metrics/figures/13_7A_metrics_barplot.pdf",
  width = 8,
  height = 8,
  device = cairo_pdf
)

# Plot them together using cowplot 
all_7a_plots <- plot_grid(
  p_7a_1,
  p_7a_2,
  p_7a_3,
  labels = "A",
  label_size = 12,
  rel_widths = c(1, 1, 3),
  ncol = 1
)
save_plot(
  "outs/balanced_metrics/figures/13_7A_all_plots.pdf",
  all_7a_plots,
  base_asp = 1,
  base_height = 8,
  base_width = 12
)

### Fig 7B Analysis - Second use case on simulated data  

# Plot the class coordinate results 
p_7b_1 <- ggplot(data = bal_7B_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = class)) +
  labs(
    color = "Class",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
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
p_7b_1
ggsave(
  "outs/balanced_metrics/figures/13_7B_trial_class_coordinates.pdf",
  width = 6,
  height = 6
)

# Plot the cluster coordinate results 
p_7b_2 <- ggplot(data = bal_7B_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(kmeans))) +
  labs(
    color = "Cluster",
    x = "x1",
    y = "x2"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_classic() +
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
p_7b_2
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
p_7b_3 <- ggplot(
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
  theme_classic() +
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
p_7b_3
ggsave(
  "outs/balanced_metrics/figures/13_7B_metrics_barplot.pdf",
  width = 8,
  height = 8,
  device = cairo_pdf
)

# Plot all 7B plots together 
all_7b_plots <- plot_grid(
  p_7b_1,
  p_7b_2,
  p_7b_3,
  labels = "B",
  label_size = 12,
  rel_widths = c(1, 1, 3),
  ncol = 1
)
save_plot(
  "outs/balanced_metrics/figures/13_7B_all_plots.pdf",
  all_7b_plots,
  base_asp = 1,
  base_height = 8,
  base_width = 12
)



### Fig 7C Analysis - First use case on single-cell data 
palette_7C_celltype <- kev_palette[1:length(unique(bal_7C_cluster_df$celltype))]
palette_7C_cluster <- kev_palette[1:length(unique(bal_7C_cluster_df$cluster))]

# Plot the celltype coordinate results 
ggplot(data = bal_7C_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = celltype), size = 1) +
  labs(
    color = "Celltype",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = palette_7C_celltype) +
  theme_classic() +
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
  geom_point(aes(color = factor(cluster)), size = 1) +
  labs(
    color = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = palette_7C_cluster) +
  theme_classic() +
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
# Focus on just the ARI and AMI results 
bal_7C_metrics_df$Metric_bare <- str_split_fixed(
  bal_7C_metrics_df$Metric,
  pattern = " ",
  n = 2
)[, 1]
bal_7C_metrics_df_sub <- bal_7C_metrics_df[
  bal_7C_metrics_df$Metric_bar %in% c("ARI", "AMI")
]
ggplot(
  data = bal_7C_metrics_df_sub, 
  aes(
    x = factor(
      bal_7C_metrics_df_sub$Metric_bare,
      levels = unique(bal_7C_metrics_df_sub$Metric_bare)
    ), 
    y = Value
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = factor(
        bal_7C_metrics_df_sub$Type, 
        levels = c("imbalanced", "balanced")
      )
    ), 
    position = position_dodge2()
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
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
  theme(aspect.ratio = 1) 
ggsave(
  "outs/balanced_metrics/figures/13_7C_metrics_barplot.pdf",
  width = 6,
  height = 6,
  device = cairo_pdf
)

### Fig 7D Analysis - Second use case on single-cell data 

# Plot the celltype coordinate results 
p_7d_1 <- ggplot(data = bal_7D_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = celltype), size = 1) +
  labs(
    color = "Celltype",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") +
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
p_7d_1
ggsave(
  "outs/balanced_metrics/figures/13_7D_trial_celltype_coordinates.pdf",
  width = 7,
  height = 7
)

# Plot the cluster coordinate results 
p_7d_2 <- ggplot(data = bal_7D_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(cluster)), size = 1) +
  labs(
    color = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set2") +
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
p_7d_2
ggsave(
  "outs/balanced_metrics/figures/13_7D_trial_cluster_coordinates.pdf",
  width = 7,
  height = 7
)

# Plot the comparisons of the balanced and imbalanced metrics for 7D
# Focus on just the ARI and Homogeneity results 
bal_7D_metrics_df$Metric_bare <- str_split_fixed(
  bal_7D_metrics_df$Metric,
  pattern = " ",
  n = 2
)[, 1]
bal_7D_metrics_df_sub <- bal_7D_metrics_df[
  bal_7D_metrics_df$Metric_bar %in% c("ARI", "Homogeneity")
]
p_7d_3 <- ggplot(
  data = bal_7D_metrics_df_sub, 
  aes(
    x = factor(
      bal_7D_metrics_df_sub$Metric_bare,
      levels = unique(bal_7D_metrics_df_sub$Metric_bare)
    ), 
    y = Value
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = factor(
        bal_7D_metrics_df_sub$Type, 
        levels = c("imbalanced", "balanced")
      )
    ), 
    position = position_dodge2()
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_classic() +
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
  coord_fixed() +
  theme(aspect.ratio = 1) 
p_7d_3
ggsave(
  "outs/balanced_metrics/figures/13_7D_metrics_barplot.pdf",
  width = 6,
  height = 6,
  device = cairo_pdf
)

# Plot all fig 7d results together 
all_7d_plots <- plot_grid(
  p_7d_1,
  p_7d_2,
  p_7d_3,
  labels = "D",
  label_size = 12,
  rel_widths = c(2.2, 1.5, 2),
  ncol = 3
)
save_plot(
  "outs/balanced_metrics/figures/13_7D_all_plots.pdf",
  all_7d_plots,
  base_asp = 1,
  base_height = 4,
  base_width = 14
)

### Fig 7E Analysis - Third use case on single-cell data - integration based

# Plot the celltype coordinate results per method 
p_7e_1 <- ggplot(data = bal_7E_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = celltype), size = 1) +
  facet_wrap(.~Subset, scales = "free") + 
  labs(
    color = "Celltype",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme_few() +
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
p_7e_1
ggsave(
  "outs/balanced_metrics/figures/13_7E_trial_celltype_coordinates_facet.pdf",
  width = 10,
  height = 8
)

# Plot the cluster coordinate results per method 
palette_7E_cluster <- kev_palette[1:length(unique(bal_7E_cluster_df$cluster))]
p_7e_2 <- ggplot(data = bal_7E_cluster_df, aes(x = x, y = y)) +
  geom_point(aes(color = factor(cluster)), size = 1) +
  facet_wrap(.~Subset, scales = "free") + 
  labs(
    color = "Cluster",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_manual(values = palette_7E_cluster) +
  theme_few() +
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
p_7e_2
ggsave(
  "outs/balanced_metrics/figures/13_7E_trial_cluster_coordinates.pdf",
  width = 10,
  height = 8
)

# Aggregate scores for both balanced and imbalanced metrics and rank the methods
bal_7E_metrics_df_bal_sub <- bal_7E_metrics_df[
  bal_7E_metrics_df$Type == "balanced"
]
bal_7E_metrics_df_imbal_sub <- bal_7E_metrics_df[
  bal_7E_metrics_df$Type == "imbalanced"
]

bal_7E_score_agg_imbal <- bal_7E_metrics_df_imbal_sub %>% 
  group_by(Subset) %>%
  summarise("Metric_avg" = mean(Value)) %>%
  arrange(desc(Metric_avg)) %>%
  as.data.frame()

bal_7E_score_agg_bal <- bal_7E_metrics_df_bal_sub %>% 
  group_by(Subset) %>%
  summarise("Metric_avg" = mean(Value)) %>%
  arrange(desc(Metric_avg)) %>%
  as.data.frame()

# Format and plot comparison of the mean results and rankings 
bal_7E_score_agg_imbal$Type <- "imbalanced"
bal_7E_score_agg_bal$Type <- "balanced"

bal_7E_score_agg_imbal$Rank <- seq(1, 5, 1)
bal_7E_score_agg_bal$Rank <- seq(1, 5, 1)

bal_7E_score_agg_concat <- rbind(
  bal_7E_score_agg_imbal,
  bal_7E_score_agg_bal
)

p_7e_3 <- ggplot(
  data = bal_7E_score_agg_imbal, 
  aes(
    x = factor(
      bal_7E_score_agg_imbal$Subset,
      levels = rev(bal_7E_score_agg_imbal$Subset)
    ), 
    y = Metric_avg
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = Rank
    )
  ) +
  scale_fill_viridis_c() +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_classic() +
  labs(
    fill = "Rank",
    x = "Method",
    y = "Average metric score - imbalanced metrics"
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
p_7e_3
ggsave(
  "outs/balanced_metrics/figures/13_7E_trial_avg_metric_score_imbal.pdf",
  width = 7,
  height = 7
)

p_7e_4 <- ggplot(
  data = bal_7E_score_agg_bal, 
  aes(
    x = factor(
      bal_7E_score_agg_bal$Subset,
      levels = rev(bal_7E_score_agg_bal$Subset)
    ), 
    y = Metric_avg
  )
) +
  geom_bar(
    stat = "identity", 
    aes(
      fill = Rank
    )
  ) +
  scale_fill_viridis_c() + 
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_classic() +
  labs(
    fill = "Rank",
    x = "Method",
    y = "Average metric score - balanced metrics"
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
p_7e_4
ggsave(
  "outs/balanced_metrics/figures/13_7E_trial_avg_metric_score_bal.pdf",
  width = 7,
  height = 7
)

# Plot all the figure 7E plots together 
p_7e_all <- plot_grid(
  p_7e_1,
  p_7e_3,
  p_7e_2,
  p_7e_4,
  labels = "E",
  label_size = 12,
  rel_widths = c(1, 0.5, 1, 0.5),
  ncol = 2,
  nrow = 2
)
p_7e_all
save_plot(
  "outs/balanced_metrics/figures/13_7E_all_plots.pdf",
  p_7e_all,
  base_asp = 1,
  base_height = 12,
  base_width = 16
)
