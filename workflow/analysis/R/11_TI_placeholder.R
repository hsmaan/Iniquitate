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

# Load in and concatenate TI inference results files
setwd("../../../results/lowcap_modified/paga_ti_scores/")
ti_scores_files <- list.files()
ti_scores_loaded <- lapply(ti_scores_files, fread)
ti_scores_concat <- Reduce(rbind, ti_scores_loaded)

# Load in and concatenate imbalance summary files 
setwd("../paga_imbalance_summaries")
ti_imba_files <- list.files()
ti_imba_loaded <- lapply(ti_imba_files, fread)
ti_imba_concat <- Reduce(rbind, ti_imba_loaded)
gc()

# Change to top level dir
setwd("../../../")

# Merge together the imbalance summary and TI results 
ti_imba_scores_merged <- merge(
  ti_scores_concat,
  ti_imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Indicate which panels are control and which ones are ablations or downsampling
ti_imba_scores_merged$type <- ifelse(
  ti_imba_scores_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    ti_imba_scores_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Plot the comparisons based on celltypes downsampled, method utilized and 
# pearson correlations 
ggplot(
  data = ti_imba_scores_merged, 
  aes(
    x = factor(
      ti_imba_scores_merged$`Downsampled celltypes`,
      levels = rev(
        c(
          unique(ti_imba_scores_merged$`Downsampled celltypes`)
        )
      )
    ),
    y = `Pearson correlations`
  )
) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(
    .~`Method`
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "Pearson correlation"
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
  theme(legend.text = element_text(size = 16)) +
  coord_flip()
ggsave(
  paste0(
    "outs/lowcap_modified/figures/",
    "11_chondrocyte_ti_pearson_corr_ds_celltype_method.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Plot the comparisons based on celltypes downsampled, method utilized and 
# spearman correlations
ggplot(
  data = ti_imba_scores_merged, 
  aes(
    x = factor(
      ti_imba_scores_merged$`Downsampled celltypes`,
      levels = rev(
        c(
          unique(ti_imba_scores_merged$`Downsampled celltypes`)
        )
      )
    ),
    y = `Spearman correlations`
  )
) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(
    .~`Method`
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "Spearman correlation"
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
  theme(legend.text = element_text(size = 16)) +
  coord_flip()
ggsave(
  paste0(
    "outs/lowcap_modified/figures/",
    "11_chondrocyte_ti_spearman_corr_ds_celltype_method.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

# Plot the comparisons based on celltypes downsampled, method utilized and 
# kendall correlations
ggplot(
  data = ti_imba_scores_merged, 
  aes(
    x = factor(
      ti_imba_scores_merged$`Downsampled celltypes`,
      levels = rev(
        c(
          unique(ti_imba_scores_merged$`Downsampled celltypes`)
        )
      )
    ),
    y = `Kendall correlations`
  )
) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  facet_wrap(
    .~`Method`
  ) +
  labs(
    fill = "Type",
    x = "Celltype downsampled",
    y = "Kendall correlation"
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
  theme(legend.text = element_text(size = 16)) +
  coord_flip()
ggsave(
  paste0(
    "outs/lowcap_modified/figures/",
    "11_chondrocyte_ti_kendall_corr_ds_celltype_method.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)
