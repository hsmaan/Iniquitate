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
setwd("../../../results/control_ti_only/paga_ti_scores/")
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

# Create results directory if it doesn't exist
if (!dir.exists("outs/control_ti_only/figures/")) {
  dir.create("outs/control_ti_only/figures/", recursive = TRUE)
}

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
    "outs/control_ti_only/figures/",
    "19_chondrocyte_ti_pearson_corr_ds_celltype_method.pdf"
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
    "outs/control_ti_only/figures/",
    "19_chondrocyte_ti_spearman_corr_ds_celltype_method.pdf"
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
    "outs/control_ti_only/figures/",
    "19_chondrocyte_ti_kendall_corr_ds_celltype_method.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

#############################################################################
### REDO ALL FIGURES ABOVE WITHOUT LIGER ###
#############################################################################
rm(list = ls())
gc()

`%ni%` <- Negate(`%in%`)
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
setwd("results/control_ti_only/paga_ti_scores/")
ti_scores_files <- list.files()
ti_scores_loaded <- lapply(ti_scores_files, fread)
ti_scores_concat <- Reduce(rbind, ti_scores_loaded)
ti_scores_concat <- ti_scores_concat[
  ti_scores_concat$Method != "liger"
]

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

# Remove the underscores from the celltype names
ti_imba_scores_merged$`Downsampled celltypes` <- gsub(
  "_", " ", ti_imba_scores_merged$`Downsampled celltypes`
)

# Change 'Intermediate Mesoderm' to 'Intermediate mesoderm'
ti_imba_scores_merged$`Downsampled celltypes` <- gsub(
  "Intermediate Mesoderm", "Intermediate mesoderm", 
  ti_imba_scores_merged$`Downsampled celltypes`
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
    "outs/control_ti_only/figures/",
    "19_chondrocyte_ti_pearson_corr_ds_celltype_method_no_liger.pdf"
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
    "outs/control_ti_only/figures/",
    "19_chondrocyte_ti_spearman_corr_ds_celltype_method_no_liger.pdf"
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
    "outs/control_ti_only/figures/",
    "19_chondrocyte_ti_kendall_corr_ds_celltype_method_no_liger.pdf"
  ),
  width = 14,
  height = 8,
  device = cairo_pdf
)

### Plot the embeddings for the Cao et al. data from both batches

# Change to Cao et al balanced data dir  
setwd("resources/h5ad_files/int_paga_datasets/cao_organ_dev_sublin_balanced_no_meta_2_batch/")

# Convert balanced data files to h5seurat format
Convert(
  "cao_dev_lite_balanced_no_meta_batch_1.h5ad", 
  dest = "h5seurat", 
  overwrite = TRUE
)
Convert(
  "cao_dev_lite_balanced_no_meta_batch_2.h5ad",
  dest = "h5seurat",
  overwrite = TRUE
)

# Load h5ad files for both balanced batches 
dev_data_1 <- SeuratDisk::LoadH5Seurat("cao_dev_lite_balanced_no_meta_batch_1.h5seurat")
dev_data_2 <- SeuratDisk::LoadH5Seurat("cao_dev_lite_balanced_no_meta_batch_2.h5seurat")

# Change to top level dir 
setwd("../../../../")

### Plotting of balanced dataset  

# Process both datasets independantly and together
dev_combined <- merge(dev_data_1, dev_data_2)

dev_data_1 <- NormalizeData(dev_data_1)
dev_data_1 <- FindVariableFeatures(
  dev_data_1, selection.method = "vst", nfeatures = 2000
)
all.genes <- rownames(dev_data_1)
dev_data_1 <- ScaleData(dev_data_1, features = all.genes)
dev_data_1 <- RunPCA(dev_data_1, features = VariableFeatures(object = dev_data_1))
dev_data_1 <- FindNeighbors(dev_data_1, dims = 1:20)
dev_data_1 <- FindClusters(dev_data_1, resolution = 0.5)
dev_data_1 <- RunUMAP(dev_data_1, dims = 1:20)

dev_data_2 <- NormalizeData(dev_data_2)
dev_data_2 <- FindVariableFeatures(
  dev_data_2, selection.method = "vst", nfeatures = 2000
)
all.genes <- rownames(dev_data_2)
dev_data_2 <- ScaleData(dev_data_2, features = all.genes)
dev_data_2 <- RunPCA(dev_data_2, features = VariableFeatures(object = dev_data_2))
dev_data_2 <- FindNeighbors(dev_data_2, dims = 1:20)
dev_data_2 <- FindClusters(dev_data_2, resolution = 0.5)
dev_data_2 <- RunUMAP(dev_data_2, dims = 1:20)

dev_combined <- NormalizeData(dev_combined)
dev_combined <- FindVariableFeatures(
  dev_combined, selection.method = "vst", nfeatures = 2000
)
all.genes <- rownames(dev_combined)
dev_combined <- ScaleData(dev_combined, features = all.genes)
dev_combined <- RunPCA(dev_combined, features = VariableFeatures(
  object = dev_combined
))
dev_combined <- FindNeighbors(dev_combined, dims = 1:20)
dev_combined <- FindClusters(dev_combined, resolution = 0.5)
dev_combined <- RunUMAP(dev_combined, dims = 1:20)

### Create plots of celltype embeddings for the two datasets combined

# Substitute the batch names with more informative names 
dev_combined$batch <- ifelse(
  dev_combined$batch == "batch_1", 
  "Batch 1 (Day 1)", 
  "Batch 2 (Day 5)"
)

# Plot combined data results as they are, for both celltype and batch 
DimPlot(dev_combined, reduction = "umap", group.by = "celltype") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Dark2") +  
  theme(aspect.ratio = 1)
ggsave(
  "outs/control_ti_only/figures/19_dev_balanced_combined_celltypes.pdf",
  width = 6, 
  height = 6
)
DimPlot(dev_combined, reduction = "umap", group.by = "batch") +
  theme(plot.title = element_blank()) + 
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  scale_color_brewer(palette = "Set1") + 
  theme(aspect.ratio = 1)
ggsave(
  "outs/control_ti_only/figures/19_dev_balanced_combined_batch.pdf",
  width = 6,
  height = 6
)