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

# Change to results dir for uamp results data
setwd("../../../results/umap/")

# Load color palette 
kev_palette <- c("dodgerblue2", 
                 "#E31A1C",
                 "green4",
                 "#6A3D9A", 
                 "#FF7F00", 
                 "black",
                 "gold1",
                 "skyblue2",
                 "#FB9A99", 
                 "palegreen2",
                 "#CAB2D6", 
                 "#FDBF6F", 
                 "gray70", 
                 "khaki2",
                 "maroon",
                 "orchid1",
                 "deeppink1",
                 "blue1",
                 "steelblue4",
                 "darkturquoise",
                 "green1",
                 "yellow4",
                 "yellow3",
                 "darkorange4",
                 "brown")

##### Analysis of PBMC 2 batch balanced data - baseline #####

# Load in the umap plot results 
setwd("umap_plots/")
umap_files <- list.files()
umap_files <- grep(
  ".tsv",
  umap_files,
  value = TRUE
)
umap_files <- grep(
  "pbmc_2_batch_base_balanced",
  umap_files,
  value = TRUE
)
umap_loaded <- lapply(umap_files, fread)
umap_names <- str_split_fixed(umap_files, fixed(".tsv"), 2)[,1]
names(umap_loaded) <- umap_names

setwd("../../..")

# Create directory for umap results if it doesn't exist
if (!dir.exists("outs/umap/results")) {
  dir.create("outs/umap/results", recursive = TRUE)
}
if (!dir.exists("outs/umap/figures")) {
  dir.create("outs/umap/figures")
}

# Create function to loop over the umap files and return the results 
umap_plot <- function(df, save_prefix) {
  # Format celltype names
  df$Clustering <- plyr::mapvalues(
    df$Clustering,
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
  
  # Format batch names 
  df$Clustering <- plyr::mapvalues(
    df$Clustering,
    from = c(
      "batch_1",
      "batch_2"
    ),
    to = c(
      "Batch 1",
      "Batch 2"
    )
  )
  
  unique_cluster_len <- length(unique(df$Clustering))
  if (unique_cluster_len > 8) {
    ggplot(data = df, aes(x = `UMAP 1`, y = `UMAP 2`)) +
      geom_point(
        aes(
          color = factor(
            as.numeric(Clustering),
            levels = sort(as.numeric(unique(df$Clustering)))
          )
        ),
        size = 0.25
      ) +
      facet_wrap(
        .~Subset, 
        scales = "free"
      ) +
      labs(
        color = "",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      scale_color_manual(
        name = "",
        values = kev_palette[1:unique_cluster_len]
      ) + 
      guides(color = guide_legend(override.aes = list(size=2))) + 
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
      paste0(
        "outs/umap/figures/",
        save_prefix,
        ".pdf"
      ),
      width = 16,
      height = 8,
      device = cairo_pdf
    )
  } else {
    if (any(grepl("Batch", df$Clustering))) {
      pal = "Set1"
    } else {
      pal = "Dark2"
    }
    ggplot(data = df, aes(x = `UMAP 1`, y = `UMAP 2`)) +
      geom_point(
        aes(
          color = factor(Clustering),
        ),
        size = 0.5
      ) +
      facet_wrap(
        .~Subset, 
        scales = "free"
      ) +
      labs(
        color = "",
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
      guides(color = guide_legend(override.aes = list(size=2))) +
      scale_color_brewer(palette = pal) +  
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
      paste0(
        "outs/umap/figures/",
        save_prefix,
        ".pdf"
      ),
      width = 16,
      height = 8,
      device = cairo_pdf
    )
  }
}

# Iterate over the umap files and names and save the results 
mapply(
  umap_plot,
  df = umap_loaded, 
  save_prefix = umap_names
)
