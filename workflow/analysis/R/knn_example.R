# Load the necessary libraries - please note that the analysis environment
# (found in envs/analysis.yaml) should be used to run this script
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

# Change to results dir for the custom data - 
# change this directory as necessary
setwd("../../../results/custom")

# Load in and concatenate imbalance summary files
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)

# Change to top level dir
setwd("../../..")

# Create directory for output of results and figures
if (!dir.exists("outs/custom/figures")) {
  dir.create("outs/custom/figures", recursive = TRUE)
}
if (!dir.exists("outs/custom/results")) {
  dir.create("outs/custom/results", recursive = TRUE)
}

### Results of celltype downsampling and ablation on  
### KNN classification scores 

# Merge imbalance and knn classification results together - several 
# other analysis can be done with this data now, but we'll only 
# highlight the main KNN analysis that was done in the paper 
imba_knn_merged <- merge(
  imba_concat,
  knn_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_knn_merged <- distinct(imba_knn_merged)

# Subset for only cases where the celltype downsampled is equal to the 
# celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype$type <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Create function to format facet labels (downsampled celltypes) and plot 
# the results (this figure will be in the same format as the paper)
ds_celltype_labelled <- function(variable,value){
  return(paste0("Cell-type affected = ", value))
}

ggplot(data = imba_knn_merged_celltype, aes(x = `Method`, y = `F1-score`)) +
  geom_boxplot(
    aes(
      fill = factor(`type`, levels = c("Control", "Downsampled", "Ablated")),
    ),
    notch = FALSE,
    alpha = 0.8 
  ) +
  ylim(0, 1) +
  facet_wrap(
    .~Celltype, 
    scales = "free_x", 
    labeller = ds_celltype_labelled,
    ncol = 2
  ) +
  labs(
    fill = "Type",
    x = "Method",
    y = "Affected celltype F1-classification score post-integration"
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
  theme(legend.text = element_text(size = 16))
# Change the name and extension of file as necessary
ggsave(
  "outs/custom/figures/ds_ablate_allmethod_knn_f1_score.pdf",
  width = 12,
  height = 14,
  device = cairo_pdf
)