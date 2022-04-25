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

##### Analysis of PBMC 2 batch imbalanced data vs balanced data ##### - 
#### with respect to ARI and other integration metrics 

# Change to results dir for lowcap modified data 
setwd("../../../results/lowcap_modified/")

# Load in and concatenate celltype imbalance summary files
setwd("celltype_imbalance_summaries")
cimba_files <- list.files()
pbmc_2_batch_imba_cimba_files <- grep(
  "pbmc_2_batch",
  cimba_files,
  value = TRUE
)
pbmc_2_batch_imba_cimba_loaded <- lapply(
  pbmc_2_batch_imba_cimba_files, 
  fread
)
pbmc_2_batch_imba_cimba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_imba_cimba_loaded
)
gc()

# Load in and concatenate full imbalance summary files 
setwd("../imbalance_summaries/")
imba_files <- list.files()
pbmc_2_batch_imba_imba_files <- grep(
  "pbmc_2_batch",
  imba_files,
  value = TRUE
)
pbmc_2_batch_imba_imba_loaded <- lapply(
  pbmc_2_batch_imba_imba_files, 
  fread
)
pbmc_2_batch_imba_imba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_imba_imba_loaded
)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
pbmc_2_batch_imba_clus_files <- grep(
  "pbmc_2_batch",
  clus_files,
  value = TRUE
)
pbmc_2_batch_imba_clus_loaded <- lapply(
  pbmc_2_batch_imba_clus_files, 
  fread
)
pbmc_2_batch_imba_clus_concat <- Reduce(
  rbind, 
  pbmc_2_batch_imba_clus_loaded
)
gc()

# Load in relatedness metric summary 
setwd("../relatedness_results/")
relate_files <- list.files()
pbmc_2_batch_imba_relate_file <- grep(
  "pbmc_2_batch",
  relate_files,
  value = TRUE
)
pbmc_2_batch_imba_relate_loaded <- fread(pbmc_2_batch_imba_relate_file)

# Change to results dir for control data 
setwd("../../control/")

# Load in and concatenate celltype imbalance summary files
setwd("celltype_imbalance_summaries")
cimba_files <- list.files()
pbmc_2_batch_bal_cimba_files <- grep(
  "pbmc_2_batch_base_balanced",
  cimba_files,
  value = TRUE
)
pbmc_2_batch_bal_cimba_loaded <- lapply(
  pbmc_2_batch_bal_cimba_files, 
  fread
)
pbmc_2_batch_bal_cimba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_bal_cimba_loaded
)
gc()

# Load in and concatenate full imbalance summary files 
setwd("../imbalance_summaries/")
imba_files <- list.files()
pbmc_2_batch_bal_imba_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_files,
  value = TRUE
)
pbmc_2_batch_bal_imba_loaded <- lapply(
  pbmc_2_batch_bal_imba_files, 
  fread
)
pbmc_2_batch_bal_imba_concat <- Reduce(
  rbind, 
  pbmc_2_batch_bal_imba_loaded
)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
pbmc_2_batch_bal_clus_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_files,
  value = TRUE
)
pbmc_2_batch_bal_clus_loaded <- lapply(
  pbmc_2_batch_bal_clus_files, 
  fread
)
pbmc_2_batch_bal_clus_concat <- Reduce(
  rbind, 
  pbmc_2_batch_bal_clus_loaded
)
gc()

# Load in relatedness metric summary 
setwd("../relatedness_results/")
relate_files <- list.files()
pbmc_2_batch_bal_relate_file <- grep(
  "pbmc_2_batch_base_balanced",
  relate_files,
  value = TRUE
)
pbmc_2_batch_bal_relate_loaded <- fread(pbmc_2_batch_bal_relate_file)

# Change to top level dir 
setwd("../../..")

# Make lowcap figure directory if it does not exist
if (!dir.exists("outs/lowcap_modified/figures")) {
  dir.create("outs/lowcap_modified/figures", recursive = TRUE)
}

### Fig 4A) - Analysis of PBMC 2 batch results with respect to integration
### and measures of relatedness 

# Format the relatedness metric results for both the balanced and imbalanced
# PBMC 2 batch datasets 
pbmc_2_batch_bal_relate_formatted <- pbmc_2_batch_bal_relate_loaded %>%
  group_by(`Celltype 1`, `Celltype 2`) %>%
  summarize(`Average PCA cosine dist` = mean(`PCA cosine dist`)) %>%
  as.data.frame() 

pbmc_2_batch_imba_relate_formatted <- pbmc_2_batch_imba_relate_loaded %>%
  group_by(`Celltype 1`, `Celltype 2`) %>%
  summarize(`Average PCA cosine dist` = mean(`PCA cosine dist`)) %>%
  as.data.frame() 

# Merge and subset balanced clustering data to only include control experiments 
imba_clus_merged_bal <- merge(
  pbmc_2_batch_bal_clus_concat,
  pbmc_2_batch_bal_imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)
imba_clus_merged_bal$type <- ifelse(
  imba_clus_merged_bal$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged_bal$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)
imba_clus_merged_bal_control <- imba_clus_merged_bal[
  imba_clus_merged_bal$type %in% c("Control")
]

# Merge imbalanced clustering data
imba_clus_merged_imba <- merge(
  pbmc_2_batch_imba_clus_concat,
  pbmc_2_batch_imba_imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Create dataframe of ARI celltype and batch values for both datasets 
imba_clus_imba_sub <- imba_clus_merged_imba[,
  c("Method", "Celltype ARI Imbalanced", "Batch ARI")
]
imba_clus_imba_sub$Dataset <- "PBMC full \ndataset imbalanced"

imba_clus_bal_sub <- imba_clus_merged_bal_control[,
  c("Method", "Celltype ARI Imbalanced", "Batch ARI")
]
imba_clus_bal_sub$Dataset <- "PBMC control \ndataset balanced"

imba_clus_imba_bal_sub_merged <- rbind(
  imba_clus_imba_sub, imba_clus_bal_sub
)

# Plot comparison of ARI - celltype and batch - values between the two datasets
ggplot(
  data = imba_clus_imba_bal_sub_merged, 
  aes(
    x = Method, 
    y = `Celltype ARI Imbalanced`
  )
) +
  geom_point(aes(color = Dataset), position = position_jitterdodge()) +
  geom_boxplot(
    aes(fill = Dataset), 
    position = "dodge2",
    alpha = 0.8,
    notch = TRUE
  ) +
  ylim(c(0, 1)) +
  labs(
    x = "Method",
    y = "Celltype ARI"
  ) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/lowcap_modified/figures/",
    "10_pbmc_full_dataset_vs_control_dataset_celltype_ari_no_ds.pdf"
  ),
  width = 9,
  height = 6,
  device = cairo_pdf
)  

ggplot(
  data = imba_clus_imba_bal_sub_merged, 
  aes(
    x = Method, 
    y = `Batch ARI`
  )
) +
  geom_point(aes(color = Dataset), position = position_jitterdodge()) +
  geom_boxplot(
    aes(fill = Dataset), 
    position = "dodge2",
    alpha = 0.8,
    notch = TRUE
  ) +
  labs(
    x = "Method",
    y = "(1 - batch) ARI"
  ) +
  ylim(c(0, max(imba_clus_imba_bal_sub_merged$`Batch ARI`))) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  theme_few() +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  theme(axis.text.y = element_text(size = 14)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))
ggsave(
  paste0(
    "outs/lowcap_modified/figures/",
    "10_pbmc_full_dataset_vs_control_dataset_batch_ari_no_ds.pdf"
  ),
  width = 9,
  height = 6,
  device = cairo_pdf
)

##### Plotting of relatedness metric for both PBMC 2 batch balanced and ##### 
##### imbalanced datasets #####
ggplot(
  data = pbmc_2_batch_bal_relate_formatted,
  aes(
    x = `Celltype 1`,
    y = `Celltype 2`
  ) 
) + 
  geom_tile(
    aes(fill = `Average PCA cosine dist`)
  ) +
  scale_fill_gradient(
    low = "firebrick2",
    high = "white"
  ) +
  scale_x_discrete(
    limits = rev(levels(factor(pbmc_2_batch_bal_relate_formatted$`Celltype 1`)))
  ) + 
  scale_y_discrete(
    limits = rev(levels(factor(pbmc_2_batch_bal_relate_formatted$`Celltype 2`)))
  ) + 
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
  theme(aspect.ratio = 1) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  labs(
    fill = "Average PCA distance \n across batches",
    x = "Celltype 1",
    y = "Celltype 2"
  )
ggsave(
  "outs/lowcap_modified/figures/10_pbmc_2_batch_bal_celltype_relatedness.pdf",
  width = 7,
  height = 7,
  device = cairo_pdf
)

ggplot(
  data = pbmc_2_batch_imba_relate_formatted,
  aes(
    x = `Celltype 1`,
    y = `Celltype 2`
  ) 
) + 
  geom_tile(
    aes(fill = `Average PCA cosine dist`)
  ) +
  scale_fill_gradient(
    low = "firebrick2",
    high = "white"
  ) +
  scale_x_discrete(
    limits = rev(levels(factor(pbmc_2_batch_imba_relate_formatted$`Celltype 1`)))
  ) + 
  scale_y_discrete(
    limits = rev(levels(factor(pbmc_2_batch_imba_relate_formatted$`Celltype 2`)))
  ) + 
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
  theme(aspect.ratio = 1) + 
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 14)) +
  labs(
    fill = "Average PCA distance \n across batches",
    x = "Celltype 1",
    y = "Celltype 2"
  )
ggsave(
  "outs/lowcap_modified/figures/10_pbmc_2_batch_imba_celltype_relatedness.pdf",
  width = 7,
  height = 7,
  device = cairo_pdf
)

##### Analysis of relatedness metric and cell number with respect to #####
##### KNN clasification results for each dataset analyzed #####

# Define datasets and keywords to grep
datasets <- c(
  "pbmc_2_batch",
  "pbmc_4_batch",
  "mouse_hindbrain_6_batch",
  "peng_pdac_8_batch"
)

# Load relatedness metrics for all four datasets
setwd("results/lowcap_modified/relatedness_results/")
relatedness_files <- list.files()
relatedness_loaded <- lapply(
  datasets, function(x) {
    file_grep <- grep(x, relatedness_files, value = TRUE)
    file_loaded <- fread(file_grep)
    return(file_loaded)
  }
)
names(relatedness_loaded) <- datasets

# Load in KNN classification results for all four datasets
setwd("../knn_classification_reports/")
knn_class_files <- list.files()
knn_class_files_loaded <- lapply(
  datasets, function(x) {
    files_grep <- grep(x, knn_class_files, value = TRUE)
    files_loaded <- lapply(
      files_grep, fread
    )
    files_concat <- Reduce(rbind, files_loaded)
    return(files_concat)
  }
)
names(knn_class_files_loaded) <- datasets

# Load in celltype imbalance results for all four datasets 
setwd("../celltype_imbalance_summaries/")
celltype_imba_files <- list.files()
celltype_imba_files_loaded <- lapply(
  datasets, function(x) {
    files_grep <- grep(x, celltype_imba_files, value = TRUE)
    files_loaded <- lapply(
      files_grep, fread
    )
    files_concat <- Reduce(rbind, files_loaded)
    return(files_concat)
  }
)
names(celltype_imba_files_loaded) <- datasets

# Change to top level dir
setwd("../../..")

# Create function to plot concordance of KNN classification performance
# with respect to relatedness metrics 
knn_relatedness_plot <- function(
    dataset, relatedness_df, knn_class_df, plot_height, plot_width
  ) {
  # Format relatedness df to take the minimum (distance) of each celltype 
  # to another celltype in the dataset, while excluding self-matches 
  relatedness_df_sub <- relatedness_df[
    relatedness_df$`Celltype 1` != relatedness_df$`Celltype 2`
  ]
  relatedness_df_min <- relatedness_df_sub %>%
    group_by(`Celltype 1`) %>% 
    summarize(var = min(`PCA cosine dist`)) %>%
    as.data.frame()
  colnames(relatedness_df_min) <- c("Celltype", "Min PCA cosine dist")
  
  # Merge together relatedness and knn class df
  relatedness_knn_class_merged <- merge(
    relatedness_df_min,
    knn_class_df,
    by = c(
      "Celltype"
    )
  )
  
  # Plot and save the relationship between the min celltype relatedness and 
  # f1-scores
  ggplot(
    data = relatedness_knn_class_merged, 
    aes(
      x = `Celltype`,
      y = `F1-score`
    )
  ) +
    scale_x_discrete(limits = rev(
      levels(
        factor(relatedness_knn_class_merged$Celltype)
        )
      )
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = (
        c(
          min(relatedness_knn_class_merged$`F1-score`), 
          max(relatedness_knn_class_merged$`F1-score`)
        )
      ),
      oob = scales::squish
    ) +
    facet_wrap(.~Method, scales = "free_x") +
    geom_jitter(aes(color = `Min PCA cosine dist`)) +
    geom_boxplot(aes(fill = `Min PCA cosine dist`)) +
    scale_fill_gradient(low = "blue", high = "red", na.value = NA) +
    scale_color_gradient(low = "blue", high = "red", na.value = NA) +
    labs(
      fill = "Minimum \npairwise \ncelltype \ndistance",
      x = "Celltype",
      y = "F1-classification score post-integration"
    ) +
    guides(
      color = "none"
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
    paste0(
      "outs/lowcap_modified/figures/10_",
      dataset,
      "_knn_results_celltype_relatedness_comp.pdf"
    ),
    width = plot_width,
    height = plot_height,
    device = cairo_pdf
  )

}

# Create lists of heights and widths to save plots
heights <- list(
  7, 7, 20, 7
)
widths <- list(
  16, 16, 16, 16
)

# Iterate over datasets, relatedness_dfs, knn_class_dfs and plot
mapply(
  knn_relatedness_plot,
  dataset = datasets,
  relatedness_df = relatedness_loaded,
  knn_class_df = knn_class_files_loaded,
  plot_height = heights,
  plot_width = widths
)

# Create function to plot concordance of KNN classification performance
# with respect to support for the given celltypes

knn_support_plot <- function(dataset, knn_class_df, plot_height, plot_width) {
  # Add log of support to knn class df
  knn_class_df$Log_support <- log2(knn_class_df$Support)
  
  # Plot the relationship between celltype support and f1-scores and save
  ggplot(
    data = knn_class_df, 
    aes(
      x = `Celltype`,
      y = `F1-score`
    )
  ) +
    scale_x_discrete(
      limits = rev(levels(factor(knn_class_df$Celltype)))
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = (
        c(
          min(knn_class_df$`F1-score`), 
          max(knn_class_df$`F1-score`)
        )
      ),
      oob = scales::squish
    ) +
    facet_wrap(.~Method, scales = "free_x") +
    scale_fill_gradient(low = "goldenrod1", high = "firebrick2", na.value = NA) +
    scale_color_gradient(low = "goldenrod1", high = "firebrick2", na.value = NA) +
    geom_jitter(aes(color = Log_support)) +
    geom_boxplot(aes(fill = Log_support)) +
    labs(
      fill = "Log2 \ncelltype \nnumber",
      x = "Celltype",
      y = "F1-classification score post-integration"
    ) +
    guides(
      color = "none"
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
    paste0(
      "outs/lowcap_modified/figures/10_",
      dataset,
      "_knn_results_celltype_num_comp.pdf"
    ),
    width = plot_width,
    height = plot_height,
    device = cairo_pdf
  )
  
}

# Create lists of heights and widths to save plots
heights <- list(
  7, 7, 20, 7
)
widths <- list(
  16, 16, 16, 16
)

# Iterate over datasets, knn_class_dfs and plot and save 
mapply(
  knn_support_plot,
  dataset = datasets,
  knn_class_df = knn_class_files_loaded, 
  plot_height = heights,
  plot_width = widths
)
  
