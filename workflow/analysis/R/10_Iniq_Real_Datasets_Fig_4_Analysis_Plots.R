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

### Fig 4A) - Analysis of PBMC 2 batch results with respect to integration
### and measures of relatedness 

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



