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

# Change to results dir for lowcap modified data 
setwd("../../../results/lowcap_modified/")

##### Analysis of PBMC 2 batch imbalanced data #####

# Load in and concatenate celltype imbalance summary files
setwd("celltype_imbalance_summaries")
cimba_files <- list.files()
cimba_files <- grep(
  "pbmc_2_batch",
  cimba_files,
  value = TRUE
)
cimba_loaded <- lapply(cimba_files, fread)
cimba_concat <- Reduce(rbind, cimba_loaded)
gc()

# Load in and concatenate full imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch",
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
  "pbmc_2_batch",
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
  "pbmc_2_batch",
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
  "pbmc_2_batch",
  dge_files,
  value = TRUE
)
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)
gc()


# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_files <- grep(
  "pbmc_2_batch",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
gc()

# Change to top level dir 
setwd("../../..")

