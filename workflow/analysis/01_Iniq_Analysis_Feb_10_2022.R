library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)

# Change to top level dir
setwd("../../results/")

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)

# Merge imbalance and clustering summary results
colnames(clus_concat)[2] <- "Number of batches downsampled"
clus_imba_merged <- merge(
  imba_concat,
  clus_concat,
  by = c(
    "Dataset",
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Replicate"
  ),
  allow.cartesian = TRUE
)

# Plot initial results for celltype ARI/AMI/Homogeneity
ggplot(data = clus_imba_merged, aes(
  x = `Celltype intersection ratio`, 
  y = `Celltype ARI`
)) + 
  theme_few() +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  labs(
    fill = "Density",
    x = "Celltype intersection ratio (Jaccard distance)"
  )

ggplot(data = clus_imba_merged, aes(
  x = `Mean proportion cosine distance`, 
  y = `Celltype ARI`
)) + 
  theme_few() +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  labs(
    fill = "Density",
    x = "Mean proportion cosine distance"
  )

ggplot(data = clus_imba_merged, aes(
  x = `Length coeff var`, 
  y = `Celltype ARI`
)) + 
  theme_few() +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  labs(
    fill = "Density",
    x = "Batch size coefficient of variation"
  )

# Plot Celltype ARI as barplot based on method 
clus_imba_merged_hindbrain <- clus_imba_merged[
  clus_imba_merged$Dataset %in% c("pbmc_4_batch")
]

ggplot(data = clus_imba_merged_hindbrain, aes(
  x = factor(`Proportion downsampled.x`), 
  y = `Celltype ARI`
)) + 
  theme_few() +
  geom_boxplot(aes(fill = factor(`Number of celltypes downsampled`))) 
