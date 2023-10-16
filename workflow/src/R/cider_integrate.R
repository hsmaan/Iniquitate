library(CIDER)
library(Seurat)
library(SeuratDisk)
library(reticulate)
library(parallel)
library(doParallel)
library(plyr)

# Read in matrix for full data, including last batch column 
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
filename <- args[2]

# Load anndata and scanpy 
ad <- import("anndata")
sc <- import("scanpy")

# Load h5ad object through reticulate and create Seurat object
temp_adata <- ad$read_h5ad(file) 
exprs <- t(temp_adata$X$todense())
colnames(exprs) <- temp_adata$obs_names$to_list()
rownames(exprs) <- temp_adata$var_names$to_list()
seur_obj <- CreateSeuratObject(exprs)
seur_obj <- SetAssayData(seur_obj, "data", exprs)
seur_obj <- AddMetaData(seur_obj, temp_adata$obs)

# Format the batch names 
batch_names <- seur_obj$batch
batch_names_unique <- unique(batch_names)
batch_names_encoded <- as.numeric(factor(batch_names)) 
seur_obj$Batch <- batch_names_encoded

# Use the high-level CIDER pipeline
seur_obj <- initialClustering(seur_obj, nfeatures = 2500, resolution = 1.0, dims = 1:20)
ider <- getIDEr(seur_obj, downsampling.size = 35, use.parallel = TRUE, n.cores = 8, verbose = FALSE)
seur_obj <- finalClustering(seur_obj, ider, cutree.h = 0.35) # final clustering

# Extract the cluster labels and append to a dataframe 
cluster_labels <- as.character(seur_obj$CIDER_cluster)
batch_labels <- as.character(seur_obj$batch)
cluster_df <- data.frame(
    "batch" = batch_labels,
    "cider_cluster" = cluster_labels
)

# Write out cluster labels to file
write.table(
    cluster_df,
    paste0("./tmp/", filename, "_cider_out.txt"),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
)                
