library(data.table)
library(liger)
library(Seurat)
library(SeuratDisk)
library(reticulate)

# Read in matrix for full data, including last batch column 
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
filename <- args[2]

# Load anndata, scanpy, and scipy sparse through reticulate
ad <- import("anndata")
sc <- import("scanpy")
sp_sparse <- import("scipy.sparse")

# Convert h5ad file to h5seurat 
Convert(file, dest = "h5seurat", overwrite = TRUE)

# Create Seurat object and split by batch information
seur_obj <- LoadH5Seurat(paste0("./tmp/", filename, ".h5seurat"))
seur_obj_list <- SplitObject(
    seur_obj,
    split.by = "batch"
)

# Get matrices of rna data for each batch and name by batch 
seur_counts_list <- lapply(seur_obj_list, function(x) {
    return(x@assays$RNA@counts)
})
seur_counts_list_names <- lapply(seur_obj_list, function(x) {
    return(unique(x@meta.data$batch))
})
names(seur_counts_list) <- seur_counts_list_names

# Create Liger object from seurat list of matrices
liger_obj <- createLiger(seur_counts_list, remove.missing = FALSE)

# Add norm information and HVG information (same as raw 
# object information)
liger_obj@norm.data <- liger_obj@raw.data
liger_obj@var.genes <- rownames(liger_obj@raw.data[[1]])

# Scale data, perform iNFM and quantile normalization
liger_obj <- scaleNotCenter(liger_obj)
liger_obj <- optimizeALS(liger_obj, k = 20)
liger_obj <- quantile_norm(liger_obj)

# Extract normalized cell loadings, save as h5seurat object,
# and convert to h5ad
liger_norm_h <- liger_obj@H.norm
rownames(liger_norm_h) <- colnames(seur_obj)
colnames(liger_norm_h) <- paste0(
    "h_norm_comp_", seq(1:ncol(liger_norm_h))
)
norm_cell_loadings <- CreateSeuratObject(counts = liger_norm_h)
SaveH5Seurat(
    object = norm_cell_loadings, 
    filename = paste0("./tmp/", filename, "_liger_out.h5seurat"),
    overwrite = TRUE,
    verbose = TRUE
)

# Convert tempfile to h5ad object
Convert(
    paste0("./tmp/", filename, "_liger_out.h5seurat"),
    dest = "h5ad"
)