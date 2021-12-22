library(Seurat)
library(SeuratDisk)
library(reticulate)

# Read in matrix for full data, including last batch column 
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
filename <- args[2]
int_type <- args[3] # Integration type

# Convert h5ad file to h5seurat 
Convert(file, dest = "h5seurat", overwrite = TRUE)

# Load anndata and scanpy 
ad <- import("anndata")
sc <- import("scanpy")

# Create Seurat object and split by batch information
seur_obj <- LoadH5Seurat(paste0("./tmp/", filename, ".h5seurat"))
seur_obj_list <- SplitObject(
    seur_obj,
    split.by = "batch"
)

# Iterate over batches and find highly variable genes 
for (i in 1:length(seur_obj_list)) {
    seur_obj_list[[i]] <- FindVariableFeatures(
        seur_obj_list[[i]], 
        selection.method = "mean.var.plot",
        nfeatures = 2500, 
        verbose = TRUE
    )
}

# Determine type of integration to perform (CCA or RPCA)
if (int_type == "CCA") {
    int_anchors <- FindIntegrationAnchors(
        object.list = seur_obj_list,
        dims = 1:20,
        anchor.features = 2500
    )
    batches_integrated <- IntegrateData(anchorset = int_anchors, dims = 1:20)
} else if (int_type == "RPCA") {
    int_features <- SelectIntegrationFeatures(object.list = seur_obj_list)
    for (i in 1:length(seur_obj_list)) {
        x <- seur_obj_list[[i]]
        x <- ScaleData(x, features = int_features, verbose = TRUE)
        x <- RunPCA(x, features = int_features, verbose = TRUE)
        seur_obj_list[[i]] <- x
    }
    int_anchors <- FindIntegrationAnchors(
        object.list = seur_obj_list,
        reduction = "rpca",
        dims = 1:20
    )
    batches_integrated <- IntegrateData(anchorset = int_anchors, dims = 1:20)
} else {
    stop(
        "Please indicate either 'CCA' or 'RPCA' for the integration type option"
    )
}

# Return integrated adata object as hda5 file -> tempfile
SaveH5Seurat(
    object = batches_integrated, 
    filename = paste0("./tmp/", filename, "_seur_out.h5Seurat"),
    overwrite = TRUE,
    verbose = TRUE
)

# Convert tempfile to h5ad object
Convert(
    paste0("./tmp/", filename, "_seur_out.h5Seurat"),
    dest = "h5ad"
)
