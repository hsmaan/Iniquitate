library(Seurat)
library(SeuratDisk)
library(reticulate)

# Read in arguments  
args <- commandArgs(trailingOnly = TRUE)
ref_file <- args[1]
temp_adata_file <- args[2]
temp_adata_filename <- args[3]
outfile_name <- args[4] 

# Read in h5seurat reference data 
ref_data <- LoadH5Seurat(ref_file) 

# Load anndata and scanpy 
ad <- import("anndata")
sc <- import("scanpy")

# Convert h5ad anndata temp file 
temp_adata <- ad$read_h5ad(temp_adata_file)

# Create Seurat object and split by batch information - use anndata import
exprs <- t(temp_adata$X$todense())
colnames(exprs) <- temp_adata$obs_names$to_list()
rownames(exprs) <- temp_adata$var_names$to_list()
query_obj <- CreateSeuratObject(exprs)
query_obj <- SetAssayData(query_obj, "data", exprs)
query_obj <- AddMetaData(query_obj, temp_adata$obs)
query_obj_list <- SplitObject(
    query_obj,
    split.by = "batch"
) 

# Normalize query batches using scTransform  
query_obj_list <- lapply(X = query_obj_list, FUN = SCTransform, verbose = FALSE)

# Get anchors between each query batches and the reference 
anchors <- list()
for (i in 1:length(query_obj_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = ref_data,
    query = query_obj_list[[i]],
    reference.reduction = "spca", 
    dims = 1:50
  )
}

# Note - this may not be ideal to simulate effects of downsampling 
# as each batch is being mapped individually here and not separately
# Map each of the query batches individually 
for (i in 1:length(query_obj_list)) {
  query_obj_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query_obj_list[[i]],
    reference = ref_data, 
    refdata = list(
        celltype.l1 = "celltype.l1",
        celltype.l2 = "celltype.l2",
        predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}

# Set default assay of all to SCT (for outputting scTransformed counts)
for (i in 1:length(query_obj_list)) {
  DefaultAssay(query_obj_list[[i]]) <- "SCT"
}

# Remerge the batches into one object - reset default assay as failsafe 
query_ref_mapped_obj <- Reduce(merge, query_obj_list)
DefaultAssay(query_ref_mapped_obj) <- "SCT"

# Return reference mapped adata object as hda5 file -> tempfile
SaveH5Seurat(
    object = query_ref_mapped_obj, 
    filename = paste0(outfile_name, ".h5Seurat"),
    overwrite = TRUE,
    verbose = TRUE
)

# Convert tempfile to h5ad object
Convert(
    paste0(outfile_name, ".h5Seurat"),
    dest = "h5ad"
)

# Remove h5seurat file 
file.remove(
  paste0(outfile_name, ".h5Seurat")
)