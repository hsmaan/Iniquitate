library(data.table)
library(liger)
library(reticulate)

# Read in matrix for full data, including last batch column 
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
filename <- args[2]
conda_env <- args[3] # Conda env location

# Set conda env 
use_condaenv(condaenv = conda_env)

# Load anndata, scanpy, and scipy sparse through reticulate
ad <- import("anndata")
sc <- import("scanpy")
sp_sparse <- import("scipy.sparse")

# Load matrix, format, and subset batch data
full_mat <- read.csv(
    file,
    sep = "\t"
)
rownames(full_mat) <- full_mat[, 1]
full_mat <- full_mat[, -1]
batch_names <- full_mat[, ncol(full_mat)]
full_mat <- full_mat[, -ncol(full_mat)]
names(batch_names) <- rownames(full_mat)

# Split matrix into list of matrices by batch information 
mat_split <- split(full_mat, batch_names)

# Transpose and create LIGER object
mat_split_t <- lapply(mat_split, t)
liger_obj <- createLiger(mat_split_t, remove.missing = FALSE)

# Add norm information and HVG information (same as raw 
# object information)
liger_obj@norm.data <- liger_obj@raw.data
liger_obj@var.genes <- rownames(liger_obj@raw.data[[1]])

# Scale data, perform iNFM and quantile normalization
liger_obj <- scaleNotCenter(liger_obj)
liger_obj <- optimizeALS(liger_obj, k = 20)
liger_obj <- quantile_norm(liger_obj)

# Extra normalized cell loadings, convert to datatable and save
norm_loadings <- as.data.table(liger_obj@H.norm)
fwrite(
    norm_loadings,
    file = paste0(
        "./tmp/",
        filename,
        "_liger_out.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)
