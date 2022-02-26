import numpy as np
import scanpy as sc
import anndata as ann

def downsample(adata, num_celltypes = None, celltype_names = None, proportion = 0.5):
    # Initialize random number generator
    rng = np.random.default_rng()
    
    # For the given number of celltypes, select num_celltypes
    # randomly, unless non random indicated by celltype_names
    if celltype_names is not None:
        celltypes_sample = celltype_names
    else:
        if num_celltypes is None:
            raise ValueError(
                "num_celltypes and celltype_names cannot both be None"
            )
        if num_celltypes == 0:
            celltypes_sample = "None"
            return adata, celltypes_sample
        unique_celltypes = np.unique(adata.obs["celltype"].__array__())
        rng.shuffle(unique_celltypes)
        celltypes_sample = rng.choice(unique_celltypes, num_celltypes, replace = False)
        
    # Save the original batch label for later 
    adata.obs["batch_orig"] = adata.obs["batch"]
    
    # Downsample selected celltypes by given proportion
    for celltype in celltypes_sample:
        adata_celltype = adata[adata.obs["celltype"] == celltype]
        adata_noncelltype = adata[adata.obs["celltype"] != celltype]
        if proportion == 0:
            adata = adata_noncelltype
            continue
        adata_celltype_ds = sc.pp.subsample(
            adata_celltype, 
            fraction = proportion,
            random_state = None,
            copy = True
        )
        adata = ann.AnnData.concatenate(adata_noncelltype, adata_celltype_ds)

    # Replace batch column with batch original and drop batch_orig
    adata.obs["batch"] = adata.obs["batch_orig"]
    adata.obs.drop("batch_orig", axis = 1, inplace = True)
        
    # Return downsampled data + sampled celltypes   
    return adata, celltypes_sample