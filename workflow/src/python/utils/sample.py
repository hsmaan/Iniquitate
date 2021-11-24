import numpy as np
import anndata as ann

def downsample(adata, num_celltypes = None, celltype_names = None, proportion = 0.5):
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
            return adata
        unique_celltypes = np.unique(adata.obs["celltype"].__array__())
        celltypes_sample = np.random.choice(unique_celltypes, num_celltypes)
    
    # Downsample selected celltypes by given proportion
    for celltype in celltypes_sample:
        adata_celltype = adata[adata.obs["celltype"] == celltype]
        adata_noncelltype = adata[adata.obs["celltype"] != celltype]
        adata_celltype_indices_ds = np.random.choice(
            [i for i in range(len(adata_celltype))],
            int(round(len(adata_celltype)*proportion, 0))
        )
        adata_celltype_ds = adata_celltype[adata_celltype_indices_ds]
        adata = ann.AnnData.concatenate(adata_noncelltype, adata_celltype_ds)
        
    # Return downsampled data
    return adata