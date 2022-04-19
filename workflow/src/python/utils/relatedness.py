import numpy as np
import pandas as pd
import scipy.spatial as sp
import scanpy as sc 

def relatedness_score(adata, pca_performed = True):
    """Computes the relatedeness between celltypes using cosine distance in PCA space
    
    Args:
        adata (AnnData): AnnData object containing relevant count information with celltype
            and batch observations.
        pca_performed (bool): True or False value indicating whether PCA decomposition steps
            have been performed already for AnnData object. Default is True.
    """
    # Perform PCA if not already performed
    if pca_performed is False:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2500)
        sc.pp.pca(adata)
        
    # Get the batch and celltype information from AnnData object
    batch_vals = np.unique(adata.obs["batch"].__array__())
    if len(batch_vals) > 1:
        raise ValueError("More than one batch found in AnnData object")
    batch = batch_vals[0]
    celltypes = np.unique(adata.obs["celltype"].__array__())
    
    # Utilize the cosine distance between the average PCA embedding for celltype i and 
    # average PCA embedding for celltype j
    pca_top_20 = adata.obsm["X_pca"][:, 0:20]
    top_20_pc_weights = adata.uns["pca"]["variance_ratio"][0:20]
    celltype_is = []
    celltype_js = []
    pca_cosine_dists = []
    for celltype_i in celltypes:
        for celltype_j in celltypes:
            celltype_is.append(celltype_i)
            celltype_js.append(celltype_j)
            pca_celltype_i = pca_top_20[
                adata.obs.celltype == celltype_i
            ]
            pca_celltype_j = pca_top_20[
                adata.obs.celltype == celltype_j
            ]
            pca_celltype_i_avg = np.sum(pca_celltype_i, axis = 0)/len(pca_celltype_i)
            pca_celltype_j_avg = np.sum(pca_celltype_j, axis = 0)/len(pca_celltype_j)
            pca_cosine_dist = sp.distance.cosine(
                pca_celltype_i_avg,
                pca_celltype_j_avg,
                w = top_20_pc_weights
            )
            pca_cosine_dists.append(pca_cosine_dist)
            
    # Gather the cosine distance results in a dataframe and return  
    dist_results_df = pd.DataFrame({
        "Celltype 1": celltype_is,
        "Celltype 2": celltype_js,
        "PCA cosine dist": pca_cosine_dists,
        "Batch": batch
    })
    return dist_results_df