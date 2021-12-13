import numpy as np
import pandas as pd
from sklearn import metrics

def cluster_concordance(adata):
    # Pull out clustering values per subset
    methods = ["harmony", "scvi", "bbknn", "scanorama", "seurat", "liger"]
    adata_subsets = []
    for method in methods:
        adata_subsets.append(adata[adata.obs["integration_method"] == method])
    cluster_subsets = [
        adata_subset.obs["leiden"].__array__() for adata_subset in adata_subsets 
    ]
        
    # Get ARI values across subsets 
    ari_vals_mat = np.zeros((len(cluster_subsets), len(cluster_subsets)))
    for i, clusters_i in enumerate(cluster_subsets):
        for j, clusters_j in enumerate(cluster_subsets):
            ari_val = metrics.adjusted_rand_score(clusters_i, clusters_j)
            ari_vals_mat[i, j] = ari_val
        
    # Create dataframe of values
    ari_val_df = pd.DataFrame(ari_vals_mat)
    ari_val_df.index = methods
    ari_val_df.columns = methods
    
    # Get global ARI value (median)
    ari_vals_mat_nan_diag = np.zeros((len(cluster_subsets), len(cluster_subsets)))
    for i, clusters_i in enumerate(cluster_subsets):
        for j, clusters_j in enumerate(cluster_subsets):
            if i == j:
                ari_vals_mat_nan_diag[i, j] = np.nan
            else:
                ari_val = metrics.adjusted_rand_score(clusters_i, clusters_j)
                ari_vals_mat_nan_diag[i, j] = ari_val
    
    ari_vals_mat_no_diag = ari_vals_mat_nan_diag[~np.isnan(ari_vals_mat_nan_diag)]
    median_ari = np.median(ari_vals_mat_no_diag)
    
    # Convert concordance dataframe to long format 
    ari_val_df_long = ari_val_df.melt(ignore_index = False)
    ari_val_df_long = ari_val_df_long.reset_index()
    ari_val_df_long.columns = ["Method 1", "Method 2", "ARI"]

    # Append median ARI value to dataframe and return
    ari_val_df_long["Median ARI"] = median_ari
    
    return ari_val_df_long