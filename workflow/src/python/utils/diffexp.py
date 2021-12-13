import numpy as np
import pandas as pd 
import scanpy as sc 
import anndata as ann

def diffexp(adata, groupby, **kwargs):
    """
    Perform differential expression analysis on an AnnData object.

    Args:
        adata (AnnData):Annotated data matrix object.
        groupby (str): The column name of the dataframe to group by.
        **kwargs: Keyword arguments to be passed to scanpy.tl.rank_genes_groups.

    Returns:
        adata (AnnData): Annotated data matrix object with differential expression analysis results.
    """
    sc.tl.rank_genes_groups(adata, groupby = groupby, **kwargs)
    return adata

def dge_top_n(adata, n, obs_group):
    """
    Return clusters and genes with the top n differential expression.

    Args:
        adata (AnnData): Annotated data matrix object.
        n (int): The number of top differentially expressed genes to return.
        obs_group (str): The column name in obs of adata object to group by.
    Returns:
        data (DataFrame): A dataframe with the top n differentially expressed genes in each cluster.
    """
    unique_groups = np.sort(np.unique(adata.obs[obs_group].__array__()))
    unique_group_top_n_dges = []
    for group in unique_groups:
        score_df = sc.get.rank_genes_groups_df(adata, group = group)
        score_df_sorted = score_df.sort_values(["pvals_adj"], ascending = True)
        top_n_dges = score_df_sorted[0:n]["names"].__array__()
        unique_group_top_n_dges.append(top_n_dges)
    unique_groups_long = np.repeat(unique_groups, n)
    
    group_dges_df_n = pd.DataFrame({
        "Cluster": unique_groups_long,
        "Top {n} DGEs".format(n = n): np.concatenate(unique_group_top_n_dges)
    })
    return group_dges_df_n

def set_concordance(*args):
    """Determines number of overlapping elements between n sets.

    Args:
        *args: A list of n sets.

    Returns:
        concordance (int): The number of overlapping elements between n sets.
    """
    concordance = len(set.intersection(*args))
    return concordance
