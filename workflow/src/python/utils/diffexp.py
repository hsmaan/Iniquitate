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

def dge_top_50(adata):
    """
    Return clusters and genes with the top 50 differential expression.

    Args:
        adata (AnnData): Annotated data matrix object.

    Returns:
        data (DataFrame): A dataframe with the top 50 differentially expressed genes in each cluster.
    """
    unique_groups = np.sort(np.unique(adata.obs.leiden.__array__()))
    unique_group_top_50_dges = []
    for group in unique_groups:
        score_df = sc.get.rank_genes_groups_df(adata, group = group)
        score_df_sorted = score_df.sort_values(["pvals_adj"], ascending = True)
        top_50_dges = score_df_sorted[0:50]["names"].__array__()
        unique_group_top_50_dges.append(top_50_dges)
    unique_groups_long = np.repeat(unique_groups, 50)
    
    group_dges_df_50 = pd.DataFrame({
        "Cluster": unique_groups_long,
        "Top 50 DGEs": np.concatenate(unique_group_top_50_dges)
    })
    return group_dges_df_50

def set_concordance(*args):
    """Determines number of overlapping elements between n sets.

    Args:
        *args: A list of n sets.

    Returns:
        concordance (int): The number of overlapping elements between n sets.
    """
    concordance = len(set.intersection(*args))
    return concordance