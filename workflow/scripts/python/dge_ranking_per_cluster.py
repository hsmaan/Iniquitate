import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc

from utils import diffexp, dge_top_n

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]

    # Subset adatas based on method for integration and store lognorm counts in raw
    # attribute for diffexp testing
    methods = ["harmony", "scvi", "scanorama", "seurat", "liger", "bbknn"]
    method_adatas = []
    for method in methods:
        adata_copy = adata.copy()
        adata_subset = adata_copy[adata_copy.obs["integration_method"] == method]
        adata_subset.X = adata_subset.layers["raw"] # Unlogged, unnorm counts
        sc.pp.normalize_total(
            adata_subset,
            target_sum = 1e4
        )
        sc.pp.log1p(adata_subset)
        adata_subset.raw = adata_subset # Freeze for DGE test 
        method_adatas.append(adata_subset)
        
    # Extract list of all DGEs for all leiden clusters in each method 
    method_dge_dfs = []
    for adata_method_subset in method_adatas:
        adata_method_subset = diffexp(
            adata_method_subset, 
            groupby = "leiden",
            use_raw = True,
            method = "wilcoxon"
        )
        dge_results = dge_top_n(
            adata_method_subset, 
            n = len(adata.var),
            obs_group = "leiden"
        )
        method_dge_dfs.append(dge_results)
    
    # For each method, compute the ranking metrics for all genes in the dataset
    # based on each cluster   
    all_genes = np.sort(adata.var.index.values)
    method_adata_result_dfs = []
    for method_adata, method_dge_df in zip(method_adatas, method_dge_dfs):
        method_clusters = np.unique(method_dge_df["Cluster"].__array__())
        method_name = np.unique(method_adata.obs["integration_method"].__array__())
        cluster_ranks = []
        cluster_celltype = []
        cluster_number = []
        for cluster in method_clusters:
            cluster_celltype_unique = np.unique(
                method_adata["celltype"][method_adata.obs["leiden"] == cluster],
                return_counts = True
            )
            celltype_most_prev = cluster_celltype_unique[0][
                np.argmax(cluster_celltype_unique[1])
            ]    
            cluster_sub = method_dge_df[method_dge_df["Cluster"] == cluster]
            gene_ranks_sorted_50 = np.argsort(cluster_sub.iloc[:, 1].__array__())[0:50]
            cluster_ranks.append(gene_ranks_sorted_50)
            cluster_celltype.append(np.repeat(celltype_most_prev, 50))
            cluster_number.append(np.repeat(cluster, 50))
        cluster_ranks_full = np.concatenate(cluster_ranks)
        cluster_celltypes_full = np.concatenate(cluster_celltype)
        cluster_numbers_full = np.concatenate(cluster_number)
        method_adata_result = pd.DataFrame({
            "Cluster markers ranked (top 50)": cluster_ranks_full,
            "Cluster celltype (majority)": cluster_celltypes_full,
            "Cluster number": cluster_numbers_full,
            "Method": method_name
        })
        method_adata_result_dfs.append(method_adata_result)
        
    # Concatenate all results into one dataframe
    method_adata_result_df = pd.concat(method_adata_result_dfs)
    
    # Add all of the summary statistics to the dataframe and save
    method_adata_result_df["Dataset"] = dataset_name
    method_adata_result_df["Replicate"] = rep
    method_adata_result_df["Number of batches downsampled"] = num_batches_ds
    method_adata_result_df["Number of celltypes downsampled"] = num_celltypes_ds
    method_adata_result_df["Proportion downsampled"] = prop_ds

    method_adata_result_df.to_csv(save_loc, sep = "\t", index = False)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for dge ranking summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving dge ranking statistics of integrated h5ad file"
    )
    parser.add_argument(
        "--dataset",
        type = str,
        help = "Name of dataset"
    )
    parser.add_argument(
        "--rep",
        type = int,
        help = "Repetition number"
    )
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        save_loc = args.outfile,
        dataset_name = args.dataset,
        rep = args.rep
    )