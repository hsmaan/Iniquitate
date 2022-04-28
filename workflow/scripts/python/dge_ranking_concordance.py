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
    all_genes = np.sort(adata.var.index.values)
    all_genes_tiled = np.tile(all_genes, (len(methods)))
    methods_repeat = np.repeat(methods, len(all_genes))
    gene_max_imp_per_method = []
    gene_min_imp_per_method = []
    for method_dge_df in method_dge_dfs:
        method_clusters = np.unique(method_dge_df["Cluster"].__array__())
        cluster_ranks = []
        for cluster in method_clusters:
            cluster_sub = method_dge_df[method_dge_df["Cluster"] == cluster]
            gene_ranks_sorted = np.argsort(cluster_sub.iloc[:, 1].__array__())
            cluster_ranks.append(gene_ranks_sorted)
        cluster_ranks_stack = np.stack(cluster_ranks, axis = 0)
        # Min for max because lowest number for ranking corresponds to highest importance
        gene_max_imp_per_method.append(
            np.min(cluster_ranks_stack, axis = 0)
        )
        gene_min_imp_per_method.append(
            np.max(cluster_ranks_stack, axis = 0)
        )
        
    # Concatenate max and min ranks for each method 
    gene_max_imp_per_method_concat = np.concatenate(gene_max_imp_per_method)
    gene_min_imp_per_method_concat = np.concatenate(gene_min_imp_per_method)
    
    # Create summary df of all genes and their ranking metrics
    dge_ranking_summary_df = pd.DataFrame({
        "Dataset": dataset_name,
        "Number of batches downsampled": num_batches_ds,
        "Number of celltypes downsampled": num_celltypes_ds,
        "Proportion downsampled": prop_ds,
        "Replicate": rep,
        "Method": methods_repeat,
        "Gene": all_genes_tiled,
        "Max rank": gene_max_imp_per_method_concat,
        "Min rank": gene_min_imp_per_method_concat
    })
    dge_ranking_summary_df.to_csv(save_loc, sep = "\t", index = False)
    
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