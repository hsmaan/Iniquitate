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

    # Store lognorm counts in raw attribute for diffexp testing
    adata.raw = adata.copy()

    # Extract list of all DGEs for all leiden clusters 
    adata = diffexp(
        adata,
        groupby = "leiden",
        use_raw = True,
        method = "wilcoxon"
    )
    dge_results = dge_top_n(
        adata,
        n = len(adata.var),
        obs_group = "leiden"
    )
    
    # Compute the ranking metrics for all genes in the dataset
    all_genes = np.sort(adata.var.index.values)
    gene_max_imp = []
    gene_min_imp = []
    clusters = np.unique(dge_results["Cluster"].__array__())
    cluster_ranks = []
    for cluster in clusters:
        cluster_sub = dge_results[dge_results["Cluster"] == cluster]
        gene_ranks_sorted = np.argsort(cluster_sub.iloc[:, 1].__array__())
        cluster_ranks.append(gene_ranks_sorted)
    cluster_ranks_stack = np.stack(cluster_ranks, axis = 0)
    gene_max_imp = np.min(cluster_ranks_stack, axis = 0)
    gene_min_imp = np.max(cluster_ranks_stack, axis = 0)
    
    # Create summary df of all genes and their ranking metrics
    dge_ranking_summary_df = pd.DataFrame({
        "Dataset": dataset_name,
        "Number of batches downsampled": num_batches_ds,
        "Number of celltypes downsampled": num_celltypes_ds,
        "Proportion downsampled": prop_ds,
        "Replicate": rep,
        "Gene": all_genes,
        "Max rank": gene_max_imp,
        "Min rank": gene_min_imp
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