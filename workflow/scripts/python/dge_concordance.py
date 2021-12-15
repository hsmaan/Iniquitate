import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc

from utils import faiss_kmeans, dge_top_n, diffexp

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    
    # Define method subsets and iterate over them until the same number of k clusters is found
    k = 10
    k_initial = k
    methods = ["harmony", "scvi", "bbknn", "scanorama", "seurat", "liger"]
    method_dge_adatas = []
    i = 0
    while i < len(methods):
        # Define method subset
        adata_subset = adata[adata.obs["integration_method"] == methods[i]]
        
        # Perform HVG selection on raw (unnormalized, unlogged) data
        adata_subset.X = adata_subset.layers["raw"]
        sc.pp.normalize_total(
            adata_subset,
            target_sum = 1e4
        )
        sc.pp.log1p(adata_subset)
        sc.pp.highly_variable_genes(
            adata_subset,
            n_top_genes = 2500,
            flavor = "seurat"
        )
        
        # Store raw attribute (lognorm counts) for later differential expression analysis
        adata_subset.raw = adata_subset
        
        # Perform faiss kmeans clustering
        adata_subset, k_method = faiss_kmeans(adata_subset, k)
        
        # Test concordance of k values and either append or reset
        if k_method != k:
            k = k_method
            i = 0 
            method_dge_adatas.clear()
            continue 
        else:
            i += 1
            method_dge_adatas.append(adata_subset)
            
    # Extract top 50 DGEs for each cluster in each method 
    method_dge_dfs = []
    for adata in method_dge_adatas:
        adata = diffexp(
            adata, 
            groupby = "kmeans_faiss",
            use_raw = True,
            method = "wilcoxon"
        )
        dge_results = dge_top_n(
            adata, 
            n = 50,
            obs_group = "kmeans_faiss"
        )
        method_dge_dfs.append(dge_results)
        
    # Concatenate DGE results from each method 
    method_dge_dfs_concat = pd.concat(method_dge_dfs, axis = 0)
    
    # Create long form array for methods 
    methods_long = np.repeat(np.array(methods), 50*k)
    
    # Create and save summary dataframe for DGE results
    dge_summary_df = pd.DataFrame({
        "Dataset": dataset_name,
        "Batches downsampled": num_batches_ds,
        "Number of celltypes downsampled": num_celltypes_ds,
        "Proportion downsampled": prop_ds,
        "Replicate": rep,
        "Cluster number before convergence": k_initial,
        "Cluster number after convergence": k,
        "Method": methods_long,
        "Cluster": method_dge_dfs_concat["Cluster"].__array__(),
        "Differentially expressed genes": method_dge_dfs_concat["Top 50 DGEs"].__array__()
    })
    dge_summary_df.to_csv(save_loc, sep = "\t", index = False)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for downsampling summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving downsampling statistics of integrated h5ad file"
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