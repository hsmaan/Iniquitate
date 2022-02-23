import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc

from utils import dge_top_n, diffexp

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    k_initial = adata.uns["kmeans_stats"]["kmeans_initial_k"]
    k_final = adata.uns["kmeans_stats"]["kmeans_final_k"]
    
    # Check if k_final is 1 and if so, skip DGE
    if k_final == 1:
        # Create and save summary dataframe for DGE results
        dge_int_summary_df = pd.DataFrame({
                "Dataset": dataset_name,
                "Batches downsampled": num_batches_ds,
                "Number of celltypes downsampled": num_celltypes_ds,
                "Proportion downsampled": prop_ds,
                "Replicate": rep,
                "Cluster number before convergence": k_initial,
                "Cluster number after convergence": k_final,
                "Method 1": "NA",
                "Method 2": "NA",
                "DGE Set Intersection Ratio": "NA",
                "Median DGE Set Intersection Ratio": "NA"
            },
            index = [0]
        )
        dge_int_summary_df.to_csv(save_loc, sep = "\t", index = False)
    else:
        # Subset adatas based on method for integration and store lognorm counts in raw
        # attribute for diffexp testing
        methods = ["harmony", "scvi", "scanorama", "seurat", "liger"]
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
                
        # Extract top 50 DGEs for each cluster in each method 
        method_dge_dfs = []
        for adata_method_subset in method_adatas:
            adata_method_subset = diffexp(
                adata_method_subset, 
                groupby = "kmeans_faiss",
                use_raw = True,
                method = "wilcoxon"
            )
            dge_results = dge_top_n(
                adata_method_subset, 
                n = 50,
                obs_group = "kmeans_faiss"
            )
            method_dge_dfs.append(dge_results)
            
        # Concatenate DGE results from each method 
        method_dge_dfs_concat = pd.concat(method_dge_dfs, axis = 0)
        
        # Create long form array for methods 
        methods_long = np.repeat(np.array(methods), 50*k_final)
        
        # Create summary dataframe for DGE results
        dge_summary_df = pd.DataFrame({
            "Dataset": dataset_name,
            "Batches downsampled": num_batches_ds,
            "Number of celltypes downsampled": num_celltypes_ds,
            "Proportion downsampled": prop_ds,
            "Replicate": rep,
            "Cluster number before convergence": k_initial,
            "Cluster number after convergence": k_final,
            "Method": methods_long,
            "Cluster": method_dge_dfs_concat["Cluster"].__array__(),
            "Differentially expressed genes": method_dge_dfs_concat["Top 50 DGEs"].__array__()
        })
        
        # Determine DGE concordance through set intersection in a pairwise manner 
        method_concordance_mat = np.zeros((len(methods), len(methods)))
        for i, method_i in enumerate(methods):
            for j, method_j in enumerate(methods):
                dge_sub_1 = dge_summary_df[dge_summary_df["Method"] == method_i]
                dge_sub_2 = dge_summary_df[dge_summary_df["Method"] == method_j]
                dge_sub_1_genes = dge_sub_1["Differentially expressed genes"].values
                dge_sub_2_genes = dge_sub_2["Differentially expressed genes"].values
                set_int = np.intersect1d(dge_sub_1_genes, dge_sub_2_genes)
                int_ratio = len(set_int)/len(dge_sub_1_genes)
                method_concordance_mat[i, j] = int_ratio
        
                
        # Create dataframe of values
        method_int_df = pd.DataFrame(method_concordance_mat)
        method_int_df.index = methods
        method_int_df.columns = methods
        
        # Convert to long format 
        method_int_df_long = method_int_df.melt(ignore_index = False)
        method_int_df_long = method_int_df_long.reset_index()
        method_int_df_long.columns = ["Method 1", "Method 2", "DGE Set Intersection Ratio"] 

        # Get median of DGE concordance
        method_int_df_no_self = method_int_df_long[method_int_df_long["Method 1"] != method_int_df_long["Method 2"]]
        median_set_int_ratio = np.median(method_int_df_no_self["DGE Set Intersection Ratio"])
        method_int_df_long["Median DGE Set Intersection Ratio"] = median_set_int_ratio
        
        # Create and save summary dataframe for DGE intersection results
        dge_int_summary_df = pd.DataFrame({
            "Dataset": dataset_name,
            "Batches downsampled": num_batches_ds,
            "Number of celltypes downsampled": num_celltypes_ds,
            "Proportion downsampled": prop_ds,
            "Replicate": rep,
            "Cluster number before convergence": k_initial,
            "Cluster number after convergence": k_final,
            "Method 1": method_int_df_long["Method 1"].__array__(),
            "Method 2": method_int_df_long["Method 2"].__array__(),
            "DGE Set Intersection Ratio": method_int_df_long["DGE Set Intersection Ratio"].__array__(),
            "Median DGE Set Intersection Ratio": method_int_df_long["Median DGE Set Intersection Ratio"].__array__()
        })
        dge_int_summary_df.to_csv(save_loc, sep = "\t", index = False)
        
    
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