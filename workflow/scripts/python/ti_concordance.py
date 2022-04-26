import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.stats as sp 

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    batches_ds = adata.uns["downsampling_stats"]["ds_batch_names"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    downsampled_celltypes = adata.uns["downsampling_stats"]["downsampled_celltypes"]
    batches = np.unique(adata.obs.batch.__array__())
    
    # Concatenate downsampled celltypes and batches - conditioned on return type
    if isinstance(downsampled_celltypes, str):
        if downsampled_celltypes == "None":  
            downsampled_celltypes = "None"
        else:
            raise ValueError("Downsampled celltypes is a str and not 'None'")
    elif isinstance(downsampled_celltypes, list):
        downsampled_celltypes = np.array(downsampled_celltypes)
    elif isinstance(downsampled_celltypes, np.ndarray):
        if len(downsampled_celltypes.shape) == 1:
            downsampled_celltypes = downsampled_celltypes
        else:
            downsampled_celltypes = np.concatenate(downsampled_celltypes)
    else:
        raise TypeError("Downsampled celltypes is not a str, list, or ndarray")
    
    # Subset h5ad based on batch-correction method used
    adata_method_sub = []
    methods = ["harmony", "scvi", "bbknn", "scanorama", "seurat", "liger"]
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        adata_method_sub.append(
            adata_sub
        )
        
    # Determine pearson, spearman, and kendall correlations between post-integration PAGA estimated pseudotime and 
    # pre-integration pseudotime for each batch-correction method
    spearman_corrs = []
    pearson_corrs = []
    kendall_corrs = []
    for adata_sub in adata_method_sub:
        # Get pre-integration pseudotime
        preint_pt = adata_sub.obs["Sub_trajectory_Pseudotime"].__array__()
        
        # Get DPT pseudotime estimates
        dpt_pt = adata_sub.obs["dpt_pseudotime"].__array__()
        
        # Get correlations between pre-integration pseudotime and DPT pseudotime
        spearman_corr = sp.spearmanr(preint_pt, dpt_pt)[0]
        pearson_corr = sp.pearsonr(preint_pt, dpt_pt)[0]
        kendall_corr = sp.kendalltau(preint_pt, dpt_pt)[0]
        
        spearman_corrs.append(spearman_corr)
        pearson_corrs.append(pearson_corr)
        kendall_corrs.append(kendall_corr)
        
    # Create a dataframe with the results
    ti_corr_df = pd.DataFrame({
        "Method" : methods,
        "Spearman correlations" : spearman_corrs,
        "Pearson correlations" : pearson_corrs,
        "Kendall correlations" : kendall_corrs
    })
    ti_corr_df["Dataset"] = dataset_name
    ti_corr_df["Number of batches downsampled"] = num_batches_ds
    ti_corr_df["Number of celltypes downsampled"] = num_celltypes_ds
    ti_corr_df["Batches downsampled"] = batches_ds
    ti_corr_df["Proportion downsampled"] = prop_ds
    ti_corr_df["Downsampled celltypes"] = downsampled_celltypes
    ti_corr_df["Replicate"] = rep
    ti_corr_df["Total batches"] = len(batches)
    
    # Save dataframe to file 
    ti_corr_df.to_csv(
        save_loc,
        index = False,
        sep = "\t"
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for trajectory inference scoring"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of PAGA integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving trajectory inference scoring results"
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