import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc
from utils import cluster_concordance 

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    
    # Get clustering concordance results 
    cluster_concordance_df = cluster_concordance(adata = adata)
    
    # Create cluster concordance summary df 
    cluster_concordance_summary_df = pd.DataFrame({
        "Dataset": dataset_name,
        "Number of batches downsampled": num_batches_ds,
        "Number of celltypes downsampled": num_celltypes_ds,
        "Proportion downsampled": prop_ds,
        "Replicate": rep,
        "Method 1": cluster_concordance_df["Method 1"].__array__(),
        "Method 2": cluster_concordance_df["Method 2"].__array__(),
        "ARI": cluster_concordance_df["ARI"].__array__(),
        "Median ARI": cluster_concordance_df["Median ARI"].__array__()
    }) 
        
    # Save clustering concordance dataframe to tsv
    cluster_concordance_summary_df.to_csv(
        save_loc,
        index=False,
        sep="\t"
    )

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