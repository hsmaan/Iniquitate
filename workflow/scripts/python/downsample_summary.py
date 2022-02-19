import argparse 
import os 
import sys 

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc 

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    batches_ds = adata.uns["downsampling_stats"]["ds_batch_names"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    downsampled_celltypes = adata.uns["downsampling_stats"]["downsampled_celltypes"]
    
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
    
    # Extract data from just one integration method subset - for getting unique batches
    int_method_select = np.random.choice(
        np.unique(adata.obs.integration_method.__array__())
    )
    adata_select = adata[adata.obs.integration_method == int_method_select]
    
    # Create downsampling summary df 
    ds_summary_df = pd.DataFrame(
        {
            "Dataset": dataset_name,
            "Number of batches downsampled": num_batches_ds,
            "Batches downsampled": batches_ds,
            "Number of celltypes downsampled": num_celltypes_ds,
            "Proportion downsampled": prop_ds,
            "Downsampled celltypes": downsampled_celltypes,
            "Replicate": rep,
            "Total batches": len(np.unique(adata_select.obs["batch"]))
        }
    )
    ds_summary_df.to_csv(save_loc, index=False, sep="\t")
    
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