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
    downsampled_cells = adata.uns["downsampling_stats"]["downsampled_cells"]
    
    # Check if corresponding to 0 batch case 
    if num_batches_ds == 0:
        ds_summary_df = pd.DataFrame({
            "Dataset": dataset_name,
            "Batches downsampled": "None",
            "Number of celltypes downsampled": "None",
            "Proportion downsampled": "NA",
            "Batch label": "NA",
            "Downsampled celltypes": "NA",
            "Replicate": rep
        })
        ds_summary_df.to_csv(save_loc, index=False)
    
    else:
        # Concatenate downsampled cells and format data 
        downsampled_cells_concat = np.concatenate(downsampled_cells)
        batch_label = np.repeat(batches_ds, num_celltypes_ds)
        
        # Create downsampling summary df 
        ds_summary_df = pd.DataFrame({
            "Dataset": dataset_name,
            "Batches downsampled": num_batches_ds,
            "Number of celltypes downsampled": num_celltypes_ds,
            "Proportion downsampled": prop_ds,
            "Batch label": batch_label,
            "Downsampled celltypes": downsampled_cells_concat,
            "Replicate": rep
        })
        ds_summary_df.to_csv(save_loc, index=False)
    
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