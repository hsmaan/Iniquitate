import argparse 
import os 
import sys 

import scanpy as sc
import anndata as ann
import numpy as np

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    
    # Subset h5ad based on batch-correction method used
    adata_method_sub = []
    methods = ["harmony", "scvi", "scanorama", "seurat", "liger"] # Omitting BBKNN due to lack of embedding
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        adata_method_sub.append(
            adata_sub
        )
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for clustering results summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving clustering results of integrated h5ad file"
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