import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import scanpy as sc

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
        
    # Determine pearson, spearman, and kendall correlations between post-integration PAGA estimated pseudotime and 
    # pre-integration pseudotime for each batch-correction method
    spearman_corrs = []
    pearson_corrs = []
    kendall_corrs = []
    for adata_sub in adata_method_sub:
        pass