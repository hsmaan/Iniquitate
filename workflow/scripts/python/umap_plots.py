import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc 

from utils import Umap

def main(h5ad_loc, leiden_save_loc, celltype_save_loc, batch_save_loc):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Get the umap coordinates for all the methods and create dictionary object
    methods = ["bbknn", "harmony", "scanorama", "scvi", "seurat"]
    umap_dict = {}
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        umap_dict[method] = adata_sub.obsm["X_umap"].__array__()
        
    # Get the leiden clustering for all the methods and create dictionary object
    leiden_dict = {}
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        leiden_dict[method] = adata_sub.obs["leiden"].__array__()
        
    # Create a cell type dictionary object
    celltype_dict = {}
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        celltype_dict[method] = adata_sub.obs["celltype"].__array__()
        
    # Create a batch dictionary object
    batch_dict = {}
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        batch_dict[method] = adata_sub.obs["batch"].__array__()
        
    # Create a umap object for each subset of information
    umap_leiden = Umap(
        coords = umap_dict,
        clustering = leiden_dict,
        subset_name = "Clustering",
    )
    umap_celltype = Umap(
        coords = umap_dict,
        clustering = celltype_dict,
        subset_name = "Cell-type",
    )
    umap_batch = Umap(
        coords = umap_dict,
        clustering = batch_dict,
        subset_name = "Batch",
    )
    
    # Plot each of the unap objects
    umap_leiden.umap_df()
    umap_leiden.umap_plot(show_plot=True)
    umap_leiden.save_umap(
        save_dir=leiden_save_loc,
        dpi=300
    )
    
    umap_celltype.umap_df()
    umap_celltype.umap_plot(show_plot=True)
    umap_celltype.save_umap(
        save_dir=celltype_save_loc,
        dpi=300
    )
    
    umap_batch.umap_df()
    umap_batch.umap_plot(show_plot=True)
    umap_batch.save_umap(
        save_dir=batch_save_loc,
        dpi=300
    )
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for UMAP plot generation"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--leiden_plot_loc",
        type = str,
        help = "Filepath for saving leiden overlayed UMAP results"
    )
    parser.add_argument(
        "--celltype_plot_loc",
        type = str,
        help = "Filepath for saving celltype overlayed UMAP results"
    )
    parser.add_argument(
        "--batch_plot_loc",
        type = str,
        help = "Filepath for saving batch overlayed UMAP results"
    )
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        leiden_save_loc = args.leiden_plot_loc,
        celltype_save_loc = args.celltype_plot_loc,
        batch_save_loc = args.batch_plot_loc
    )