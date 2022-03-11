import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc

from utils import dge_top_n

def main(h5ad_dir, save_loc, top_n = 10):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f), as_sparse = "raw/X")
        adata.layers["raw"] = adata.X # Store raw counts
        adata.obs = adata.obs[["batch", "celltype"]] # Only store relevant columns
        if "gene" not in adata.var.columns:
            adata.var["gene"] = adata.var_names # Add gene names if not present
        adata.var = adata.var[["gene"]] # Only store relevant columns
        adata_loaded.append(adata)
        
    # Get the differential gene expression results for the celltypes in each batch, top n
    # for each celltype and get and return the union
    adata_dge_top_n_dfs = []
    for adata in adata_loaded:
        # Log-normalize the data
        sc.pp.normalize_total(
            adata,
            target_sum = 1e4
        )
        sc.pp.log1p(adata)
        
        # Store lognorm counts and perform DGE based on celltype
        adata.raw = adata # Freeze for DGE test - lognorm counts
        sc.tl.rank_genes_groups(
            adata, 
            groupby = "celltype",
            use_raw = True,
            method = "wilcoxin" 
        )
        
        # Get the top n degs for each celltype and append to all results 
        dge_results = dge_top_n(
            adata, 
            n = top_n,
            obs_group = "celltype"
        )
        adata_dge_top_n_dfs.append(dge_results)
        
    # Concatenate all dge dataframes and keep distinct rows 
    adata_dge_top_n_concat = pd.concat(adata_dge_top_n_dfs)
    adata_dge_top_n_concat = adata_dge_top_n_concat.drop_duplicates()
    
    # Rename columns appropriately and save 
    adata_dge_top_n_concat.columns = [
        "Celltype",
        "Top {n} marker genes (union across batches)".format(
            n = top_n
        )
    ]
    adata_dge_top_n_concat.to_csv(save_loc, sep = "\t", index = False)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for marker gene summary"
    )
    parser.add_argument(
        "--filedir",
        type = str,
        help = "Path of directory containing scRNA-seq h5ad files"
    )        
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving output from marker gene selection"
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        save_loc = args.outfile,
    )