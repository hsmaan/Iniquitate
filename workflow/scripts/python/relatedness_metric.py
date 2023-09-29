import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import scanpy as sc

from utils import relatedness_score

def main(h5ad_dir, save_loc, dataset_name):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f), as_sparse = "raw/X")
        adata_loaded.append(adata)
        
    # Get relatedness metric for each celltype within each batch 
    celltype_relatedness_dfs = []
    for adata in adata_loaded:
        celltype_relatedness_df = relatedness_score(adata, pca_performed = False)
        celltype_relatedness_dfs.append(celltype_relatedness_df)
        
    # Concatenate results, add relevant metadata and save
    celltype_relatedness_dfs_concat = pd.concat(celltype_relatedness_dfs)
    celltype_relatedness_dfs_concat["Dataset"] = dataset_name
    celltype_relatedness_dfs_concat.to_csv(save_loc, sep = "\t", index = False)
         
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for calculating relatedness metric"
    )
    parser.add_argument(
        "--filedir",
        type = str,
        help = "Path of directory containing scRNA-seq h5ad files"
    )        
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving output from relatedness metric calculation"
    )
    parser.add_argument(
        "--dataset",
        type = str,
        help = "Name of dataset"
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        save_loc = args.outfile,
        dataset_name = args.dataset
    )