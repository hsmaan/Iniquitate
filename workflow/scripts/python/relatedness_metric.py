import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import scanpy as sc

def main(h5ad_dir, save_loc, dataset_name):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f), as_sparse = "raw/X")
        
    # Get relatedness metric for each celltype within each batch 
    celltype_relatedness_dfs = []
    for adata in adata_loaded:
        batch_name = np.unique(adata.obs["batch"].__array__())[0]
        celltypes = np.unique(adata.obs["celltype"].__array__())
        celltype_i = []
        celltype_j = []
        relatedness_scores = []
        for celltype_i in celltypes:
            for celltype_j in celltypes:
                celltype_i.append(celltype_i)
                celltype_j.append(celltype_j)
                relatedness_scores.append(
                    relatedness(adata, celltype_i, celltype_j) # Placeholder for now 
                )
        celltype_relatedness_df = pd.DataFrame({
            "Celltype 1": celltype_i,
            "Celltype 2": celltype_j,
            "Relatedness": relatedness_scores,
            "Batch": batch_name
        })
        celltype_relatedness_dfs.append(celltype_relatedness_df)
        
    # Concatenate results, add relevant metadata and save
    celltype_relatedness_dfs_concat = pd.concat(celltype_relatedness_dfs)
    celltype_relatedness_dfs_concat["Dataset"] = dataset_name
    celltype_relatedness_dfs_concat.to_csv(save_loc, sep = "\t", index = False)
         
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