import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc
from sklearn import metrics

def none_or_str(value):
    if value == 'None':
        return None
    return value

def main(h5ad_loc, save_loc, annofile, dataset_name, ds_celltypes, ds_celltypes_names, 
         ds_proportions, num_batches, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Load acceptable annotations
    annos = pd.read_csv(annofile, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for annotation results summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of Seurat annotated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving annotation analysis and scoring results"
    )
    parser.add_argument(
        "--annofile",
        type = str,
        help = "Filepath for acceptable annotation results/matched with actual celltypes"
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
    parser.add_argument(
        "--ds_celltypes",
        type = int,
        help = "Number of celltypes to randomly downsample in given batch"
    )
    parser.add_argument(
        "--ds_celltype_names",
        type = none_or_str,
        nargs = "?",
        default = None,
        help = "Custom names of celltypes to downsample in given batch"
    )
    parser.add_argument(
        "--ds_proportions",
        type = float,
        help = "Proportion of downsampling per celltype in a given batch"
    )
    parser.add_argument(
        "--num_batches",
        type = int,
        help = "Number of batches to perform downsampling on"
    )
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        save_loc = args.outfile,
        annofile = args.annofile,
        dataset_name = args.dataset,
        rep = args.rep,
        ds_celltypes = args.ds_celltypes,
        ds_celltypes_names = args.ds_celltype_names,
        ds_proportions = args.ds_proportions,
        num_batches = args.num_batches  
    )