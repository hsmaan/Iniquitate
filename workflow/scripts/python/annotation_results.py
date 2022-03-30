import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc
from sklearn.metrics import accuracy_score, balanced_accuracy_score, \
    f1_score, classification_report

def none_or_str(value):
    if value == 'None':
        return None
    return value

def main(h5ad_loc, save_loc, dataset_name, ds_celltypes, ds_proportions, 
         num_batches, rep):
    # Load h5ad file for query to reference mapping results
    adata = sc.read_h5ad(h5ad_loc)
    
    # Get the classification results as a dataframe 
    class_results = pd.DataFrame({
        "Real celltype": adata.obs["celltype"],
        "Predicted L1": adata.obs["predicted.celltype.l1"],
        "Predicted L2": adata.obs["predicted.celltype.l2"],
        "Control predicted L1": adata.obs["baseline.knn.l1"],
        "Control predicted L2": adata.obs["baseline.knn.l2"]
    }) 
    
    # Append information on dataset to results
    class_results["Dataset"] = dataset_name
    class_results["Number of batches downsampled"] = num_batches
    class_results["Number of celltypes downsampled"] = ds_celltypes
    class_results["Proportion downsampled"] = ds_proportions
    class_results["Replicate"] = rep
    
    # Save results to file
    class_results.to_csv(save_loc, index=False, sep="\t")

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
        help = "Filepath for saving annotation  results"
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
        dataset_name = args.dataset,
        rep = args.rep,
        ds_celltypes = args.ds_celltypes,
        ds_proportions = args.ds_proportions,
        num_batches = args.num_batches  
    )