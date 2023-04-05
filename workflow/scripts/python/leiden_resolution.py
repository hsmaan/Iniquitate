import argparse 
import os 
import sys 
sys.path.append("src/python/")
os.environ['CUDA_VISIBLE_DEVICES'] = "0, 1"

import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd

from utils import Integration

def none_or_str(value):
    if value == 'None':
        return None
    return value

def main(h5ad_dir, save_loc):
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
        
    # Create ranges of leiden resolution to test for each method
    leiden_ranges = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    # Store batch name separately for each anndata object
    for adata in adata_loaded:
        adata.obs["batch_name"] = adata.obs["batch"]

    # Concatenate files (assume data is raw counts)
    adata_concat = ann.AnnData.concatenate(*adata_loaded)
    adata_concat.obs_names = range(len(adata_concat.obs_names))
    adata_concat.obs_names_make_unique()
    adata_concat.obs["batch"] = adata_concat.obs["batch_name"]
    adata_concat.obs.drop("batch_name", axis = 1, inplace = True)
    
    # Loop over leiden ranges and perform integration using custom ranges
    # for 20 iterations 
    iters = 0
    optimal_resolutions = [] 
    while iters < 20:
        diffs_list = []
        for leiden_range in leiden_ranges:
            # Integrate the data in the given range
            leiden_resolutions = [leiden_range for i in range(6)]
            integration = Integration(adata = adata_concat.copy(), custom_resolutions = leiden_resolutions)

            scvi_integrated = integration.scvi_integrate()            
            harmony_integrated = integration.harmony_integrate()
            bbknn_integrated = integration.bbknn_integrate()
            scanorama_integrated = integration.scanorama_integrate()
            seurat_integrated = integration.seurat_integrate()
            liger_integrated = integration.liger_integrate()
            
            # determine the difference between the actual number of cell-types 
            n_celltypes = len(np.unique(adata_concat.obs["celltype"].values))
            n_clusters_scvi = len(np.unique(scvi_integrated.obs["leiden"].values))
            n_clusters_harmony = len(np.unique(harmony_integrated.obs["leiden"].values))
            n_clusters_bbknn = len(np.unique(bbknn_integrated.obs["leiden"].values))
            n_clusters_scanorama = len(np.unique(scanorama_integrated.obs["leiden"].values))
            n_clusters_seurat = len(np.unique(seurat_integrated.obs["leiden"].values))
            n_clusters_liger = len(np.unique(liger_integrated.obs["leiden"].values))
            
            # Get the absolute difference between the number of cell-types and the number of clusters
            diff_scvi = abs(n_celltypes - n_clusters_scvi)
            diff_harmony = abs(n_celltypes - n_clusters_harmony)
            diff_bbknn = abs(n_celltypes - n_clusters_bbknn)
            diff_scanorama = abs(n_celltypes - n_clusters_scanorama)
            diff_seurat = abs(n_celltypes - n_clusters_seurat)
            diff_liger = abs(n_celltypes - n_clusters_liger)
            
            # Create a list of the values and append to the list of differences
            diffs = [diff_scvi, diff_harmony, diff_bbknn, diff_scanorama, diff_seurat, diff_liger]
            diffs_list.append(diffs)
        
        # Get the resolution that minimizes the difference between the number of cell-types and
        # the number of clusters for each method
        diffs_arr = np.transpose(np.array(diffs_list))
        min_diffs = np.argmin(diffs_arr, axis = 0)
        resolutions = [leiden_ranges[i] for i in min_diffs]
        optimal_resolutions.append(resolutions)
        
    # Take the median of the optimal resolutions for each method 
    optimal_res_median = np.median(np.array(optimal_resolutions), axis = 0)
    
    # Create a dataframe of the optimal resolutions for each method
    methods = ["scvi", "harmony", "bbknn", "scanorama", "seurat", "liger"]
    opt_res_df = pd.DataFrame({"method": methods, "resolution": optimal_res_median})
    opt_res_df.to_csv(save_loc, index=False, sep="\t")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for scRNA-seq integration"
    )
    parser.add_argument(
        "--filedir",
        type = str,
        help = "Path of directory containing scRNA-seq h5ad files"
    )
    parser.add_argument(
        "--outfile",
       type = str,
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        save_loc = args.outfile
    )