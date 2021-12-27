import argparse 
import os 
import sys 
sys.path.append("src/python/")
os.environ['CUDA_VISIBLE_DEVICES'] = "0, 1"

import scanpy as sc
import anndata as ann
import numpy as np

from utils import Integration, downsample, faiss_kmeans

def main(h5ad_dir, save_loc, ds_celltypes, ds_proportions, num_batches):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f), as_sparse = "raw/X")
        adata.layers["raw"] = adata.X # Store raw counts
        adata_loaded.append(adata)
    
    # Downsample loaded h5ad files based on params 
    if num_batches == 0:
        pass
    else:
        selected_indices = np.random.choice(
            len(adata_loaded), num_batches, replace = False
        )
        adata_selected = [adata_loaded[i] for i in selected_indices]
        adata_unselected = [adata_loaded[i] for i in range(len(adata_loaded)) if i not in selected_indices]
        adata_downsampled = []
        selected_celltypes_downsampled = []
        for adata in adata_selected:
            adata_ds, selected_celltypes_ds = downsample(
                adata = adata, 
                num_celltypes = ds_celltypes, 
                proportion = ds_proportions
            )
            adata_downsampled.append(adata_ds)
            selected_celltypes_downsampled.append(selected_celltypes_ds)
        adata_loaded = adata_unselected + adata_downsampled

    # Store batch name separately for each anndata object
    for adata in adata_loaded:
        adata.obs["batch_name"] = adata.obs["batch"]

    # Concatenate files (assume data is raw counts)
    adata_concat = ann.AnnData.concatenate(*adata_loaded)
    adata_concat.obs_names = range(len(adata_concat.obs_names))
    adata_concat.obs_names_make_unique()
    adata_concat.obs["batch"] = adata_concat.obs["batch_name"]
    adata_concat.obs.drop("batch_name", axis = 1, inplace = True)
    
    # Create integration class instance 
    integration = Integration(adata = adata_concat, gpu = False)
    
    # Integrate across subsets
    harmony_integrated = integration.harmony_integrate()
    bbknn_integrated = integration.bbknn_integrate()
    scanorama_integrated = integration.scanorama_integrate()
    seurat_integrated = integration.seurat_integrate()
    liger_integrated = integration.liger_integrate()
    
    # Add integration type to each subset and concatenate
    harmony_integrated.obs["integration_method"] = "harmony" 
    bbknn_integrated.obs["integration_method"] = "bbknn"
    scanorama_integrated.obs["integration_method"] = "scanorama"
    seurat_integrated.obs["integration_method"] = "seurat"
    liger_integrated.obs["integration_method"] = "liger"
    
    integrated_concat = ann.concat([
        harmony_integrated,
        bbknn_integrated,
        scanorama_integrated,
        seurat_integrated,
        liger_integrated
    ])
    integrated_concat.obs_names = range(len(integrated_concat.obs_names))
    integrated_concat.obs_names_make_unique()
    
    # Add placeholder in entire obs dataframe for kmeans clustering
    integrated_concat.obs["kmeans_faiss"] = np.zeros(len(integrated_concat.obs_names))
    
    # Perform kmeans clustering on integrated data 
    # Define method subsets and iterate over them until the same number of k clusters is found
    k = 10
    k_initial = k # Integers are immutable 
    methods = ["harmony", "scanorama", "seurat", "liger"]
    method_kmeans_adatas = []
    i = 0
    while i < len(methods):
        # Create a copy of adata to avoid overwriting the original
        adata_copy = integrated_concat.copy()
        
        # Define method subset
        adata_subset = adata_copy[adata_copy.obs["integration_method"] == methods[i]]
        
        # Perform HVG selection on raw (unnormalized, unlogged) data
        adata_subset.X = adata_subset.layers["raw"]
        sc.pp.normalize_total(
            adata_subset,
            target_sum = 1e4
        )
        sc.pp.log1p(adata_subset)
        sc.pp.highly_variable_genes(
            adata_subset,
            n_top_genes = 2500,
            flavor = "seurat"
        )
        
        # Store raw attribute (lognorm counts) for later differential expression analysis
        adata_subset.raw = adata_subset
        
        # Perform faiss kmeans clustering
        adata_subset, k_method = faiss_kmeans(adata_subset, k)
        
        # Test concordance of k values and either append or reset
        if k_method != k:
            k = k_method
            i = 0 
            method_kmeans_adatas.clear()
            continue 
        else:
            i += 1
            method_kmeans_adatas.append(adata_subset)
    
    # Append kmeans cluster info to integrated data
    for method, method_kmeans_adata in zip(methods, method_kmeans_adatas):
        method_kmeans_clusters = method_kmeans_adata.obs["kmeans_faiss"].__array__().astype('str')
        integrated_concat.obs.loc[
            integrated_concat.obs["integration_method"] == method,
            "kmeans_faiss"
        ] = method_kmeans_clusters
        
    # Add placeholder for bbknn kmeans clustering
    integrated_concat.obs.loc[
        integrated_concat.obs["integration_method"] == "bbknn",
        "kmeans_faiss"
    ] = "NA"
    
    # Append information about kmeans faiss clusters to .uns of adata_concat
    integrated_concat.uns["kmeans_stats"] = {
        "kmeans_initial_k": k_initial,
        "kmeans_final_k": k
    }

    # Add data about downsampling to .uns of adata_concat
    if num_batches == 0:
        integrated_concat.uns["downsampling_stats"] = {
            "num_batches": 0,
            "num_celltypes_downsampled": 0,
            "ds_batch_names": "None",
            "proportion_downsampled": 1,
            "downsampled_celltypes": "None"
        }
    else:
        integrated_concat.uns["downsampling_stats"] = {
            "num_batches": num_batches,
            "num_celltypes_downsampled": ds_celltypes,
            "ds_batch_names": np.concatenate([np.unique(adata.obs["batch"].__array__()) for adata in adata_downsampled]),
            "proportion_downsampled": ds_proportions,
            "downsampled_celltypes": selected_celltypes_downsampled
        }
        
    # Save integrated h5ad object
    integrated_concat.write_h5ad(
        filename = save_loc,
        compression = "gzip"
    )
    
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
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving output from scRNA-seq integration"
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        save_loc = args.outfile,
        ds_celltypes = args.ds_celltypes,
        ds_proportions = args.ds_proportions,
        num_batches = args.num_batches        
    )
