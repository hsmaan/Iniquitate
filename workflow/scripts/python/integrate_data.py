import argparse 
import os 
import sys 
sys.path.append("src/python/")

import scanpy as sc
import anndata as ann
import numpy as np

from utils import Integration, downsample

def main(h5ad_dir, save_loc, ds_celltypes, ds_proportions, num_batches):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f))
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
        for adata in adata_selected:
            adata_ds = downsample(
                adata = adata, 
                num_celltypes = ds_celltypes, 
                proportion = ds_proportions
            )
            adata_downsampled.append(adata_ds)
        adata_loaded = adata_unselected + adata_downsampled

    # Concatenate files (assume data is raw counts)
    adata_concat = ann.AnnData.concatenate(*adata_loaded)
    
    # Create integration class instance 
    integration = Integration(adata = adata_concat)
    
    # Integrate across subsets
    harmony_integrated = integration.harmony_integrate()
    scvi_integrated = integration.scvi_integrate()
    bbknn_integrated = integration.bbknn_integrate()
    scanorama_integrated = integration.scanorama_integrate()
    seurat_integrated = integration.seurat_integrate()
    # liger_integrated = integration.liger_integrate()
    
    # Add integration type to each subset and concatenate and save
    harmony_integrated.obs["integration_method"] = "harmony" 
    scvi_integrated.obs["integration_method"] = "scvi"
    bbknn_integrated.obs["integration_method"] = "bbknn"
    scanorama_integrated.obs["integration_method"] = "scanorama"
    seurat_integrated.obs["integration_method"] = "seurat"
    # liger_integrated.obs["integration_method"] = "liger"
    
    integrated_concat = ann.concat([
        harmony_integrated,
        scvi_integrated,
        bbknn_integrated,
        scanorama_integrated,
        seurat_integrated,
        # liger_integrated
    ])
    
    integrated_concat.write_h5ad(
        filename = save_loc
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