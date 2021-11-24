import argparse 
import os 
import sys 
sys.path.append("src/python/")

import scanpy as sc
import anndata as ann
import pandas as pd 

from utils import Integration

def main(h5ad_dir, save_loc):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f))
        adata_loaded.append(adata)
    
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
        "--outfile",
        type = str,
        help = "Filepath for saving output from scRNA-seq integration"
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        save_loc = args.outfile
    )