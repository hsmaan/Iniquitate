import argparse 
import os 
import sys 
sys.path.append("src/python/")
os.environ['CUDA_VISIBLE_DEVICES'] = "0, 1"

import scanpy as sc
import anndata as ann
import numpy as np

from utils import IntegrationPAGA, downsample

def none_or_str(value):
    if value == 'None':
        return None
    return value

def main(h5ad_dir, root_celltype, save_loc, ds_celltypes, ds_celltypes_names, ds_proportions, 
         num_batches):
    # Load h5ad files 
    files_list = os.listdir(h5ad_dir)
    adata_loaded = []
    for f in files_list:
        adata = sc.read_h5ad(os.path.join(h5ad_dir, f), as_sparse = "raw/X")
        adata.layers["raw"] = adata.X # Store raw counts
        if "gene" not in adata.var.columns:
            adata.var["gene"] = adata.var_names # Add gene names if not present
        adata.var = adata.var[["gene"]] # Only store relevant columns
        adata.obs.celltype = [
            c.replace(" ", "_") for c in adata.obs.celltype
        ] # Remove spaces from celltype names - for Snakemake wildcard matching
        adata_loaded.append(adata)
        
    # Parse celltypes names if not none 
    if ds_celltypes_names is not None:
        ds_celltypes_names = ds_celltypes_names.split(", ")
    
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
                celltype_names = ds_celltypes_names,
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
    
    # Create PAGA integration class instance 
    integration_paga = IntegrationPAGA(
        adata = adata_concat,
        root_celltype = root_celltype
    )
    
    # Integrate across subsets (including unintegrated)
    unintegrated = integration_paga.unintegrated()
    harmony_integrated = integration_paga.harmony_integrate()
    scvi_integrated = integration_paga.scvi_integrate()
    bbknn_integrated = integration_paga.bbknn_integrate()
    scanorama_integrated = integration_paga.scanorama_integrate()
    seurat_integrated = integration_paga.seurat_integrate()
    liger_integrated = integration_paga.liger_integrate()
    
    # Add integration type to each subset and concatenate
    unintegrated.obs["integration_method"] = "unintegrated"
    harmony_integrated.obs["integration_method"] = "harmony" 
    scvi_integrated.obs["integration_method"] = "scvi"
    bbknn_integrated.obs["integration_method"] = "bbknn"
    scanorama_integrated.obs["integration_method"] = "scanorama"
    seurat_integrated.obs["integration_method"] = "seurat"
    liger_integrated.obs["integration_method"] = "liger"
    
    integrated_concat = ann.concat([
        unintegrated,
        harmony_integrated,
        scvi_integrated,
        bbknn_integrated,
        scanorama_integrated,
        seurat_integrated,
        liger_integrated
    ])
    integrated_concat.obs_names = range(len(integrated_concat.obs_names))
    integrated_concat.obs_names_make_unique()
    
    # Add data about downsampling to .uns of adata_concat
    if num_batches == 0:
        integrated_concat.uns["downsampling_stats"] = {
            "num_batches": 0,
            "num_celltypes_downsampled": ds_celltypes,
            "ds_batch_names": "None",
            "proportion_downsampled": ds_proportions,
            "downsampled_celltypes": "None"
        }
    else:
        integrated_concat.uns["downsampling_stats"] = {
            "num_batches": num_batches,
            "num_celltypes_downsampled": ds_celltypes,
            "ds_batch_names": np.concatenate([np.unique(adata.obs["batch"].__array__()) for adata in adata_downsampled]),
            "proportion_downsampled": ds_proportions,
            "downsampled_celltypes": np.array(selected_celltypes_downsampled)
        }
        
    # Save integrated h5ad object
    integrated_concat.write_h5ad(
        filename = save_loc,
        compression = "gzip"
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for scRNA-seq PAGA integration"
    )
    parser.add_argument(
        "--filedir",
        type = str,
        help = "Path of directory containing scRNA-seq h5ad files"
    )
    parser.add_argument(
        "--root_celltype",
        type = str,
        help = "Root celltype to utilize for diffusion pseudotime estimation"
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
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving output from scRNA-seq integration"
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        root_celltype  = args.root_celltype,
        save_loc = args.outfile,
        ds_celltypes = args.ds_celltypes,
        ds_celltypes_names = args.ds_celltype_names,
        ds_proportions = args.ds_proportions,
        num_batches = args.num_batches        
    )
