import argparse 
import os 
import sys 
sys.path.append("src/python/")
os.environ['CUDA_VISIBLE_DEVICES'] = "0, 1"

import scanpy as sc
import anndata as ann
import numpy as np

from utils import downsample

def none_or_str(value):
    if value == 'None':
        return None
    return value

def main(h5ad_dir, save_loc, ds_celltypes, ds_proportions, num_batches):
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
    
    # Downsample first batch from loaded h5ad files based on params 
    rng = np.random.default_rng()
    adata_selected = adata_loaded[0]
    
    # Downsample the selected celltypes in the first batch
    if num_batches == 0:
        selected_celltypes_downsampled = "None" # Placeholder - not used
        batches_ds = "None" # Placeholder - not used
        adata_ds = adata_selected
    else:
        celltypes_all = np.unique(np.concatenate([adata.obs["celltype"].__array__() for adata in adata_selected]))
        rng.shuffle(celltypes_all)
        celltypes_selected = rng.choice(celltypes_all, ds_celltypes, replace = False)
        selected_celltypes_downsampled = np.array(celltypes_selected)
        adata_ds, selected_celltypes_ds = downsample(
            adata = adata_selected, 
            num_celltypes = None,
            celltype_names = celltypes_selected,
            proportion = ds_proportions
        )
    
    # Perform downstream steps to get a reduced representation of the data followed by 
    # clustering 
    sc.pp.normalize_total(
        adata_ds,
        target_sum = 1e4
    )
    sc.pp.log1p(adata_ds)
    sc.pp.highly_variable_genes(
        adata_ds,
        n_top_genes = 2500,
        flavor = "seurat"
    )
    sc.pp.pca(adata_ds, svd_solver = "arpack")
    sc.pp.neighbors(
        adata_ds,
        n_neighbors = 15,
        n_pcs = 20,
        use_rep = "X_pca"
    )
    sc.tl.leiden(adata_ds)
    sc.tl.umap(adata_ds)
    
    # If downsampled celltypes are of array length greater than one, combine them 
    if len(selected_celltypes_downsampled) > 1:
        selected_celltypes_downsampled = np.array(",".join(selected_celltypes_downsampled))

    # Add data about downsampling to .uns of selected batch
    if num_batches == 0:
        adata_ds.uns["downsampling_stats"] = {
            "num_batches": 0,
            "num_celltypes_downsampled": ds_celltypes,
            "ds_batch_names": "None",
            "proportion_downsampled": ds_proportions,
            "downsampled_celltypes": "None"
        }
    else:
        adata_ds.uns["downsampling_stats"] = {
            "num_batches": num_batches,
            "num_celltypes_downsampled": ds_celltypes,
            "ds_batch_names": "Placeholder due to h5py bug",
            "proportion_downsampled": ds_proportions,
            "downsampled_celltypes": selected_celltypes_downsampled
        }
        
    # Save downsampled and processed h5ad object
    adata_ds.write_h5ad(
        filename = save_loc,
        compression = "gzip"
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for scRNA-seq control downsampling"
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
