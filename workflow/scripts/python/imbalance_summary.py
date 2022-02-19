import argparse 

import numpy as np
import pandas as pd
import scipy.spatial as sp
import scanpy as sc 

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata_full = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata_full.uns["downsampling_stats"]["num_batches"]
    batches_ds = adata_full.uns["downsampling_stats"]["ds_batch_names"]
    num_celltypes_ds = adata_full.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata_full.uns["downsampling_stats"]["proportion_downsampled"]
    downsampled_celltypes = adata_full.uns["downsampling_stats"]["downsampled_celltypes"]
    
    # Concatenate downsampled celltypes and batches - conditioned on return type
    if isinstance(downsampled_celltypes, str):
        if downsampled_celltypes == "None":  
            downsampled_celltypes = "None"
        else:
            raise ValueError("Downsampled celltypes is a str and not 'None'")
    elif isinstance(downsampled_celltypes, list):
        downsampled_celltypes = np.array(downsampled_celltypes)
    elif isinstance(downsampled_celltypes, np.ndarray):
        if len(downsampled_celltypes.shape) == 1:
            downsampled_celltypes = downsampled_celltypes
        else:
            downsampled_celltypes = np.concatenate(downsampled_celltypes)
    else:
        raise TypeError("Downsampled celltypes is not a str, list, or ndarray")
    
    if isinstance(batches_ds, str):
        if batches_ds == "None":
            batches_ds = "None"
        else:
            raise ValueError("Downsampled batches is a str and not 'None'")
    elif isinstance(batches_ds, list):
        batches_ds = np.array(batches_ds)
    elif isinstance(batches_ds, np.ndarray):
        if len(batches_ds.shape) == 1:
            batches_ds = batches_ds
        else:
            batches_ds = np.concatenate(batches_ds)
    else:
        raise TypeError("Downsampled batches is not a str, list, or ndarray")
    
    # Repeat batches downsampled for each celltype downsampled
    batch_label = np.repeat(batches_ds, num_celltypes_ds)

    # Extract data from just one integration method subset 
    int_method_select = np.random.choice(
        np.unique(adata_full.obs.integration_method.__array__())
    )
    
    # Subset data for only one method and split datasets by batch
    adata_select = adata_full[adata_full.obs.integration_method == int_method_select]
    adata_list = []
    batches = np.unique(adata_select.obs.batch.__array__())
    for batch in batches:
        adata_batch_select = adata_select[adata_select.obs.batch == batch]
        adata_list.append(adata_batch_select)
    
    # Get union of cell types across all batches
    celltype_union = np.unique(np.concatenate([adata.obs.celltype.__array__() for adata in adata_list]))
    
    # Get intersection of cell types across all batches
    celltype_intersection = set.intersection(*[set(adata.obs.celltype.__array__()) for adata in adata_list])
    
    # Get proportion vector of cells in each batch
    celltype_props = []
    for adata in adata_list:
        celltype_prop = np.zeros(len(celltype_union))
        for idx, celltype in enumerate(celltype_union):
            celltype_prop[idx] = np.sum(adata.obs.celltype.__array__() == celltype)/len(adata)
        celltype_props.append(celltype_prop)
        
    # Get cosine distances across celltype proportions
    cos_distances = []
    for celltype_prop in celltype_props:
        for celltype_prop_other in celltype_props:
            cos_distances.append(sp.distance.cosine(celltype_prop, celltype_prop_other))
    
    # Get mean cosine distance across batches
    cos_dist_mean = np.mean(cos_distances)
    
    # Get ratio of unique celltypes over intersection (Jaccard index)
    celltype_unique_ratio = len(celltype_intersection) / len(celltype_union)
    
    # Get adata lens stdev proportional to total cells 
    adata_lens = [len(adata) for adata in adata_list]
    adata_lens_stdev = np.std(adata_lens)
    adata_lens_mean = np.mean(adata_lens)
    adata_coeff_var = adata_lens_stdev / adata_lens_mean
    
    # Return dataset imbalance summary stats 
    imba_summary_df = pd.DataFrame(
        {
            "Dataset": dataset_name,
            "Number of batches downsampled": num_batches_ds,
            "Batches downsampled": batch_label,
            "Number of celltypes downsampled": num_celltypes_ds,
            "Proportion downsampled": prop_ds,
            "Downsampled celltypes": downsampled_celltypes,
            "Replicate": rep,
            "Total batches": len(batches),
            "Celltype intersection ratio": celltype_unique_ratio,
            "Mean proportion cosine distance": cos_dist_mean,
            "Length coeff var": adata_coeff_var
        },
        index=[0]
    )
    imba_summary_df.to_csv(save_loc, index=False, sep="\t")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for imbalance summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving imbalance statistics of h5ad file"
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
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        save_loc = args.outfile,
        dataset_name = args.dataset,
        rep = args.rep
    )