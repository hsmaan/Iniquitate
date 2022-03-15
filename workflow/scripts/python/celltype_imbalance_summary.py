import argparse 
import functools

import numpy as np
import pandas as pd
import scanpy as sc 

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata_full = sc.read_h5ad(h5ad_loc)
    
    # Extract data from just one integration method subset 
    int_method_select = np.random.choice(
        np.unique(adata_full.obs.integration_method.__array__())
    )
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata_full.uns["downsampling_stats"]["num_batches"]
    batches_ds = adata_full.uns["downsampling_stats"]["ds_batch_names"]
    num_celltypes_ds = adata_full.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata_full.uns["downsampling_stats"]["proportion_downsampled"]
    downsampled_celltypes = adata_full.uns["downsampling_stats"]["downsampled_celltypes"]
    
    # Format downsampled celltypes and batches to correspond to a single item
    if isinstance(downsampled_celltypes, str):
        if downsampled_celltypes == "None":  
            downsampled_celltypes = "None"
        else:
            raise ValueError("Downsampled celltypes is a str and not 'None'")
    elif isinstance(downsampled_celltypes, np.ndarray):
        if downsampled_celltypes.shape == (1,):
            downsampled_celltypes = downsampled_celltypes[0]
        else:
            downsampled_celltypes = np.concatenate(downsampled_celltypes).flatten()
            downsampled_celltypes = ", ".join(downsampled_celltypes)
    else:
        raise TypeError("Downsampled celltypes is not a str or ndarray")
    
    if isinstance(batches_ds, str):
        if batches_ds == "None":
            batches_ds = "None"
        else:
            raise ValueError("Downsampled batches is a str and not 'None'")
    elif isinstance(batches_ds, np.ndarray):
        if batches_ds.shape == (1,):
            batches_ds = batches_ds[0]
        else:
            batches_ds = np.concatenate(batches_ds).flatten()
            batches_ds = ", ".join(batches_ds)
    else:
        raise TypeError("Downsampled batches is not a str or ndarray")
    
    # Subset data for only one method and split datasets by batch
    adata_select = adata_full[adata_full.obs.integration_method == int_method_select]
    adata_list = []
    batches = np.unique(adata_select.obs.batch.__array__())
    for batch in batches:
        adata_batch_select = adata_select[adata_select.obs.batch == batch]
        adata_list.append(adata_batch_select)
        
    # Get celltype value counts for each batch
    val_counts_dfs = []
    for idx, adata in enumerate(adata_list):
        val_counts_df = pd.DataFrame(adata.obs.celltype.value_counts())
        val_counts_df = val_counts_df.reset_index()
        val_counts_df.columns = ["celltype", "celltype_count_batch_{}".format(idx)]
        val_counts_dfs.append(val_counts_df)
        
    # Concatenate all celltype value counts results 
    merge = functools.partial(pd.merge, on = ["celltype"], how = "outer")
    val_counts_merged = functools.reduce(merge, val_counts_dfs)
    
    # Replace NAs with 0 and add downsampling information
    val_counts_merged.iloc[:, 1:] = val_counts_merged.iloc[:, 1:].fillna(0)
    val_counts_merged["Dataset"] = dataset_name
    val_counts_merged["Number of batches downsampled"] = num_batches_ds
    val_counts_merged["Batches downsampled"] = batches_ds
    val_counts_merged["Number of celltypes downsampled"] = num_celltypes_ds
    val_counts_merged["Proportion downsampled"] = prop_ds
    val_counts_merged["Downsampled celltypes"] = downsampled_celltypes
    val_counts_merged["Replicate"] = rep
    val_counts_merged["Total batches"] = len(batches)
    val_counts_merged.to_csv(save_loc, index=False, sep="\t")
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for celltype imbalance summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving celltype imbalance statistics of h5ad file"
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