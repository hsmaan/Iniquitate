import argparse 
import sys 
sys.path.append("src/python/")

import numpy as np
import scanpy as sc 
import anndata as ann

from utils import cross_data_knn

def main(h5ad_loc, ref_h5_loc, save_loc):
    # Load the seurat reference mapped h5ad file 
    query_h5ad_full = sc.read_h5ad(h5ad_loc)
    query_h5ad = query_h5ad_full.raw.to_adata() # Use sctransformed counts
    
    # Load the reference h5ad file
    ref_h5ad = sc.read_h5ad(ref_h5_loc)
    
    # Get the intersection of the genes in the query and reference h5ad files
    query_genes = set(query_h5ad.var._index.__array__())
    ref_genes = set(ref_h5ad.var._index.__array__())
    common_genes_list = list(ref_genes.intersection(query_genes))

    # Get indices of the common genes in both query and ref h5ad files
    query_sct_gene_indices = np.where(
        np.isin(
            query_h5ad.var._index.__array__(), 
            common_genes_list
        )
    )[0]
    ref_sct_gene_indices = np.where(
        np.isin(
            ref_h5ad.var._index.__array__(),
            common_genes_list
        )
    )[0]

    # Subset anndata objects for the common genes
    query_h5ad_sub = query_h5ad[:, query_sct_gene_indices]
    ref_h5ad_sub = ref_h5ad[:, ref_sct_gene_indices]

    # Ensure genes are equal between query and reference 
    if not np.array_equal(
        query_h5ad_sub.var._index.__array__(), 
        ref_h5ad_sub.var._index.__array__()
    ):
        raise ValueError(
            "Genes not equal between query and reference h5ad files after intersection"
        )
    
    # Get highly variable gene subsets of both jointly and return indices 
    query_ref_sub_concat = ann.AnnData.concatenate(query_h5ad_sub, ref_h5ad_sub)
    sc.pp.highly_variable_genes(query_ref_sub_concat, flavor="seurat", n_top_genes=2500)
    hvg_indices = query_ref_sub_concat.var["highly_variable"].__array__()

    # Get the SCTransformed data subsets for both the query and reference (common gene subsets)
    query_sct = query_h5ad_sub.X.toarray()
    ref_sct = ref_h5ad_sub.X.toarray()

    # Subset for the indices of the hvg genes
    query_sct_hvg_subset = query_sct[:, hvg_indices]
    ref_sct_hvg_subset = ref_sct[:, hvg_indices]
    
    # Get the (1) nearest neighbors for the reference data within the query data
    query_1_nn = cross_data_knn(query_sct_hvg_subset, ref_sct_hvg_subset, 1)
    
    # Get the celltypes (both l1 and l2) corresponding to the nearest neighbors for the reference data
    ref_celltypes_l1 = ref_h5ad.obs["celltype.l1"][query_1_nn.flatten()].__array__()
    ref_celltypes_l2 = ref_h5ad.obs["celltype.l2"][query_1_nn.flatten()].__array__()
    
    # Append the celltypes to the query h5ad file
    query_h5ad.obs["baseline.knn.l1"] = ref_celltypes_l1
    query_h5ad.obs["baseline.knn.l2"] = ref_celltypes_l2

    # Change colnames of query var to not collide with h5ad writing in anndata
    query_h5ad.var["gene_name"] = query_h5ad.var._index 
    query_h5ad.var = query_h5ad.var.drop(query_h5ad.var.columns[0], axis=1)
    
    # Save the query h5ad file with baseline annotations
    query_h5ad.write_h5ad(
        filename = save_loc,
        compression = "gzip"
    )
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for query to reference mapping"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of seurat reference mapped h5ad file"
    )
    parser.add_argument(
        "--ref_file",
        type = str,
        help = "Path of reference h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving seurat and control reference mapped and annotated h5ad file"
    )
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        ref_h5_loc = args.ref_file,
        save_loc = args.outfile,
    )
    