import numpy as np
import scanpy as sc
import anndata as ann

def cluster_num(adata):
    clusters_unique = np.unique(adata.obs.leiden.__array__())
    clusters_len = len(clusters_unique)
    return clusters_len

def leiden_clip(adata, num_clusters, step_size = 0.05):
    counter = 0
    leiden_resolution = 1
    while cluster_num(adata) != num_clusters:
        if cluster_num(adata) < num_clusters:
            leiden_resolution += step_size
            sc.tl.leiden(adata, resolution = leiden_resolution)
        elif cluster_num(adata) > num_clusters:
            leiden_resolution -= step_size
            sc.tl.leiden(adata, resolution = leiden_resolution)
        counter += 1
        if counter > 100:
            raise Exception(
                "Attempted more than 100 iterations - convergence not possible, set lower step size"
            )
    return adata

