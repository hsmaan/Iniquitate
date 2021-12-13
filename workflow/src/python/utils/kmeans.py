import numpy as np
import faiss 

# Perform k-means clustering using Faiss
def faiss_kmeans(adata, k, niter = 100, nredo = 10, 
                 min_points_per_centroid = 10):
    """Function to perform k-means clustering using FAISS on AnnData objects.

    Args:
        adata (AnnData): An object of AnnData class with highly variable gene selection
            performed.
        k (int): Number of clusters to form.
        niter (int): Number of iterations to run k-means. Defaults to 100
        nredo (int): Number of times to run k-means - selects best result. 
            Defaults to 10.
        min_points_per_centroid (int): Minimum number of points per k-means
            centroid. Defaults to 10.
    """
    # Subset data to highly variable genes
    hvg_sub = adata.X.toarray()[:, adata.var['highly_variable']]
    
    # Run k-means using faiss given options 
    kmeans_faiss = faiss.Kmeans(
        d = hvg_sub.shape[1], 
        k = k, 
        niter = niter, 
        nredo = nredo, 
        min_points_per_centroid = min_points_per_centroid
    )
    kmeans_faiss.train(np.ascontiguousarray(hvg_sub, dtype = np.float32))
    kmeans_faiss_labels = np.concatenate(
        kmeans_faiss.index.search(
            np.ascontiguousarray(hvg_sub, dtype = np.float32), 1
        )[1]
    )
    
    # Append kmeans labels to AnnData object
    kmeans_faiss_labels_str = kmeans_faiss_labels.astype("str")
    adata.obs['kmeans_faiss'] = kmeans_faiss_labels_str
    
    return adata