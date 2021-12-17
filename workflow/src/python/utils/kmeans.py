import numpy as np
import faiss 

# Perform k-means clustering using Faiss
def faiss_kmeans(adata, k, niter = 300, nredo = 10, 
                 min_points_per_centroid = 5):
    """Function to perform k-means clustering using FAISS on AnnData objects.

    Args:
        adata (AnnData): An object of AnnData class with highly variable gene selection
            performed.
        k (int): Number of clusters to form.
        niter (int): Number of iterations to run k-means. Defaults to 300.
        nredo (int): Number of times to run k-means - selects best result. 
            Defaults to 10.
        min_points_per_centroid (int): Minimum number of points per k-means
            centroid. Defaults to 5.
    """
    # Subset data to kmeans reduction to utilize for clustering
    reduction_sub = adata.obsm["X_kmeans"]
    
    # Run k-means using faiss given options 
    kmeans_faiss = faiss.Kmeans(
        d = reduction_sub.shape[1], 
        k = k, 
        niter = niter, 
        nredo = nredo, 
        min_points_per_centroid = min_points_per_centroid
    )
    kmeans_faiss.train(np.ascontiguousarray(reduction_sub, dtype = np.float32))
    kmeans_faiss_labels = np.concatenate(
        kmeans_faiss.index.search(
            np.ascontiguousarray(reduction_sub, dtype = np.float32), 1
        )[1]
    )
    
    # Check if any clusters have less than the min required members and redo clustering with less
    unique_labels, counts = np.unique(kmeans_faiss_labels, return_counts = True)
    while any(counts < min_points_per_centroid):
        k -= 1
        kmeans_faiss = faiss.Kmeans(
            d = reduction_sub.shape[1], 
            k = k, 
            niter = niter, 
            nredo = nredo, 
            min_points_per_centroid = min_points_per_centroid
        )
        kmeans_faiss.train(np.ascontiguousarray(reduction_sub, dtype = np.float32))
        kmeans_faiss_labels = np.concatenate(
            kmeans_faiss.index.search(
                np.ascontiguousarray(reduction_sub, dtype = np.float32), 1
            )[1]
        )
        unique_labels, counts = np.unique(kmeans_faiss_labels, return_counts = True)
            
    # Append kmeans labels to AnnData object
    kmeans_faiss_labels_str = kmeans_faiss_labels.astype("str")
    adata.obs['kmeans_faiss'] = kmeans_faiss_labels_str
    
    # Return AnnData object and kmeans number
    return adata, k