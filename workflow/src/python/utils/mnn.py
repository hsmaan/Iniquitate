from itertools import combinations

import faiss
import numpy as np

def find_mutual_nn(data_list, k = 15):
    """Gets mutual nearest neighbors pairs across all datasets.

    Using each data subset in data_list, gets mutual nearest neighbors by
    considering the intersection of MNNs across all datasets.

    Args:
        data_list (list): List of numpy arrays corresponding to data subsets. Datasets
            must have an equivalent number of features.
        k (integer): Positive integer value indicating how many neighbors to consider
            in the mutual-nearest-neighbors algorithm. Default value is 15.

    Returns:
        mnn_1_concat (array): Array of values corresponding to query MNN indices (MNN_1) that can be in
            any dataset, and are indexed based on the concatenated representation of all datasets in
            dataset_list.
        mnn_2_concat (array): Array of values corresponding to query-value MNN indices (MNN_2) that can
            be in any dataset, and are indexed based on the concatenated representation of all datasets in
            dataset_list.
    """
    # Get lengths of all datasets for reindexing
    data_lens = [len(dataset) for dataset in data_list]

    # Get all indices
    indices = [i for i in range(len(data_list))]

    # Create combinations for all indices
    index_combo_iter = combinations(indices, 2)
    index_combos = [combo for combo in index_combo_iter]

    # Iterate over combinations and record mnn pairs - append to both lists
    mnn_1_list = []
    mnn_2_list = []
    for combo in index_combos:
        idx_1 = combo[0]
        idx_2 = combo[1]
        data_1 = data_list[idx_1]
        data_2 = data_list[idx_2]
        mnn_1, mnn_2 = mutual_nn(data_1, data_2, k1=k, k2=k)
        if idx_1 == 0:
            pass
        else:
            len_addition = sum(data_lens[0:idx_1])
            mnn_1 = mnn_1 + len_addition
        if idx_2 == 0:
            pass
        else:
            len_addition = sum(data_lens[0:idx_2])
            mnn_2 = mnn_2 + len_addition
        mnn_concat_1 = np.concatenate((mnn_1, mnn_2))
        mnn_concat_2 = np.concatenate((mnn_2, mnn_1))
        mnn_1_list.append(mnn_concat_1)
        mnn_2_list.append(mnn_concat_2)

    # Concatenate MNNs in mnn_1 and mnn_2 into one array and return
    mnn_1_concat = np.concatenate(mnn_1_list)
    mnn_2_concat = np.concatenate(mnn_2_list)
    return mnn_1_concat, mnn_2_concat


def mutual_nn(data_1, data_2, k1, k2):
    """Given two datasets, gets and returns mutual nearest neighbors.

    Args:
        data_1 (array): Data array 1 that is used to create the graph representing
            dataset 1. Dataset 1 and 2 must have the same numbers of featres.
        data_2 (array): Data array 1 that is used to create the graph representing
            dataset 2. Dataset 1 and 2 must have the same numbers of features.
        k1 (integer): Positive integer value indicating how many neighbors to consider
            in the mutual-nearest-neighbors algorithm for dataset 1.
        k2 (integer): Positive integer value indicating how many neighbors to consider
            in the mutual-nearest-neighbors algorithm for dataset 2.

    Returns:
        mutual_1_arr (array): Array of mutual-nearest neighbors corresponding to indices in dataset 1.
        mutual_2_arr (array): Array of mutual-nearest neighbors corresponding to indices in dataset 2.
    """

    data_1 = np.ascontiguousarray(data_1, dtype = np.float32)
    data_2 = np.ascontiguousarray(data_2, dtype = np.float32)
    
    index_1 = faiss.IndexFlatL2(data_1.shape[1])
    index_2 = faiss.IndexFlatL2(data_2.shape[1])
    
    index_1.add(data_1)
    index_2.add(data_2)
    
    d_index_1, k_index_1 = index_1.search(data_2, k1)
    d_index_2, k_index_2 = index_2.search(data_1, k2)
    
    mutual_1 = []
    mutual_2 = []
    for index_2 in range(data_2.shape[0]):
        for index_1 in k_index_1[index_2]:
            if index_2 in k_index_2[index_1]:
                mutual_1.append(index_1)
                mutual_2.append(index_2)
    mutual_1_arr = np.asarray(mutual_1)
    mutual_2_arr = np.asarray(mutual_2)
    return mutual_1_arr, mutual_2_arr

def cross_data_knn(data_1, data_2, k):
    """Given two datasets, gets and returns KNN of dataset 1 in dataset 2.
    
    data_1 (array): Data array 1 that is used to create the graph representing
        dataset 1. Dataset 1 and 2 must have the same numbers of featres.
    data_2 (array): Data array 1 that is used to create the graph representing
        dataset 2. Dataset 1 and 2 must have the same numbers of features.
    k (integer): Positive integer value indicating how many neighbors to consider
        in the cross-data knn lookup.

    Returns:
        knn_arr (array): Array of k-nearest neighbors corresponding to indices in dataset 1, 
            with indices in dataset 2.
    """
    data_1 = np.ascontiguousarray(data_1, dtype = np.float32)
    data_2 = np.ascontiguousarray(data_2, dtype = np.float32)
    
    index_2 = faiss.IndexFlatL2(data_2.shape[1])
    d_index_2, k_index_2 = index_2.search(data_1, k)
    
    return k_index_2


def find_knn(data_list, k = 15):
    """Gets k nearest-neighbors for all datasets.

    Using each data subset in data_list, gets k-nearest-neighbors for each subset and
    returns concatenated k-nearest-neighbors for each index across all datasets.

    Args:
        data_list (list): List of numpy arrays corresponding to data subsets. Datasets
            must have an equivalent number of features.
        k (integer): Positive integer value indicating how many neighbors to consider
            in the k-nearest-neighbors algorithm. Default value is 15.

    Returns:
        knn_concat (array): Array of values corresponding to query KNN indices for all
            datasets in data_list, indexed based on their concatenation. Will return
            k nearest-neighbors at query position for each index from input data.
    """
    # Get lengths of all datasets for reindexing
    data_lens = [len(dataset) for dataset in data_list]

    # Get all indices
    indices = [i for i in range(len(data_list))]

    # Get knn pairs for each dataset and reindex as necessary
    knn_list = []
    for idx in indices:
        dataset = data_list[idx]
        dataset = np.ascontiguousarray(dataset, dtype = np.float32)
        index = faiss.IndexFlatL2(dataset.shape[1])
        index.add(dataset)
        knn_vals, knn = index.search(dataset, k)
        if idx == 0:
            knn_corrected = knn
            knn_list.append(knn_corrected)
        else:
            knn_corrected = []
            len_addition = sum(data_lens[0:idx])
            for i in range(len(knn)):
                knn_corrected.append(knn[i] + len_addition)
            knn_corrected_arr = np.asarray(knn_corrected)
            knn_list.append(knn_corrected_arr)

    # Concatenate all KNNs corresponding to all dataset indices and return
    knn_concat = np.concatenate(knn_list)
    return knn_concat