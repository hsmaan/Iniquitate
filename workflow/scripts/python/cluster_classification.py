import argparse 

import scanpy as sc
import anndata as ann
import numpy as np
import pandas as pd 
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
from sklearn.linear_model import LogisticRegression

def main(h5ad_loc, save_loc, dataset_name, rep):
    # Load h5ad file 
    adata = sc.read_h5ad(h5ad_loc)
    
    # Extract summary statistics from h5ad file
    num_batches_ds = adata.uns["downsampling_stats"]["num_batches"]
    num_celltypes_ds = adata.uns["downsampling_stats"]["num_celltypes_downsampled"]
    prop_ds = adata.uns["downsampling_stats"]["proportion_downsampled"]
    
    # Subset h5ad based on batch-correction method used
    adata_method_sub = []
    methods = ["harmony", "scvi", "scanorama", "bbknn", "seurat", "liger", "cider"] 
    for method in methods:
        adata_sub = adata[adata.obs["integration_method"] == method]
        adata_method_sub.append(
            adata_sub
        )

    # Determine cluster-based accuracy for each batch-correction method
    precision_scores = []
    recall_scores = []
    f1_scores = []
    supports = []
    celltypes = []
    for adata_sub in adata_method_sub:
        # Split testing and training data in stratified manner (70/30)
        X = adata_sub.obs["leiden"].__array__()
        y = adata_sub.obs["celltype"].__array__()
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, stratify=y, test_size=0.7, random_state=42
        )
        
        # Train a logistic regression classifier to predict celltype based on cluster
        X_train_num = np.array(
            [int(x) for x in X_train]
        )
        X_test_num = np.array(
            [int(x) for x in X_test]
        )
        # Reshape X test and train to be 2D arrays
        X_train_num = X_train_num.reshape(-1, 1)
        X_test_num = X_test_num.reshape(-1, 1)
        
        # Train logistic regression classifier and predict celltypes
        clf = LogisticRegression(
            random_state=42, max_iter=1000
        ).fit(X_train_num, y_train)
        y_pred = clf.predict(X_test_num)
        
        # Get classification report and subset for only relevant columns 
        class_report_dict = classification_report(
            y_test, y_pred, output_dict=True
        )
        class_report_df = pd.DataFrame(class_report_dict)
        class_report_df = class_report_df.iloc[:, :-3]
        
        # Append appropriate values to lists
        precision_scores.append(class_report_df.loc["precision"].values)
        recall_scores.append(class_report_df.loc["recall"].values)
        f1_scores.append(class_report_df.loc["f1-score"].values)
        supports.append(class_report_df.loc["support"].values)
        celltypes.append(class_report_df.columns.values)
        
    # Repeat method values to have same length as scores (one for each celltype)
    methods_repeat = np.repeat(methods, len(precision_scores[0]))
    
    # Concatenate scores and celltypes 
    precision_scores_concat = np.concatenate(precision_scores)
    recall_scores_concat = np.concatenate(recall_scores)
    f1_scores_concat = np.concatenate(f1_scores)
    supports_concat = np.concatenate(supports)
    celltypes_concat = np.concatenate(celltypes)
    
    # Create summary dataframe for classification statistics and save
    classification_summary_df = pd.DataFrame({
        "Dataset": dataset_name,
        "Number of batches downsampled": num_batches_ds,
        "Number of celltypes downsampled": num_celltypes_ds,
        "Proportion downsampled": prop_ds,
        "Replicate": rep,
        "Method": methods_repeat,
        "Celltype": celltypes_concat,
        "Precision": precision_scores_concat,
        "Recall": recall_scores_concat,
        "F1-score": f1_scores_concat,
        "Support": supports_concat,
        "Mean LogReg F1-score": np.mean(f1_scores_concat)
    })
    classification_summary_df.to_csv(
        save_loc,
        index=False,
        sep="\t"
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for cluster classification summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving cluster classification results of integrated h5ad file"
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