import argparse 
import os 
import sys 
sys.path.append("src/python/")

import numpy as np
import pandas as pd
import anndata as ann
import scanpy as sc
from sklearn.metrics import accuracy_score, balanced_accuracy_score, \
    f1_score, classification_report

def none_or_str(value):
    if value == 'None':
        return None
    return value

def main(h5ad_loc, save_loc, annofile, dataset_name, ds_celltypes, ds_proportions, 
         num_batches, rep):
    # Load h5ad file for query to reference mapping results
    adata = sc.read_h5ad(h5ad_loc)
    
    # Load acceptable annotations
    annos = pd.read_csv(annofile, sep="\t")
    
    # Get the classification results as a dataframe 
    class_results = pd.DataFrame({
        "Real celltype": adata.obs["celltype"],
        "Predicted L1": adata.obs["predicted.celltype.l1"],
        "Predicted L2": adata.obs["predicted.celltype.l2"],
        "Control predicted L1": adata.obs["baseline.knn.l1"],
        "Control predicted L2": adata.obs["baseline.knn.l2"]
    }) 
    
    # Format classification results by acceptable annotations - both L1 and L2,
    # control (baseline) and query to reference predictions
    class_results_l1 = []
    for celltype_real, celltype_pred in zip(
        class_results["Real celltype"], class_results["Predicted L1"]
    ):
        l1_accept = annos[annos["Real celltype"] == celltype_real]["Acceptable L1"].__array__()[0]
        if celltype_pred in l1_accept:
            class_results_l1.append(celltype_real)
        else:
            class_results_l1.append("Incorrect")
    
    class_results_l2 = []
    for celltype_real, celltype_pred in zip(
        class_results["Real celltype"], class_results["Predicted L2"]
    ):
        l2_accept = annos[annos["Real celltype"] == celltype_real]["Acceptable L2"].__array__()[0]
        if celltype_pred in l2_accept:
            class_results_l2.append(celltype_real)
        else:
            class_results_l2.append("Incorrect")
    
    baseline_results_l1 = []
    for celltype_real, celltype_pred in zip(
        class_results["Real celltype"], class_results["Control predicted L1"]
    ):
        l1_accept = annos[annos["Real celltype"] == celltype_real]["Acceptable L1"].__array__()[0]
        if celltype_pred in l1_accept:
            baseline_results_l1.append(celltype_real)
        else:
            baseline_results_l1.append("Incorrect")
    
    baseline_results_l2 = []
    for celltype_real, celltype_pred in zip(
        class_results["Real celltype"], class_results["Control predicted L2"]
    ):
        l2_accept = annos[annos["Real celltype"] == celltype_real]["Acceptable L2"].__array__()[0]
        if celltype_pred in l2_accept:
            baseline_results_l2.append(celltype_real)
        else:
            baseline_results_l2.append("Incorrect")
    
    # Append formatted results to class_results
    class_results["Predicted L1 Formatted"] = class_results_l1
    class_results["Predicted L2 Formatted"] = class_results_l2
    class_results["Baseline L1 Formatted"] = baseline_results_l1
    class_results["Baseline L2 Formatted"] = baseline_results_l2
    
    # Compute scores for predicted L1 results
    l1_true = class_results["Real celltype"].__array__()
    l1_pred = class_results["Predicted L1 Formatted"].__array__()
    l1_accuracy = accuracy_score(l1_true, l1_pred)
    l1_bal_accuracy = balanced_accuracy_score(l1_true, l1_pred)
    l1_f1 = f1_score(l1_true, l1_pred, average = "micro")
    l1_class_report = pd.DataFrame(
        classification_report(l1_true, l1_pred, output_dict = True)
    )
    l1_class_report["Score type"] = ["precision", "recall", "f1-score", "support"]
    l1_class_report = l1_class_report[["Score type"] + [col for col in l1_class_report.columns if col != "Score type"]]
    l1_class_report["Overall accuracy"] = l1_accuracy
    l1_class_report["Overall balanced accuracy"] = l1_bal_accuracy
    l1_class_report["Overall F1-score"] = l1_f1
    l1_class_report["Subset"] = "L1"
    
    # Compute scores for predicted L2 results 
    l2_true = class_results["Real celltype"].__array__()
    l2_pred = class_results["Predicted L2 Formatted"].__array__()
    l2_accuracy = accuracy_score(l2_true, l2_pred)
    l2_bal_accuracy = balanced_accuracy_score(l2_true, l2_pred)
    l2_f1 = f1_score(l2_true, l2_pred, average = "micro")
    l2_class_report = pd.DataFrame(
        classification_report(l2_true, l2_pred, output_dict = True)
    )
    l2_class_report["Score type"] = ["precision", "recall", "f1-score", "support"]
    l2_class_report = l2_class_report[["Score type"] + [col for col in l2_class_report.columns if col != "Score type"]]
    l2_class_report["Overall accuracy"] = l2_accuracy
    l2_class_report["Overall balanced accuracy"] = l2_bal_accuracy
    l2_class_report["Overall F1-score"] = l2_f1
    l2_class_report["Subset"] = "L2"
    
    # Compute scores for baseline L1 results
    l1_true = class_results["Real celltype"].__array__()
    l1_baseline_pred = class_results["Baseline L1 Formatted"].__array__()
    l1_baseline_accuracy = accuracy_score(l1_true, l1_baseline_pred)
    l1_baseline_bal_accuracy = balanced_accuracy_score(l1_true, l1_baseline_pred)
    l1_baseline_f1 = f1_score(l1_true, l1_baseline_pred, average = "micro")
    l1_baseline_class_report = pd.DataFrame(
        classification_report(l1_true, l1_baseline_pred, output_dict = True)
    )
    l1_baseline_class_report["Score type"] = ["precision", "recall", "f1-score", "support"]
    l1_baseline_class_report = l1_baseline_class_report[
        ["Score type"] + [col for col in l1_baseline_class_report.columns if col != "Score type"]
    ]
    l1_baseline_class_report["Overall accuracy"] = l1_baseline_accuracy
    l1_baseline_class_report["Overall balanced accuracy"] = l1_baseline_bal_accuracy
    l1_baseline_class_report["Overall F1-score"] = l1_baseline_f1
    l1_baseline_class_report["Subset"] = "L1 baseline"
    
    # Compute scores for baseline L2 results 
    l2_true = class_results["Real celltype"].__array__()
    l2_baseline_pred = class_results["Baseline L2 Formatted"].__array__()
    l2_baseline_accuracy = accuracy_score(l2_true, l2_baseline_pred)
    l2_baseline_bal_accuracy = balanced_accuracy_score(l2_true, l2_baseline_pred)
    l2_baseline_f1 = f1_score(l2_true, l2_baseline_pred, average = "micro")
    l2_baseline_class_report = pd.DataFrame(
        classification_report(l2_true, l2_baseline_pred, output_dict = True)
    )
    l2_baseline_class_report["Score type"] = ["precision", "recall", "f1-score", "support"]
    l2_baseline_class_report = l2_baseline_class_report[
        ["Score type"] + [col for col in l2_baseline_class_report.columns if col != "Score type"]
    ]
    l2_baseline_class_report["Overall accuracy"] = l2_baseline_accuracy
    l2_baseline_class_report["Overall balanced accuracy"] = l2_baseline_bal_accuracy
    l2_baseline_class_report["Overall F1-score"] = l2_baseline_f1
    l2_baseline_class_report["Subset"] = "L2 baseline"
    
    # Concatenate L1 and L2 results - predicted and baseline
    l1_l2_results_all = pd.concat(
        [
            l1_class_report, 
            l2_class_report, 
            l1_baseline_class_report, 
            l2_baseline_class_report
        ], 
        axis=0
    )
    
    # Append information on dataset to results
    l1_l2_results_all["Dataset"] = dataset_name
    l1_l2_results_all["Number of batches downsampled"] = num_batches
    l1_l2_results_all["Number of celltypes downsampled"] = ds_celltypes
    l1_l2_results_all["Proportion downsampled"] = ds_proportions
    l1_l2_results_all["Replicate"] = rep
    
    # Save results to file
    l1_l2_results_all.to_csv(save_loc, index=False, sep="\t")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for annotation scoring summary"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of Seurat annotated h5ad file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving annotation analysis and scoring results"
    )
    parser.add_argument(
        "--annofile",
        type = str,
        help = "Filepath for acceptable annotation results/matched with actual celltypes"
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
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        save_loc = args.outfile,
        annofile = args.annofile,
        dataset_name = args.dataset,
        rep = args.rep,
        ds_celltypes = args.ds_celltypes,
        ds_proportions = args.ds_proportions,
        num_batches = args.num_batches  
    )