import argparse 
import os 
import sys 
sys.path.append("src/python/")
os.environ['CUDA_VISIBLE_DEVICES'] = "0, 1"

import numpy as np
import pandas as pd 

def main(leiden_file_list, save_loc):
    # Load all of the leiden files
    leiden_results = [pd.read_csv(f, sep='\t') for f in leiden_file_list]
    
    # Concatenate all of the leiden files by row
    leiden_results_concat = pd.concat(leiden_results, axis=0)

    # For each of the iterations, find the resolution that has the minimum cluster difference
    # for each method. Append the iter number and resolution number to a list
    leiden_results_sub_iters = []
    leiden_results_sub_resolutions = []
    leiden_results_sub_methods = []
    for iteration in np.unique(leiden_results_concat["iter"]):
        for method in np.unique(leiden_results_concat["method"]):
            leiden_results_sub = leiden_results_concat.loc[
                (leiden_results_concat["iter"] == iteration) & (leiden_results_concat["method"] == method)
            ]
            leiden_results_sub_iters.append(iteration)
            leiden_results_sub_resolutions.append(
                leiden_results_sub.iloc[leiden_results_sub["diff"].argmin()]["resolution"]
            )
            leiden_results_sub_methods.append(method)


    # Subset the concatenated leiden results by the iterations, resolutions, and methods in the lists
    leiden_results_best_res_sub = []
    for iteration, resolution, method in zip(
            leiden_results_sub_iters, leiden_results_sub_resolutions, leiden_results_sub_methods
        ):
        leiden_results_best_res_sub.append(
            leiden_results_concat.loc[(leiden_results_concat["iter"] == iteration) & 
                                    (leiden_results_concat["resolution"] == resolution) & 
                                    (leiden_results_concat["method"] == method)])
        
    leiden_results_best_res_concat = pd.concat(leiden_results_best_res_sub, axis=0)

    # Get the median value for each method 
    median_vals = leiden_results_best_res_concat.groupby(["method", "resolution"])["diff"].median()

    # Get the methods and their median values by iterating through multi-index
    median_vals_index = median_vals.index.values
    median_vals_list = [val[1] for val in median_vals_index]
    median_vals_methods_list = [val[0] for val in median_vals_index]

    # Create a dataframe with the median values and the methods 
    resolutions_df = pd.DataFrame({"method": median_vals_methods_list, "resolution": median_vals_list})
    
    # Order the dataframe by method based on the order indicated in 'methods_order'
    methods_order = ["scvi", "harmony", "bbknn", "scanorama", "seurat", "liger"]
    resolutions_df["method"] = pd.Categorical(
        resolutions_df["method"],
        categories = methods_order,
        ordered = True
    )
    resolutions_df = resolutions_df.sort_values("method")
        
    # Save the dataframe
    resolutions_df.to_csv(save_loc, sep='\t', index=False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for scRNA-seq integration"
    )
    parser.add_argument(
        "--files",
        type = str,
        nargs = "*",
        help = "List of leiden files to aggregate"
    )
    parser.add_argument(
        "--outfile",
        type = str,
    )
    args = parser.parse_args()
    main(
        leiden_file_list = args.files,
        save_loc = args.outfile
    )    