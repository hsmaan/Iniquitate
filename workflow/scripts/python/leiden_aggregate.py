import argparse 
import os 
import sys 
sys.path.append("src/python/")
os.environ['CUDA_VISIBLE_DEVICES'] = "0, 1"

import numpy as np
import pandas as pd 

def main(leiden_file_list, save_loc):
    # Load all of the leiden files
    leiden_results = [pd.read_csv(f, sep='\t', index_col=0) for f in leiden_file_list]
    
    # For each of the replicates, find the resolution that has the minimum cluster difference
    # for each method
    results_sub = [] 
    for leiden_result in leiden_results:
        leiden_result_sub = leiden_result.groupby(["method", "iter", "resolution"])["diff"].idxmin() 
        leiden_result_sub = leiden_result.loc[leiden_result_sub]
        results_sub.append(leiden_result_sub)
        
    # Concatenate all of the leiden files by column
    results_sub_concat = pd.concat(results_sub, axis=1)
    
    # Get the median value for each method 
    median_vals = results_sub_concat.groupby(["method", "resolution"])["diff"].median()
    
    # Convert to a list 
    median_vals_list = median_vals.to_list()
    
    # Create a dataframe with the median values and the methods 
    methods = results_sub[0]["method"].values
    resolutions_df = pd.DataFrame({"method": methods, "resolution": median_vals_list})
    
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
        leiden_file_list = args.filedir,
        save_loc = args.outfile
    )    