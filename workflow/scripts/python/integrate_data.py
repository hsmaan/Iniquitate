import argparse 
import os 
import sys 
sys.path.append("src/python/")

import pandas as pd 

from utils import Integration

def main(h5ad_dir, save_loc):
    print("Works")
    test = pd.DataFrame({
        "A": [1, 2],
        "B": [2, 1]
    })
    test.to_csv(
        save_loc
    )
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for scRNA-seq integration"
    )
    parser.add_argument(
        "--filedir",
        type = str,
        help = "Path of directory containing scRNA-seq h5ad files"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving output from scRNA-seq integration"
    )
    args = parser.parse_args()
    main(
        h5ad_dir = args.filedir,
        save_loc = args.outfile
    )