import argparse 

import pandas as pd 

def main(dge_rank_file_loc, marker_file_loc, save_loc):
    # Read in the dge rank file and marker file
    dge_rank_df = pd.read_csv(dge_rank_file_loc, sep = "\t")
    marker_df = pd.read_csv(marker_file_loc, sep = "\t")
    
    # Subset the dge rank df by the markers in the marker df
    dataset_markers = marker_df["Top 10 marker genes (union across batches)"].__array__()
    dge_rank_df_marker_sub = dge_rank_df[dge_rank_df["Gene"].isin(dataset_markers)]
    
    # Save the marker subset dge rank df
    dge_rank_df_marker_sub.to_csv(save_loc, sep = "\t", index = False)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "Input and output files for dge concordance summary"
    )
    parser.add_argument(
        "--infile_dge_rank",
        type = str,
        help = "Path of dge rank file for given dataset"
    )
    parser.add_argument(
        "--infile_marker",
        type = str,
        help = "Path of marker gene file for given dataset"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving marker gene subset dge rank file"
    )
    args = parser.parse_args()
    main(
        dge_rank_file_loc=args.infile_dge_rank,
        marker_file_loc=args.infile_marker,
        save_loc=args.outfile
    )