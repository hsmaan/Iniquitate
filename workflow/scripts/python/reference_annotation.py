import argparse 
import sys 
sys.path.append("src/python/")

from utils import SeuratReferenceMap

def main(h5ad_loc, ref_h5_loc, save_loc):
    # Create an instance of SeuratReferenceMap
    refmap = SeuratReferenceMap(
        integrated_data_h5 = h5ad_loc, 
        reference_h5 = ref_h5_loc, 
        mapped_h5 = save_loc
    )
    
    # Run SeuratReferenceMap to save mapped query object to h5ad file
    refmap.refmap()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = "Input and output files for query to reference mapping"
    )
    parser.add_argument(
        "--infile",
        type = str,
        help = "Path of integrated h5ad file"
    )
    parser.add_argument(
        "--ref_file",
        type = str,
        help = "Path of reference h5Seurat file"
    )
    parser.add_argument(
        "--outfile",
        type = str,
        help = "Filepath for saving Seurat reference mapped and annotated h5ad file"
    )
    args = parser.parse_args()
    main(
        h5ad_loc = args.infile,
        ref_h5_loc = args.ref_file,
        save_loc = args.outfile,
    )
    