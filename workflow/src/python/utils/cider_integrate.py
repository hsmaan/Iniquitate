import os 
import sys
import subprocess
import uuid

from natsort import natsorted
import numpy as np
import pandas as pd
import scipy as sp
import anndata as ann
import scanpy as sc

class CiderIntegrate:
    """
    Class for interpolating between the Integration class and R-script 
        based integration of RNA-seq batches using the CIDER package.
        As CIDER integration returns a set of cluster annotations, the
        anndata object has the cluster annotations appended and returned
    """
    def __init__(self, adata):
        """
        Args:
            adata (object): An instance of an anndata class corresponding to the CIDER
                subset from the Integration class.
        """
        self.adata = adata.copy()
        
    def _format(self):
        # Append a column on gene names 
        self.adata.var["gene"] = self.adata.var_names
        # Remove layers and raw from AnnData object (avoid conflicts with h5seurat)
        self.adata.layers = None
        self.adata.raw = None
        # Add 'Batch' column to obs as CIDER requires this
        self.adata.obs["Batch"] = self.adata.obs["batch"].values
        
    def _output_temp_h5ad(self):
        # Check if temp exists, if not, make dir
        if not os.path.exists("tmp"):
            os.makedirs("tmp")

        # Output temporary file with data 
        self.filename = ''.join(str(uuid.uuid4()).split("-"))
        self.file = "{filename}.h5ad".format(filename = self.filename)
        self.adata.write_h5ad(os.path.join("tmp", self.file))
    
    def _cider_integrate(self):
        # Call subprocess and call R script
        tempfile_script = \
            "Rscript src/R/cider_integrate.R tmp/{tempfile} {tempfile_name} --verbose".format(
                tempfile = self.file,
                tempfile_name = self.filename
            )
            
        self.sp_integrate = subprocess.run(tempfile_script, shell = True, text = True, capture_output = True)
        if self.sp_integrate.returncode != 0:
            raise Exception(
                "Subprocess call returned nonzero exit code - call: {call} \n Output: {output}".format(
                    call = self.sp_integrate.stderr,
                    output = self.sp_integrate.stdout
                )
            )
        
    def _return_integrated(self):
        # Get cider output file
        self.cider_outfile = "{filename}_cider_out.txt".format(filename = self.filename)

        # Read in as pandas dataframe
        integrated_cluster_annotations = pd.read_csv(os.path.join("tmp", self.cider_outfile), sep = "\t")
        
        # Append the cluster annotations to the original AnnData object - store as leiden
        cider_annotations = integrated_cluster_annotations["cider_cluster"].values
        self.adata.obs["leiden"] = pd.Categorical(
            values=cider_annotations.astype('U'),
            categories=natsorted(map(str, np.unique(cider_annotations))),
        ) 
        
        # Return integrated AnnData object
        return self.adata
    
    def _clean_files(self):
        # Remove temporary python and cider files
        tmp_files = os.listdir("tmp")
        tmp_files_instance = [f for f in tmp_files if self.filename in f]
        for f in tmp_files_instance:
            os.remove(os.path.join("tmp", f))
        
        # Check if all files related to the filename are removed from folder
        tmp_files = os.listdir("tmp")
        tmp_files_instance = [f for f in tmp_files if self.filename in f]
        if len(tmp_files_instance) > 0:
            raise Exception(
                "Temporary file cleanup incomplete - files remain in folder"
            )
            
    def integrate(self):
        # Perform workflow and return integrated anndata object
        self._format()
        self._output_temp_h5ad()
        self._cider_integrate()
        integrated_adata = self._return_integrated()
        self._clean_files()
        
        return integrated_adata