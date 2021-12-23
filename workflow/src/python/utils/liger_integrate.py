import os 
import sys
import subprocess
import uuid

import numpy as np
import pandas as pd
import scipy as sp
import anndata as ann
import scanpy as sc

class LigerIntegrate:
    """
    Class for interpolating between the Integration class and R-script 
        based integration of RNA-seq batches using the LIGER R package.
        Integration is done on data output to a temporary file from the 
        Integration class through and RScript, which then outputs a temporary
        file reread into python code and used to substitute the unintegrated data.
        Integration is performed using 20 latent components in the NMF factorization,
        or 20 "metagenes".
    """
    def __init__(self, adata):
        """
        Args:
            adata (object): An instance of an anndata class corresponding to the liger
                subset from the Integration class.
        """
        self.adata = adata.copy()
        self.adata_copy = adata.copy() # Keep copy for later referencing
        
    def _format(self):
        # Append a column on gene names 
        self.adata.var["gene"] = self.adata.var_names
        # Remove layers and raw from AnnData object (avoid conflicts with h5seurat)
        self.adata.layers = None
        self.adata.raw = None
        
    def _output_temp_h5ad(self):
        # Check if temp exists, if not, make dir
        if not os.path.exists("tmp"):
            os.makedirs("tmp")

        # Output temporary file with data 
        self.filename = ''.join(str(uuid.uuid4()).split("-"))
        self.file = "{filename}.h5ad".format(filename = self.filename)
        self.adata.write_h5ad(os.path.join("tmp", self.file))
    
    def _liger_integrate(self):
        # Call subprocess and call R script
        tempfile_script = \
            "Rscript src/R/liger.R tmp/{tempfile} {tempfile_name} --verbose".format(
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
        # Get liger output file and read it into memory as anndata object
        self.liger_outfile = "{filename}_liger_out.h5ad".format(filename = self.filename)
        adata_liger = sc.read_h5ad(
            os.path.join("tmp", self.liger_outfile)
        )
        
        # Read in cell-specific loadings from anndata object and convert to array
        norm_loadings_arr = adata_liger.X.toarray()
        
        # Add normalized loadings to anndata object (original copy)
        self.adata_copy.obsm["X_liger"] = norm_loadings_arr
        
        # Return integrated AnnData object
        return self.adata_copy
    
    def _clean_files(self):
        # Remove temporary python and liger files
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
        self._liger_integrate()
        integrated_adata = self._return_integrated()
        self._clean_files()
        
        return integrated_adata