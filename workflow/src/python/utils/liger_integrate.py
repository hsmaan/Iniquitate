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
        # Extract and format integration dataframe object 
        if isinstance(self.adata.X, sp.sparse.csr.csr_matrix):
            self.adata.X = self.adata.X.toarray()
        integration_mat = self.adata.X
        self.integration_df = pd.DataFrame(integration_mat)
        gene_names = self.adata.var["gene"].__array__()
        self.integration_df.columns = gene_names
        self.integration_df.columns = self.integration_df.columns.astype(str) # Clip columns to str not cat
        self.integration_df.index = [str(i) + "_bc" for i in range(len(self.integration_df))]
        self.integration_df["batch"] = self.adata.obs["batch"].__array__()
        
    def _output_temp_df(self):
        # Check if temp exists, if not, make dir
        if not os.path.exists("tmp"):
            os.makedirs("tmp")

        # Output temporary file with data 
        self.filename = ''.join(str(uuid.uuid4()).split("-"))
        self.file = "{filename}.tsv".format(filename = self.filename)
        self.integration_df.to_csv(
            os.path.join("tmp", self.file),
            sep = "\t",
            index = True
        )
    
    def _liger_integrate(self):
        # Call subprocess and call R script
        tempfile_script = \
            "Rscript R/liger.R tmp/{tempfile} {tempfile_name} --verbose".format(
                tempfile = self.file,
                tempfile_name = self.filename
            )
            
        self.sp_integrate = subprocess.run(tempfile_script, shell = True, text = True, capture_output = True)
        if self.sp_integrate.returncode != 0:
            raise Exception(
                "Subprocess call returned nonzero exit code - call: {call}".format(
                    call = self.sp_integrate.stderr
                )
            )
        
    def _return_integrated(self):
        # Get liger output file
        self.liger_outfile = "{filename}_liger_out.tsv".format(filename = self.filename)

        # Read in cell-specific loadings as dataframe and convert to array
        norm_loadings_df = pd.read_csv(
            os.path.join("tmp", self.liger_outfile),
            sep = "\t"
        )
        norm_loadings_arr = norm_loadings_df.to_numpy()
        
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
        self._output_temp_df()
        self._liger_integrate()
        integrated_adata = self._return_integrated()
        self._clean_files()
        
        return integrated_adata