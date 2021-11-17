import os 
import sys
import subprocess
import uuid

import numpy as np
import pandas as pd
import scipy as sp
import anndata as ann
import scanpy as sc

class SeuratIntegrate:
    """
    Class for interpolating between the Integration class and R-script 
        based integration of RNA-seq batches using the Seurat4.0 package.
        Integration is done on data output to a temporary file from the 
        Integration class through and RScript, which then outputs a temporary
        file reread into python code and used to substitute the unintegrated data.
        Uses the first 20 components of CCA/RPCA space to integrate and get correction
        vectors for correcting whole data matrix.
    """
    def __init__(self, adata, conda_source_loc, conda_env_loc, int_type = "CCA"):
        """
        Args:
            adata (object): An instance of an anndata class corresponding to the seurat
                subset from the Integration class.
            conda_source_loc (string): A string indicating the path of the conda config
                file that has to be sourced to activate the necessary environment.
            conda_env_loc (string): A string indicating the path of conda environment 
                that has the R and Seurat dependencies needed to execute the subprocess
                script.
            int_type (string): Either "CCA" or "RPCA", indicating which Seurat workflow to 
                utilize for integration - canonical correlation analysis and reciprocal PCA,
                respectively. RPCA should be used for larger datasets to avoid out-of-memory
                exceptions. Details on both workflows can be found at:
                https://satijalab.org/seurat/articles/integration_introduction.html 
                https://satijalab.org/seurat/articles/integration_rpca.html
        """
        self.adata = adata.copy()
        self.conda_source_loc = conda_source_loc
        self.conda_env_loc = conda_env_loc
        self.int_type = int_type
        
    def _format(self):
        # Extract and format integration dataframe object 
        if isinstance(self.adata.X, sp.sparse.csr.csr_matrix):
            self.adata.X = self.adata.X.toarray()
        integration_mat = self.adata.X
        self.integration_df = pd.DataFrame(integration_mat)
        gene_names = self.adata.var["gene"]
        self.integration_df.columns = gene_names
        self.integration_df.columns = self.integration_df.columns.astype(str) # Clip columns to str not cat
        self.integration_df.index = [str(i) + "_bc" for i in range(len(self.integration_df))]
        self.integration_df["batch"] = self.adata.obs["batch_name"].values
        
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
    
    def _seurat_integrate(self):
        # Call subprocess, activate conda env, and call R script
        tempfile_script = \
            "bash -c 'source {conda_source}".format(conda_source = self.conda_source_loc) + \
            " && conda activate {conda_env}".format(conda_env = self.conda_env_loc) + \
            " && Rscript R/seurat.R tmp/{tempfile} {tempfile_name} {env_loc} {int_type} --verbose".format(
                tempfile = self.file,
                tempfile_name = self.filename,
                env_loc = self.conda_env_loc,
                int_type = self.int_type
            ) + \
            " && conda deactivate'"
            
        self.sp_integrate = subprocess.run(tempfile_script, shell = True, text = True, capture_output = True)
        if self.sp_integrate.returncode != 0:
            raise Exception(
                "Subprocess call returned nonzero exit code - call: {call}".format(
                    call = self.sp_integrate.stderr
                )
            )
        
    def _return_integrated(self):
        # Get seurat output file
        self.seur_outfile = "{filename}_seur_out.h5ad".format(filename = self.filename)

        # Read in as AnnData object 
        integrated_adata = sc.read_h5ad(filename = os.path.join("tmp", self.seur_outfile))
        
        # Return integrated AnnData object
        return integrated_adata
    
    def _clean_files(self):
        # Remove temporary python and seurat files
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
        self._seurat_integrate()
        integrated_adata = self._return_integrated()
        self._clean_files()
        
        return integrated_adata