import os 
import sys
import subprocess
import uuid

import numpy as np
import pandas as pd
import scipy as sp
import anndata as ann
import scanpy as sc

class SeuratReferenceMap:
    """
    Class for reference to query mapping through integration of RNA-seq batches 
        using the Seurat4.0 package. Integration is done on data output to a temporary 
        file a downsampled and integrated result through an RScript, which then outputs the 
        reference mapped anndata (h5ad) file to be used for later downstream testing and 
        analysis. The reference mapping workflow follows that of: 
            https://satijalab.org/seurat/articles/multimodal_reference_mapping.html.
    """
    def __init__(self, integrated_data_h5, reference_h5, mapped_h5):
        """
        Args:
            integrated_data_h5 (str): Path to the anndata file for the integrated result after 
                downsampling and integration.
            reference_h5 (str): Path to the reference h5Seurat file that contains the data
                to be used in query to reference mapping and annotation.
            mapped_h5 (str): Path to the h5ad output file from seurat that contains the 
                mapped and annotated data.
        """
        self.integrated_data_h5 = integrated_data_h5
        self.reference_h5 = reference_h5
        self.mapped_h5 = mapped_h5
    
    def _load(self):
        # Load the integrated data and subset for the seurat results 
        self.adata = sc.read_h5ad(self.integrated_data_h5)
        self.adata = self.adata[self.adata.obs.integration_method == "seurat"]
        self.adata.obs.index = range(len(self.adata.obs)) # Reset index
    
    def _format(self):
        # Substitute in raw counts for X, remove unecessary obs, obsm, and uns info
        self.adata.X = self.adata.layers["raw"]
        self.adata.obs = self.adata.obs.drop(columns = [
            "leiden", 
            "integration_method",
            "kmeans_faiss"
        ])
        self.adata.obsm = None 
        del self.adata.uns # Uns doesn't support None for resetting 
        
        # Append a column on gene names 
        self.adata.var["gene"] = self.adata.var_names
        
        # Remove layers and raw from AnnData object (avoid conflicts with h5seurat)
        if self.adata.layers["raw"] is not None:
            del self.adata.layers["raw"]
        self.adata.raw = None
        
        # Strip mapped h5 of extension - keep only name for internal seurat h5 conversions
        self.mapped_h5_name = os.path.splitext(self.mapped_h5)[0]
        
    def _output_temp_h5ad(self):
        # Check if temp exists, if not, make dir
        if not os.path.exists("tmp"):
            os.makedirs("tmp")

        # Output temporary file with data 
        self.filename = ''.join(str(uuid.uuid4()).split("-"))
        self.file = "{filename}.h5ad".format(filename = self.filename)
        self.adata.write_h5ad(os.path.join("tmp", self.file), compression = "gzip")
    
    def _seurat_refmap(self):
        # Call subprocess and call R script
        refmap_script = \
            "Rscript src/R/seurat_reference_map.R {ref_h5} tmp/{tempfile} {tempname} {out_name} --verbose".format(
                ref_h5 = self.reference_h5,
                tempfile = self.file,
                tempname = self.filename,
                out_name = self.mapped_h5_name
            )
            
        self.sp_refmap = subprocess.run(refmap_script, shell = True, text = True, capture_output = True)
        if self.sp_refmap.returncode != 0:
            raise Exception(
                "Subprocess call returned nonzero exit code - call: {call} \n Output: {output}".format(
                    call = self.sp_refmap.stderr,
                    output = self.sp_refmap.stdout
                )
            )
        
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
            
    def refmap(self):
        # Perform workflow and return reference mapped anndata object
        self._load()
        self._format()
        self._output_temp_h5ad()
        self._seurat_refmap()
        self._clean_files()