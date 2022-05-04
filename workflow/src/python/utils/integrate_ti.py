from functools import reduce
import gc
import random 

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ann
import scvi 
import bbknn 
import torch

# Undoing scvi's random seed setting
random.seed(None)
np.random.seed(None)
torch.manual_seed(random.randint(1, 10000000000000000000))

from utils.seurat_integrate import SeuratIntegrate
from utils.liger_integrate import LigerIntegrate

class IntegrationPAGA:
    """Class for integrating scRNA-seq data and subsequently performing PAGA trajectory inference."""
    
    def __init__(self, adata, root_celltype, gpu = True):
        """
        Args:
            adata (AnnData): AnnData object to be utilized in integration methods.
                Assumes that the counts being input are unnormalized (raw counts),
                and that raw counts are stored in "counts" layer, and batch covariate
                is available.
            root_celltype (str): Celltype to be used as the root celltype for trajectory
                inference via diffusion pseudotime.
            gpu (bool): Whether or not to use GPU for scVI.
        """
        self.adata = adata
        self.root_celltype = root_celltype
        # Check anndata object 
        if not isinstance(adata, ann.AnnData):
            raise Exception("Please input an AnnData object.")
        # Check if gpu is available
        if gpu is True:
            if torch.cuda.is_available():
                self.gpu = True
            else:
                raise Exception("GPU not available. Please set gpu = False.")
        else:
            self.gpu = False
            
    def unintegrated(self, n_neighbors = 15, n_pcs = 20, num_hvgs = 2500):
        print("Performing unintegrated.." + "\n")
        aunint = self.adata.copy()
        sc.pp.normalize_total(
            aunint,
            target_sum = 1e4
        )
        sc.pp.log1p(aunint)
        sc.pp.highly_variable_genes(
            aunint,
            n_top_genes = num_hvgs,
            flavor = "seurat"
        )
        sc.pp.pca(aunint, svd_solver="arpack")
        sc.pp.neighbors(
            aunint,
            n_neighbors = n_neighbors,
            n_pcs = n_pcs,
        )
        sc.tl.leiden(aunint)
        sc.tl.umap(aunint)
        sc.tl.paga(
            aunint, 
            groups = "celltype"
        )
        aunint.uns["iroot"] = np.flatnonzero(
            aunint.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(aunint)
        print("Done!" + "\n")
        return aunint

    def scvi_integrate(self, n_neighbors = 15, n_pcs = 20):
        print("Performing scVI integration.." + "\n")
        ascvi = self.adata.copy()
        scvi.data.setup_anndata(ascvi, batch_key = "batch")
        vae = scvi.model.SCVI(ascvi)
        vae.train(use_gpu = self.gpu)
        ascvi.obsm["X_scVI"] = vae.get_latent_representation()
        sc.pp.neighbors(
            ascvi,
            n_neighbors = n_neighbors,
            n_pcs = n_pcs,
            use_rep = "X_scVI"
        )
        sc.tl.leiden(ascvi)
        sc.tl.umap(ascvi)
        sc.tl.paga(
            ascvi, 
            groups = "celltype"
        )
        ascvi.uns["iroot"] = np.flatnonzero(
            ascvi.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(ascvi)
        print("Done!" + "\n")
        return ascvi
    
    def harmony_integrate(self, n_neighbors = 15, n_pcs = 20, num_hvgs = 2500):
        print("Performing Harmony integration.." + "\n")
        aharmony = self.adata.copy()
        sc.pp.normalize_total(
            aharmony,
            target_sum = 1e4
        )
        sc.pp.log1p(aharmony)
        sc.pp.highly_variable_genes(
            aharmony,
            n_top_genes = num_hvgs,
            flavor = "seurat"
        )
        sc.pp.pca(aharmony, svd_solver="arpack")
        sc.external.pp.harmony_integrate(
            aharmony,
            key = "batch",
            random_state = None
        )
        sc.pp.neighbors(
            aharmony,
            n_neighbors = n_neighbors,
            n_pcs = n_pcs,
            use_rep = "X_pca_harmony"
        )
        sc.tl.leiden(aharmony)
        sc.tl.umap(aharmony)
        sc.tl.paga(
            aharmony, 
            groups = "celltype"
        )
        aharmony.uns["iroot"] = np.flatnonzero(
            aharmony.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(aharmony)
        print("Done!" + "\n")
        return aharmony
    
    def bbknn_integrate(self, n_pcs = 20, num_hvgs = 2500, metric = "euclidean"):
        print("Performing BBKNN integration.." + "\n")
        abbknn = self.adata.copy()
        sc.pp.normalize_total(
            abbknn,
            target_sum = 1e4
        )
        sc.pp.log1p(abbknn)
        sc.pp.highly_variable_genes(
            abbknn,
            n_top_genes = num_hvgs,
            flavor = "seurat"
        )
        sc.pp.pca(abbknn, svd_solver = "arpack")
        if metric == "euclidean":
            bbknn.bbknn(
                abbknn,
                approx = False,
                metric = "euclidean",
                batch_key = "batch",
                n_pcs = n_pcs,
                pynndescent_random_state = None
            )
        elif metric == "angular":
            bbknn.bbknn(
                abbknn,
                approx = True,
                metric = "angular",
                batch_key = "batch",
                n_pcs = n_pcs,
                pynndescent_random_state = None
            )
        else:
            raise Exception(
                "Please enter either 'euclidean' or 'angular' for 'metric'"
            )
        sc.tl.leiden(abbknn)
        sc.tl.umap(abbknn)
        sc.tl.paga(
            abbknn, 
            groups = "celltype"
        )
        abbknn.uns["iroot"] = np.flatnonzero(
            abbknn.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(abbknn)
        print("Done!" + "\n")
        return abbknn
    
    def scanorama_integrate(self, n_neighbors = 15, n_pcs = 20, num_hvgs = 2500):
        print("Performing Scanorama integration.." + "\n")
        ascanorama = self.adata.copy()
        sc.pp.normalize_total(
            ascanorama,
            target_sum = 1e4
        )
        sc.pp.log1p(ascanorama)
        sc.pp.highly_variable_genes(
            ascanorama,
            n_top_genes = num_hvgs,
            flavor = "seurat"
        )
        sc.pp.pca(ascanorama, svd_solver="arpack")
        sc.external.pp.scanorama_integrate(
            ascanorama,
            key = "batch"
        )
        sc.pp.neighbors(
            ascanorama,
            n_neighbors = n_neighbors,
            n_pcs = n_pcs,
            use_rep = "X_scanorama"
        )
        sc.tl.leiden(ascanorama)
        sc.tl.umap(ascanorama)
        sc.tl.paga(
            ascanorama, 
            groups = "celltype"
        )
        ascanorama.uns["iroot"] = np.flatnonzero(
            ascanorama.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(ascanorama)
        print("Done!" + "\n")
        return ascanorama
    
    def seurat_integrate(self,int_type = "CCA", n_neighbors = 15, n_pcs = 20):
        print("Performing Seurat integration.." + "\n")
        aseurat = self.adata.copy()
        sc.pp.normalize_total(
            aseurat,
            target_sum = 1e4
        )
        sc.pp.log1p(aseurat)
        seurat_integrate = SeuratIntegrate(
            adata = aseurat,
            int_type = int_type
        )
        aseurat_int = seurat_integrate.integrate() # Create seurat integrated anndata object
        sc.pp.pca(aseurat_int, svd_solver = "arpack")
        sc.pp.neighbors(
            aseurat_int,
            n_neighbors = n_neighbors,
            n_pcs = n_pcs
        )
        sc.tl.leiden(aseurat_int)
        sc.tl.umap(aseurat_int)
        sc.tl.paga(
            aseurat_int, 
            groups = "celltype"
        )
        aseurat_int.uns["iroot"] = np.flatnonzero(
            aseurat_int.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(aseurat_int)
        # Append seurat integrated data to original adata object
        aseurat.obs["leiden"] = aseurat_int.obs["leiden"]
        aseurat.obsm["X_pca"] = aseurat_int.obsm["X_pca"]
        aseurat.obsm["X_umap"] = aseurat_int.obsm["X_umap"]
        aseurat.obsm["seurat_hvg"] = aseurat_int.X
        aseurat.obs["dpt_pseudotime"] = aseurat_int.obs["dpt_pseudotime"]
        aseurat.uns["iroot"] = aseurat_int.uns["iroot"]
        aseurat.uns["paga"] = aseurat_int.uns["paga"]
        aseurat.obsp["distances"] = aseurat_int.obsp["distances"]
        aseurat.obsp["connectivities"] = aseurat_int.obsp["connectivities"]
        print("Done!" + "\n")
        return aseurat
        
    def liger_integrate(self, n_neighbors = 15, n_pcs = 20):
        print("Performing LIGER integration.." + "\n")
        aliger = self.adata.copy()
        # Don't normalize for LIGER (R script normalizes)
        # sc.pp.normalize_total(
        #     aliger,
        #     target_sum = 1e4
        # )
        # Don't log-transform for LIGER
        # sc.pp.log1p(aliger)   
        liger_integrate = LigerIntegrate(
            adata = aliger,
        )
        aliger_int = liger_integrate.integrate() # Create liger integrated anndata object
        sc.pp.neighbors(
            aliger_int,
            n_neighbors = n_neighbors,
            n_pcs = n_pcs,
            use_rep = "X_liger"
        )
        sc.tl.leiden(aliger_int)
        sc.tl.umap(aliger_int)
        sc.tl.paga(
            aliger_int, 
            groups = "celltype"
        )
        aliger_int.uns["iroot"] = np.flatnonzero(
            aliger_int.obs["celltype"] == self.root_celltype
        )[0]
        sc.tl.dpt(aliger_int)
        # Append Liger integrated data to original adata object
        aliger.obs["leiden"] = aliger_int.obs["leiden"]
        aliger.obsm["X_umap"] = aliger_int.obsm["X_umap"]
        aliger.obs["dpt_pseudotime"] = aliger_int.obs["dpt_pseudotime"]
        aliger.uns["iroot"] = aliger_int.uns["iroot"]
        aliger.uns["paga"] = aliger_int.uns["paga"]
        aliger.obsp["distances"] = aliger_int.obsp["distances"]
        aliger.obsp["connectivities"] = aliger_int.obsp["connectivities"]
        print("Done!" + "\n")
        return aliger