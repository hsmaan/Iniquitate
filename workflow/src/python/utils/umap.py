import pandas as pd
import numpy as np
import seaborn as sns
import colorcet as cc
from natsort import natsorted


class Umap:
    """Class for plotting results of integration experiments"""

    def __init__(self, coords, clustering, subset_name = None):
        """
        Args:
            coords (dictionary): coordinates of umap in numpy format where
                keys correspond to following integration methods -
                'bbknn', 'harmony', 'scanorama', 'seurat', and 'scvi'.
            clustering (dictionary): leiden or celltype clustering in numpy
                format of integrated where keys correspond to following
                integration methods - 'bbknn', 'harmony', 'scanorama',
                'seurat', and 'scvi'.
            subset_name (string): name of subset being utilized for clustering
                comparisons (e.g. batch, celltype).
        """
        self.clustering_harmony = clustering.get("harmony")
        self.clustering_scvi = clustering.get("scvi")
        self.clustering_bbknn = clustering.get("bbknn")
        self.clustering_scanorama = clustering.get("scanorama")
        self.clustering_seurat = clustering.get("seurat")
    
        self.umap_harmony = coords.get("harmony")
        self.umap_scvi = coords.get("scvi")
        self.umap_bbknn = coords.get("bbknn")
        self.umap_scanorama = coords.get("scanorama")
        self.umap_seurat = coords.get("seurat")
        
        if subset_name is not None:
            self.subset_name = subset_name
        else:
            self.subset_name = "Subset"
            
        sns.set_style("ticks")

    def df_get(self, subset, clustering, coords, category = None):
        df = pd.DataFrame({
            "Subset" : np.repeat(subset, len(clustering)),
            "UMAP 1" : coords[:, 0],
            "UMAP 2" : coords[:, 1]
        })
        df["Clustering"] = pd.Categorical(
            clustering, categories=category, ordered=True
        )
        return df

    def umap_df(self):
        subset_list = [
            "harmony",
            "scvi",
            "bbknn",
            "scanorama",
            "seurat"
        ]
        clustering_list = [
            self.clustering_harmony,
            self.clustering_scvi,
            self.clustering_bbknn,
            self.clustering_scanorama,
            self.clustering_seurat
        ]
        clustering_unique = natsorted(np.unique(np.concatenate(clustering_list)))
        coords_list = [
            self.umap_harmony,
            self.umap_scvi,
            self.umap_bbknn,
            self.umap_scanorama,
            self.umap_seurat
        ]
        umap_dfs = [
            self.df_get(i, j, k, category = clustering_unique) for i, j, k in zip(
                subset_list,
                clustering_list,
                coords_list
            )
        ]
        self.umap_concat = pd.concat(umap_dfs)

    def umap_plot(self, show_plot = False):
        self.umap_df()
        palette = cc.glasbey_bw[0:len(np.unique(self.umap_concat["Clustering"]))]
        self.umap_plt = sns.FacetGrid(
            self.umap_concat,
            col = "Subset",
            col_wrap = 3,
            hue = "Clustering",
            palette = palette
        )
        self.umap_plt.map(
            sns.scatterplot,
            "UMAP 1",
            "UMAP 2",
            s = 5,
            alpha = 0.5
        )
        self.umap_plt.add_legend(markerscale = 3, title = self.subset_name)
        if show_plot is True:
            return self.umap_plt                

    def save_umap(self, save_dir, dpi = 300):
        try:
            self.umap_plt.savefig(
                save_dir,
                dpi = dpi
            )
        except:
            self.umap_plot()
            self.umap_plt.savefig(
                save_dir,
                dpi = dpi
            )