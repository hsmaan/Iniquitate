# Iniquitate Pipeline <!-- omit in toc -->

This repository corresponds to the analysis and pipeline outlined in [Characterizing the impacts of dataset imbalance on single-cell data integration](https://www.nature.com/articles/s41587-023-02097-9). If you use the integration guidelines or any aspect of this pipeline, please consider [citing our work](#citation-information).

### Downsampling-based perturbation experiments for single-cell RNA sequencing integration <!-- omit in toc -->

***

## Table of contents  <!-- omit in toc -->
- [Using the imbalanced integration guidelines](#using-the-imbalanced-integration-guidelines)
- [Reproducing the paper analysis](#reproducing-the-paper-analysis)
- [Custom data perturbation configuration setup](#custom-data-perturbation-configuration-setup)
- [Citation information](#citation-information)


### Using the imbalanced integration guidelines

A separate README for the imbalanced integration guidelines, with full environment installation instructions are in the `docs` folder. 

### Reproducing the paper analysis 

Please note that `mamba` and `snakemake` are required to run the pipeline through `conda`. After installing `conda` (https://conda.io/projects/conda/en/latest/user-guide/install/index.html), please add `mamba` (https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) to the base environment, as well as `snakemake` (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) in base or a new environment:

```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

The environments necessary to run the pipeline will be automatically installed by snakemake. The only exception is the `analysis` pipeline, which is outlined in step 5 below.

The following steps are necessary to reproduce the paper analysis:


1. Clone the GitHub repository:

```
git clone https://github.com/hsmaan/Iniquitate.git
```

2. Download the resources utilized in the study, extract and move to Iniquitate directory. The data can be downloaded from Figshare or Google Drive:

    Figshare:
    https://doi.org/10.6084/m9.figshare.24625302.v1


    Google Drive Link:
    - Download the data from https://drive.google.com/file/d/1gWsYEI_u0Bn-7liar1XmvcrFqdt3IHjV/view?usp=sharing
    - Alternatively, you can use gdown (https://github.com/wkentaro/gdown) if a command-line download is needed/desired

    After downloading:
    ```
    tar -xzvf resources.tar.gz
    # If the folder name is not named `resources`, change the name via mv [file_folder] resources
    mv resources Iniquitate
    ```

4. Run the different configurations utilized in the study through the Snakemake pipeline:

    - Change the configuration option at the top of `workflow/Snakefile`. The following configs were utilized for different analyses in the study:

        - config_control 
        - config_lowcap_modified
        - config_pdac_comp

    - Run the Snakemake pipeline specific to the selected config:

    ```
    snakemake --unlock 
    snakemake -j 1000 \
        --use-conda \
        --cluster-config cluster.json \
        --cluster "sbatch \
            --mem={cluster.mem} \
            --gres=gpu:{cluster.gpu} \
            -t {cluster.time} \
            -p {cluster.partition} \
            -c {threads}" \
        --restart-times 0 \
        --latency-wait 300 \
        --keep-going \
        --rerun-incomplete 
    ```

    Note that the above Snakemake run utilizes a `workflow/cluster.json` configuration file and HPC parallelization of the various steps in the pipeline. Users will need to create a `cluster.json` file specific to their HPC setup that has resources for all of the rules in `workflow/Snakefile`. Alternatively, users can also choose to employ Snakemake profiles. Details can be found here: https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html.

    Further, currently all of the temporary integration files will by default be stored in `/tmp` based on the `TMPDIR` variable in `workflow/envs/integrate.yaml`. These files can be quite large and pile up during runtime, even though they are being deleted at the end of each run, and we therefore recommend that users set this directory to one that contains enough space - at least 100 GB. 

5. Analyze the results using the R and python scripts/notebooks:

    - Begin by installing the analysis conda environment:
    ```
    cd Iniquitate/workflows/envs
    mamba env create -f analysis.yaml
    ```

    - First run the python notebook to completion for metric analysis:
    ```
    conda activate iniq_analysis
    jupyter-notebook 01_Fig_7_Imbal_Metric_Analysis.ipynb
    ```
 
    - Run the Rscripts based on their order, through RStudio or the commandline 
    ```
    conda activate iniq_analysis
    Rscript 05_Iniq_Control_Fig_2_Analysis_Plots.R
    Rscript 06_Iniq_Control_Fig_2_Analysis_Stat_Tests.R
    ...
    ```

**It is not possible** to re-run all of the perturbation experiments and downstream analyses in a reasonable amount of time without high-performance computing (HPC). It is highly recommended that the workflow is parallelized over HPC.

It is also recommended to **Run the R and python analysis notebooks** in an HPC environment as well, because some of the steps are memory-intensive. Particularly, we don't recommend running Rscripts 08 or 09 without HPC, as they are time-intensive sampling experiments. 

***
### Custom data perturbation configuration setup 

The same software requirements for the paper analysis apply to custom data perturbation experiments (conda, mamba, snakemake). Please see the first part of **Reproducing the paper analysis** for details on how to install these requirements.

The following steps are necessary to use a custom dataset:

1) Process the batches/samples as necessary and convert to `h5ad` format. Create a folder in `resources/h5ad_files/int_datasets/` (e.g. `resources/h5ad_files/int_datasets/custom_data`) and move batches/samples to this location.

2) Modify the custom dataset configuration file (`workflow/configs/config_custom.json`). This .json file has the following parameters:

    - `config_name` - this can be left as custom, or changed to a different name, but this same name must be used when modifying the `Snakefile`
    - `int_datasets` - this is a nested dictionary of the datasets to be used in the downsampling experiments. In this case, it's best to use the name of the dataset in `resources/h5ad_files/int_datasets` as the top level name
        - `data_folder` - this value should be the same as the folder containing the batches/samples in `resources/h5ad_files/int_datasets/`
        - `ds_celltypes` - this value indicates how many cell-types to downsample and/or ablate in each run
        - `ds_proportions` - this value indicates what proportion of cell-types remain in the downsampled batch after the pertubation. In this case, 0 would indicate ablation, and 0.1 would indicate the same level of downsampling used for the main experiments in the study
        - `num_batches` - the number of batches to downsample in each perturbation. The 0 option is included here so that control experiments are possible (no perturbation to any batches).
        - `repetitions` - how many experiments based on the given grid of `ds_celltypes`, `ds_proportions`, and `num_batches` to perform. 200 is a good starting point. If the space of possible cell-types is very large (n celltypes > 20), then it may be useful to increase this value to ensure each cell-type is downsampled/ablated enough times in the total number of runs.

    - `int_ti_datasets` - if any datasets have an underlying trajectory, and PAGA-based integration is needed to be done, then they should be added here with the same options indicated in `int_datasets`
    - `query_to_reference` - A "Yes" or "No" option indicating whether or not to perform query-to-reference experiments. Currently this functionality is not available, but custom query-to-reference setups will be available soon. This should be left as "No".
    - `celltype_list` - If the user has a list of specific cell-types to downsample (and not others), they can be included her as a json list of strings based on their names. We don't recommend specifying certain cell-types, as a-priori knowledge of the effects of downsampling/perturbation may not be accurate. 

3) Modify the `Snakefile` in `workflow\Snakefile` at line 3, in reference to the name of the specific config being used. In the example given, the `configfile` line would be changed to:

    - `configfile: "configs/config_custom.json"`

4) Run snakemake:

    ```
    snakemake --unlock 
    snakemake -j 1000 \
        --use-conda \
        --cluster-config cluster.json \
        --cluster "sbatch \
            --mem={cluster.mem} \
            --gres=gpu:{cluster.gpu} \
            -t {cluster.time} \
            -p {cluster.partition} \
            -c {threads}" \
        --restart-times 0 \
        --latency-wait 300 \
        --keep-going \
        --rerun-incomplete 
    ```

    Note that the above Snakemake run utilizes a `workflow/cluster.json` configuration file and HPC parallelization of the various steps in the pipeline. Users will need to create a `cluster.json` file specific to their HPC setup that has resources for all of the rules in `workflow/Snakefile`.
    
5) Analyze the result files using the R scripts - all of the R scripts in the `workflow/analysis/R` folder can be used to analyze the results of the perturbation experiments. **Currently, an exhaustive list of custom analysis scripts is not available, but the existing scripts can be modified to suit the needs of the user**. We provide one custom script for analysis of KNN classification accuracy in the `workflow/analysis/R/knn_example.R` file. Please note that this file still needs to be modified in the appropriate input locations, which are indicated in the comments of the file.

### Citation information

Maan, H. et al. (2024) ‘Characterizing the impacts of dataset imbalance on single-cell data integration’, Nature biotechnology. Available at: https://doi.org/10.1038/s41587-023-02097-9.
