# Iniquitate Pipeline

## Downsampling-based perturbation experiments for single-cell RNA sequencing data

## Please note that this repository is still under construction
***

### Reproducing the paper analysis 

1. Clone the GitHub repository:

```
git clone https://github.com/hsmaan/Iniquitate.git
```

2. Download the resources utilized in the study, extract and move to Iniquitate directory:

```
wget [ ]
tar -xzvf resources.tar.gz 
mv resources Iniquitate
```

3. Run the different configurations utilized in the study through the Snakemake pipeline:

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


4. Analyze the results using the R scripts:

    - Begin by installing the analysis conda environment:
    ```
    cd Iniquitate/workflows/envs
    mamba env create -f analysis.yaml
    ```

    - Run the Rscripts based on their order, through RStudio or the commandline 
    ```
    conda activate iniq_analysis
    Rscript 05_Iniq_Control_Fig_2_Analysis_Plots.R
    Rscript 06_Iniq_Control_Fig_2_Analysis_Stat_Tests.R
    ...
    ```

**It is not possible** to re-run all of the perturbation experiments and downstream analyses in a reasonable amount of time without high-performance computing (HPC). It is highly recommended that the workflow is parallelized over HPC.

***
### Custom configuration setup 

**Under construction**