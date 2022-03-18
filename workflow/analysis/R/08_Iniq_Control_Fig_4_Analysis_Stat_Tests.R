library(data.table)
library(tidyverse)
library(reshape2)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Change to results dir for control data 
setwd("../../../results/control/")

##### Statistical analysis of PBMC 2 batch balanced data - baseline #####

# Load in and concatenate imbalance summary files 
setwd("imbalance_summaries/")
imba_files <- list.files()
imba_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_files,
  value = TRUE
)
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_files,
  value = TRUE
)
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)
gc()

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
clus_concord_files <- grep(
  "pbmc_2_batch_base_balanced",
  clus_concord_files,
  value = TRUE
)
clus_concord_loaded <- lapply(clus_concord_files, fread)
clus_concord_concat <- Reduce(rbind, clus_concord_loaded)
gc()

# Load in and concatenate dge concordance summaries
setwd("../dge_concord_stats/")
dge_files <- list.files()
dge_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_files,
  value = TRUE
)
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)
gc()

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("../dge_ranking_stats_marker_sub")
dge_rank_files <- list.files()
dge_rank_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_rank_files,
  value = TRUE
)
dge_rank_loaded <- lapply(dge_rank_files, fread)
dge_rank_concat <- Reduce(rbind, dge_rank_loaded)
gc()

# Load in markers and corresponding celltypes 
setwd("../marker_results/")
base_marker_genes <- fread(
  "pbmc_2_batch_base_balanced_preintegration_marker_selection.tsv"
) 

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_files <- grep(
  "pbmc_2_batch_base_balanced",
  knn_files,
  value = TRUE
)
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
gc()

# Change to top level dir 
setwd("../../..")

### Statistical tests I - cluster number concordance before and after
### downsampling 

# Merge together imbalance and cluster number results
imba_clus_merged <- merge(
  clus_concat,
  imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
)
imba_clus_merged <- distinct(imba_clus_merged)

# Indicate which samples are controls and which are real runs
imba_clus_merged$type <- ifelse(
  imba_clus_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_clus_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format column names for an ANOVA test 
colnames(imba_clus_merged) <- plyr::mapvalues(
  colnames(imba_clus_merged),
  from = c(
    "Cluster number",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "cluster_number",
    "method",
    "downsampled_celltypes"
  )
)

# Perform an ANOVA test for concordance of number of clusters vs type, after
# factoring in method and celltype downsampled - format and save the result
clus_model_fit <- lm(
  as.formula(
    paste0(
      "cluster_number", 
      "~",
      "method+",
      "downsampled_celltypes+",
      "type"
    )
  ),
  data = imba_clus_merged
)
clus_anova_result <- anova(clus_model_fit, test = "F")
clus_anova_result_dt <- as.data.table(clus_anova_result, keep.rownames = TRUE)
colnames(clus_anova_result_dt)[1] <- "Covariate"
clus_anova_result_dt$dataset_name <- "PBMC 2 batch base balanced"
clus_anova_result_dt$metric <- "Leiden cluster number"
clus_anova_result_dt$last_covariate <- "type"
fwrite(
  clus_anova_result_dt,
  paste0(
    "outs/control/results/",
    "08_pbmc_base_clus_num_leiden_aov_results_ctrl_method_celltype_ds.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

### Statistical tests II - analysis of marker gene DGE rank change 

# Merge together imbalanced and dge data 
imba_dge_merged <- merge(
  imba_concat,
  dge_rank_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
) 
imba_dge_merged <- distinct(imba_dge_merged)

# Indicate which samples are controls and which are real runs
imba_dge_merged$type <- ifelse(
  imba_dge_merged$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_dge_merged$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
)

# Format column names for an ANOVA test 
colnames(imba_dge_merged) <- plyr::mapvalues(
  colnames(imba_dge_merged),
  from = c(
    "Gene",
    "Max rank",
    "Method",
    "Downsampled celltypes"
  ),
  to = c(
    "gene",
    "max_rank",
    "method",
    "downsampled_celltypes"
  )
)

## Perform an overall ANOVA test for concordance of max rank of marker gene vs 
## type, after factoring in marker gnee, method and celltype downsampled - format 
## and save the result
dge_model_fit <- lm(
  as.formula(
    paste0(
      "max_rank", 
      "~",
      "gene+",
      "method+",
      "downsampled_celltypes+",
      "type"
    )
  ),
  data = imba_dge_merged
)
dge_anova_result <- anova(dge_model_fit, test = "F")
dge_anova_result_dt <- as.data.table(dge_anova_result, keep.rownames = TRUE)
colnames(dge_anova_result_dt)[1] <- "Covariate"
dge_anova_result_dt$dataset_name <- "PBMC 2 batch base balanced"
dge_anova_result_dt$metric <- "Max DGE rank marker gene"
dge_anova_result_dt$last_covariate <- "type"
fwrite(
  dge_anova_result_dt,
  paste0(
    "outs/control/results/",
    "08_pbmc_base_dge_rank_aov_results_ctrl_gene_method_celltype_ds.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## Perform this ANOVA test independantly for each marker gene, now, so removing
## the conditioning on the marker gene, but keeping the conditioning on the 
## method and downsampled celltypes 
marker_gene_list <- as.list(unique(imba_dge_merged$gene))

# Create ANOVA function to perform anova for each marker gene
anova_marker_gene <- function(
  dataset, 
  dataset_name,
  marker_gene
){
  # Subset dataset for the given marker
  dataset_marker_sub <- dataset[
    which(dataset$gene %in% marker_gene)
  ]
  
  # Fit ANOVA model
  model_fit <- lm(
    as.formula(
      paste0(
        "max_rank", 
        "~",
        "method+",
        "downsampled_celltypes+",
        "type"
      )
    ),
    data = dataset_marker_sub
  )
  anova_result <- anova(model_fit, test = "F")
  
  # Format results and return
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$dataset_name <- dataset_name
  anova_result_dt$marker_gene <- marker_gene
  anova_result_dt$metric <- "Max DGE rank"
  anova_result_dt$last_covariate <- "type"
  return(anova_result_dt)
}

# Iterate over marker genes, get results and concatenate and save
marker_aov_results <- lapply(marker_gene_list, function(x) {
  aov_summary_df <- anova_marker_gene(
    dataset = imba_dge_merged,
    dataset_name = "PBMC 2 batch base balanced",
    marker_gene = x
  )
})

# Concatenate results and save 
marker_aov_results_concat <- Reduce(rbind, marker_aov_results)
fwrite(
  marker_aov_results_concat,
  paste0(
    "outs/control/results/",
    "08_pbmc_base_dge_specific_rank_aov_results_ctrl_method_celltype_ds.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

