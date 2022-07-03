library(data.table)
library(tidyverse)
library(reshape2)
library(EMT)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Set seed for any sampling done 
set.seed(42)

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
clus_concat <- clus_concat[
  clus_concat$Method != "liger"
]
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
clus_concord_concat <- clus_concord_concat[
  clus_concord_concat$`Method 1` != "liger"
]
clus_concord_concat <- clus_concord_concat[
  clus_concord_concat$`Method 2` != "liger"
]
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
dge_concat <- dge_concat[dge_concat$`Method 1` != "liger"]
dge_concat <- dge_concat[dge_concat$`Method 2` != "liger"]
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
dge_rank_concat <- dge_rank_concat[
  dge_rank_concat$Method != "liger"
]
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
knn_concat <- knn_concat[knn_concat$Method != "liger"]
gc()

# Change to top level dir 
setwd("../../..")

# Format marker genes to indicate which celltypes they belong to (due to 
# duplicates)
base_marker_gene_dup_added <- base_marker_genes %>%
  group_by(`Top 10 marker genes (union across batches)`) %>%
  summarize(Celltype = paste0(unique(Celltype), collapse = ', ')) %>%
  as.data.frame
colnames(base_marker_gene_dup_added) <- c(
  "marker_gene", "celltype"
)

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

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_clus_unique_celltypes <- unique(imba_clus_merged$`Downsampled celltypes`)
imba_clus_unique_celltypes <- imba_clus_unique_celltypes[
  imba_clus_unique_celltypes %ni% "None"
]
imba_clus_none_celltype_len <- length(
  which(
    imba_clus_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_clus_none_celltype_draw <- sample(
  imba_clus_unique_celltypes,
  imba_clus_none_celltype_len,
  replace = TRUE
)
imba_clus_merged$`Downsampled celltypes`[
  which(
    imba_clus_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_clus_none_celltype_draw

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
    "09_pbmc_base_clus_num_leiden_aov_results_ctrl_method_celltype_ds_no_liger.tsv"
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

# Fill in 'None' value for downsampled celltypes with a random 
# draw from the given celltypes - to block correlation effect of 
# downsampled celltype and `type` covariate in later use 
imba_dge_unique_celltypes <- unique(imba_dge_merged$`Downsampled celltypes`)
imba_dge_unique_celltypes <- imba_dge_unique_celltypes[
  imba_dge_unique_celltypes %ni% "None"
]
imba_dge_none_celltype_len <- length(
  which(
    imba_dge_merged$`Downsampled celltypes` %in% "None"
  )
)
imba_dge_none_celltype_draw <- sample(
  imba_dge_unique_celltypes,
  imba_dge_none_celltype_len,
  replace = TRUE
)
imba_dge_merged$`Downsampled celltypes`[
  which(
    imba_dge_merged$`Downsampled celltypes` %in% "None"
  )
] <- imba_dge_none_celltype_draw

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
    "09_pbmc_base_dge_rank_aov_results_ctrl_gene_method_celltype_ds_no_liger.tsv"
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

# Concatenate results, correct for FDR, and save 
marker_aov_results_concat <- Reduce(rbind, marker_aov_results)
marker_aov_results_concat$`FDR_q` <- p.adjust(
  marker_aov_results_concat$`Pr(>F)`
)
fwrite(
  marker_aov_results_concat,
  paste0(
    "outs/control/results/",
    "09_pbmc_base_dge_specific_rank_aov_results_ctrl_method_celltype_ds_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## Get summary of linear model fit for each marker gene, after accounting for 
## method - determine which celltypes downsampled are most significantly 
## correlated for a given marker gene in this setup

# Create function that returns summary of a linear model for a 
# marker gene 
lm_marker_gene <- function(
  dataset, 
  dataset_name,
  marker_gene
){
  # Subset dataset for the given marker
  dataset_marker_sub <- dataset[
    which(dataset$gene %in% marker_gene)
  ]
  
  # Fit linear regression model
  model_fit <- lm(
    as.formula(
      paste0(
        "max_rank", 
        "~",
        "method+",
        "downsampled_celltypes"
      )
    ),
    data = dataset_marker_sub
  )
  
  # Format summary results - append relevant information
  model_fit_summary <- summary(model_fit)$coefficients
  model_fit_summary_dt <- as.data.table(
    as.data.frame(
      model_fit_summary
    ),
    keep.rownames = TRUE
  )
  colnames(model_fit_summary_dt)[1] <- "coefficient"
  model_fit_summary_dt$marker_gene <- marker_gene
  model_fit_summary_dt$metric <- "Max DGE rank"
  model_fit_summary_dt$dataset_name <- "PBMC 2 batch base balanced"
  model_fit_summary_dt_merged <- merge(
    model_fit_summary_dt,
    base_marker_gene_dup_added,
    by = "marker_gene"
  )
  
  # Indicate if the top `celltype downsampled` coefficient corresponds to the 
  # celltype the marker gene is indicative of 
  model_fit_summary_dt_merged_coeff <- model_fit_summary_dt_merged[
    grep("downsampled_celltypes", model_fit_summary_dt_merged$coefficient)
  ]
  top_downsampled_celltype_idx <- which.max(
    abs(
      model_fit_summary_dt_merged_coeff$Estimate
    )
  )
  top_downsampled_celltype_uf <- model_fit_summary_dt_merged_coeff$coefficient[
    top_downsampled_celltype_idx
  ]
  top_downsampled_celltype_formatted <- str_split_fixed(
    top_downsampled_celltype_uf,
    "downsampled_celltypes",
    2
  )[,2]
  model_fit_summary_dt_merged$top_ds_celltype <- top_downsampled_celltype_formatted
  
  marker_associated_celltypes <- base_marker_gene_dup_added[
    base_marker_gene_dup_added$marker_gene == marker_gene, 
  ]$celltype
  model_fit_summary_dt_merged$marker_assoc_celltypes <- marker_associated_celltypes
  model_fit_summary_dt_merged$same_ds_lm_celltype <- ifelse(
    grepl(
      top_downsampled_celltype_formatted,
      marker_associated_celltypes,
      ignore.case = TRUE
    ),
    "Yes",
    "No"
  )
  
  # Return summary datatable
  return(model_fit_summary_dt_merged)
}

# Iterate over marker genes and get results 
marker_lm_results <- lapply(marker_gene_list, function(x) {
  aov_summary_df <- lm_marker_gene(
    dataset = imba_dge_merged,
    dataset_name = "PBMC 2 batch base balanced",
    marker_gene = x
  )
})

# Concatenate and FDR correct the results 
marker_lm_summaries_concat <- Reduce(rbind, marker_lm_results)
marker_lm_summaries_concat$q_val <- p.adjust(
  marker_lm_summaries_concat$`Pr(>|t|)`
)

# Save results for coefficients and top associated celltypes 
fwrite(
  marker_lm_summaries_concat,
  paste0(
    "outs/control/results/",
    "09_pbmc_base_dge_rank_coeffs_celltype_ds_method_importance_ranks_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

## Get enrichment of top coefficient for celltype downsampled being the 
## celltype associated with the marker gene being perturbed - perform
## multinomial test for each downsampled celltype against all others

# Create function that takes a subset of marker associated celltypes 
# and returns a multinomial sample test for top associated downsampled
# celltype coefficient value
marker_celltype_multinomial <- function(
  celltype_name,
  dataset
) {
  # Get list of all celltypes (no-mixed) and their length
  celltypes_all <- unique(base_marker_genes$Celltype)
  celltypes_all_len <- length(celltypes_all)
  
  # Subset data for the given celltype 
  data_sub <- dataset[dataset$marker_assoc_celltypes %in% celltype_name,]
  
  # Get counts vector for all the different celltypes 
  celltype_counts_vector <- table(data_sub$top_ds_celltype)
  
  # Add any missing celltypes to vector
  missing_celltypes <- celltypes_all[
    celltypes_all %ni% names(celltype_counts_vector)
  ]
  for (celltype in missing_celltypes) {
    celltype_counts_vector[celltype] <- 0
  }
  
  # Perform multinomial test for observed celltypes in top downsampled
  # results
  prob <- rep(1/celltypes_all_len, celltypes_all_len)
  observed <- as.vector(celltype_counts_vector)
  p_val <- multinomial.test(
    observed, 
    prob, 
    ntrial = 1e6,
    MonteCarlo = TRUE
  )$p.value
  print(p_val)
  
  # Create results dataframe, sort by celltype, and return
  res_df <- as.data.frame(as.array(celltype_counts_vector))
  colnames(res_df) <- c(
    "Top_celltype_downsampled_coefficient",
    "Observed_value_for_top_celltype_downsampled"
  )
  res_df$Marker_associated_celltype <- celltype_name
  res_df$multinomial_p_val <- p_val
  res_df$dataset_name <- "PBMC 2 batch base balanced"
  res_df <- res_df[order(res_df$Top_celltype_downsampled_coefficient), ]
  print(res_df)
  return(res_df)
}

# Iterate over all celltypes (no mixing) and get multinomial test results
## THIS WILL NOT CONSIDER MARKER GENES THAT ARE MARKERS TO MORE THAN
## ONE CELLTYPE 
celltypes_list <- as.list(unique(base_marker_genes$Celltype))
multinomial_results <- lapply(celltypes_list, function(x) {
  results <- marker_celltype_multinomial(x, marker_lm_summaries_concat)
})

# Concatenate, fdr correct, and save the results
multinomial_results_concat <- Reduce(rbind, multinomial_results)
multinomial_results_concat$multinomial_q_val <- p.adjust(
  multinomial_results_concat$multinomial_p_val
)
fwrite(
  multinomial_results_concat,
  paste0(
    "outs/control/results/",
    "09_pbmc_base_ds_dge_marker_celltypes_top_ds_celltype_",
    "coeff_multinom_tests_low_rep_no_liger.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)