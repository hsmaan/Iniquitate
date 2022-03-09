library(data.table)
library(tidyverse)
library(reshape2)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Change to results dir for control data 
setwd("../../../results/control/")

# Set seed for any sampling done 
set.seed(42)

###Statistical analysis of PBMC 2 batch balanced data - baseline ### 

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

# Change to top level dir 
setwd("../../..")

### Statistical tests for ARI/AMI/etc. ### 

# Merge imbalance and clustering summary results
imba_clus_merged <- merge(
  clus_concat,
  imba_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate"
  )
)

# Format celltype names 
imba_clus_merged$`Downsampled celltypes` <- plyr::mapvalues(
  imba_clus_merged$`Downsampled celltypes`,
  from = c(
    "Monocyte_CD14",
    "Monocyte_FCGR3A",
    "CD4 T cell",
    "CD8 T cell"
  ),
  to = c(
    "CD14+ Monocyte",
    "FCGR3A+ Monocyte",
    "CD4+ T cell",
    "CD8+ T cell"
  )
)

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

# Do a wilcoxin rank sum test between all control, ablated, and downsampled
# experiment results - celltype ARI
imba_clus_control_celltype_aris <- imba_clus_merged[
  imba_clus_merged$type %in% c("Control")
]$`Celltype ARI Imbalanced`
imba_clus_downsampled_celltype_aris <- imba_clus_merged[
  imba_clus_merged$type %in% c("Downsampled")
]$`Celltype ARI Imbalanced`
imba_clus_ablated_celltype_aris <- imba_clus_merged[
  imba_clus_merged$type %in% c("Ablated")
]$`Celltype ARI Imbalanced`

# Downsample control to same number as ablated and downsampled for pairwise
# comparisons 
imba_clus_control_celltype_aris_sampled <- sample(
  imba_clus_control_celltype_aris,
  size = length(imba_clus_downsampled_celltype_aris)
)

# First test - compare between pooled downsampled and ablated to all controls
c_ari_ctrl_ds_abla_pooled_wrs <- wilcox.test(
  imba_clus_control_celltype_aris,
  c(
    imba_clus_ablated_celltype_aris,
    imba_clus_downsampled_celltype_aris
  )
)$p.value

# Second test - compare downsampled to all controls
c_ari_ctrl_ds_wrs <- wilcox.test(
  imba_clus_control_celltype_aris,
  imba_clus_downsampled_celltype_aris
)$p.value

# Third test - compare ablated to all controls 
c_ari_abla_ds_wrs <- wilcox.test(
  imba_clus_control_celltype_aris,
  imba_clus_ablated_celltype_aris
)$p.value

# Create and save summary dataframe of statistics 
c_ari_overall_results <- data.frame(
  "Test" = c(
    "Ctrl vs pooled ds and ablated",
    "Ctrl vs ds",
    "Ctrl vs ablated"
  ),
  "Metric" = "Celltype ARI",
  "Type" = "Global",
  "Wilcoxin_Rank_Sum_p" = c(
    c_ari_ctrl_ds_abla_pooled_wrs,
    c_ari_ctrl_ds_wrs,
    c_ari_abla_ds_wrs
  )
)
fwrite(
  c_ari_overall_results,
  "outs/control/results/06_celltype_ari_global_tests.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Do a wilcoxin rank sum test between all control, ablated, and downsampled
# experiment results - batch ARI
imba_clus_control_batch_aris <- imba_clus_merged[
  imba_clus_merged$type %in% c("Control")
]$`Batch ARI`
imba_clus_downsampled_batch_aris <- imba_clus_merged[
  imba_clus_merged$type %in% c("Downsampled")
]$`Batch ARI`
imba_clus_ablated_batch_aris <- imba_clus_merged[
  imba_clus_merged$type %in% c("Ablated")
]$`Batch ARI`

# First test - compare between pooled downsampled and ablated to all controls
b_ari_ctrl_ds_abla_pooled_wrs <- wilcox.test(
  imba_clus_control_batch_aris,
  c(
    imba_clus_ablated_batch_aris,
    imba_clus_downsampled_batch_aris
  )
)$p.value

# Second test - compare downsampled to all controls
b_ari_ctrl_ds_wrs <- wilcox.test(
  imba_clus_control_batch_aris,
  imba_clus_downsampled_batch_aris
)$p.value

# Third test - compare ablated to all controls 
b_ari_abla_ds_wrs <- wilcox.test(
  imba_clus_control_batch_aris,
  imba_clus_ablated_batch_aris
)$p.value

# Create and save summary dataframe of statistics 
b_ari_overall_results <- data.frame(
  "Test" = c(
    "Ctrl vs pooled ds and ablated",
    "Ctrl vs ds",
    "Ctrl vs ablated"
  ),
  "Metric" = "Batch ARI",
  "Type" = "Global",
  "Wilcoxin_Rank_Sum_p" = c(
    b_ari_ctrl_ds_abla_pooled_wrs,
    b_ari_ctrl_ds_wrs,
    b_ari_abla_ds_wrs
  )
)
fwrite(
  b_ari_overall_results,
  "outs/control/results/06_batch_ari_global_tests.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Perform the above for celltype ARI through an ANOVA test, after taking 
# into account the method and celltype that was downsampled 

# Subset for and rename relevant columns in imba_clus_merged
imba_clus_merged_sub_cari <- imba_clus_merged[
  ,c(
    "Celltype ARI Imbalanced",
    "Method",
    "Downsampled celltypes",
    "type"
  )
]
colnames(imba_clus_merged_sub_cari) <- c(
  "celltype_ari_imbalanced",
  "method",
  "downsampled_celltypes",
  "type"
)

# Perform ANOVA test for celltype ari and type, taking into account the 
# celltype that was downsampled and the method used - save result
celltype_ari_model <- lm(
  as.formula(
    paste0(
      "celltype_ari_imbalanced", 
      "~",
      "method+",
      "downsampled_celltypes+",
      "type"
    )
  ),
  data = imba_clus_merged_sub_cari
)
cari_anova_result <- anova(celltype_ari_model, test = "F")
fwrite(
  cari_anova_result,
  "outs/control/results/06_celltype_ari_global_anova.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Perform the same as above for batch ARI through an ANOVA test, after taking 
# into account the method and celltype that was downsampled 

# Subset for and rename relevant columns in imba_clus_merged
imba_clus_merged_sub_bari <- imba_clus_merged[
  ,c(
    "Batch ARI",
    "Method",
    "Downsampled celltypes",
    "type"
  )
]
colnames(imba_clus_merged_sub_bari) <- c(
  "batch_ari_imbalanced",
  "method",
  "downsampled_celltypes",
  "type"
)

# Perform ANOVA test for batch ari and type, taking into account the 
# celltype that was downsampled and the method used - save result
batch_ari_model <- lm(
  as.formula(
    paste0(
      "batch_ari_imbalanced", 
      "~",
      "method+",
      "downsampled_celltypes+",
      "type"
    )
  ),
  data = imba_clus_merged_sub_bari
)
bari_anova_result <- anova(batch_ari_model, test = "F")
fwrite(
  bari_anova_result,
  "outs/control/results/06_batch_ari_global_anova.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

