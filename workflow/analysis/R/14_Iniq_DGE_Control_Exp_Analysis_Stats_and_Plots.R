library(data.table)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggExtra)
library(ggpubr)
library(dotwhisker)
library(Seurat)
library(SeuratDisk)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)
library(networkD3)

# Helper functions
`%ni%` <- Negate(`%in%`)

# Kevin's palette for plotting many catagoricals 
kev_palette <- c(
  "dodgerblue2", "#E31A1C",
  "green4",
  "#6A3D9A", 
  "#FF7F00", 
  "black", "gold1",
  "skyblue2", "#FB9A99", 
  "palegreen2",
  "#CAB2D6", 
  "#FDBF6F", 
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

# Change to results dir for the control dge experiment
setwd("../../../results/control_dge_control/")

##### Data from of PBMC 2 batch balanced data - controlled setting #####

# In this setting, the first batch is randomly downsampled and 
# the DGE experiment is performed based on a control setup using the 
# non-downsampled batch as the reference. 

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("dge_ranking_stats_marker_sub_control")
dge_rank_control_files <- list.files()
dge_rank_control_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_rank_control_files,
  value = TRUE
)
dge_rank_control_loaded <- lapply(dge_rank_control_files, fread)
dge_rank_control_concat <- Reduce(rbind, dge_rank_control_loaded)
gc()

# Load in and concatenate imbalance summary files 
setwd("../control_imbalance_summaries/")
imba_control_files <- list.files()
imba_control_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_control_files,
  value = TRUE
)
imba_control_loaded <- lapply(imba_control_files, fread)
imba_control_concat <- Reduce(rbind, imba_control_loaded)

#### Data from PBMC 2 batch balanced perturbation experiments - power control
####

# In this setting, only 200 cells per cell-type in each batch are utilized 
# in the perturbation experiments, to control for power - otherwise the 
# integrated space for the initial PBMC perturbation experiments would 
# have twice the number of cells as the control experiments

# Change to results dir for the power control pbmc dge experiment
setwd("../../control_power/")

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("dge_ranking_stats_marker_sub")
dge_rank_power_files <- list.files()
dge_rank_power_files <- grep(
  "pbmc_2_batch_base_balanced",
  dge_rank_power_files,
  value = TRUE
)
dge_rank_power_loaded <- lapply(dge_rank_power_files, fread)
dge_rank_power_concat <- Reduce(rbind, dge_rank_power_loaded)
gc()

# Load in and concatenate imbalance summary files 
setwd("../imbalance_summaries/")
imba_power_files <- list.files()
imba_power_files <- grep(
  "pbmc_2_batch_base_balanced",
  imba_power_files,
  value = TRUE
)
imba_power_loaded <- lapply(imba_power_files, fread)
imba_power_concat <- Reduce(rbind, imba_power_loaded)

# Change to top level dir 
setwd("../../..")  

# Create output directory for figures
if (!dir.exists("outs/dge_control_exp/figures")) {
  dir.create("outs/dge_control_exp/figures", recursive = TRUE)
}

# Create output directory for results
if (!dir.exists("outs/dge_control_exp/results")) {
  dir.create("outs/dge_control_exp/results", recursive = TRUE)
}

# First merge together imba and dge dataframes for both the control and power
# subsets
imba_dge_control_merged <- merge(
  imba_control_concat,
  dge_rank_control_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
)
imba_dge_control_merged <- distinct(imba_dge_control_merged)

imba_dge_power_merged <- merge(
  imba_power_concat,
  dge_rank_power_concat,
  by = c(
    "Number of batches downsampled",
    "Number of celltypes downsampled",
    "Proportion downsampled",
    "Replicate",
    "Dataset"
  )
)
imba_dge_power_merged <- distinct(imba_dge_power_merged)

# Remove marker genes not present in both experiments
control_markers <- unique(imba_dge_control_merged$Gene)
power_markers <- unique(imba_dge_power_merged$Gene)
marker_int <- intersect(control_markers, power_markers)
imba_dge_control_merged <- imba_dge_control_merged[
  imba_dge_control_merged$Gene %in% marker_int
]
imba_dge_power_merged <- imba_dge_power_merged[
  imba_dge_power_merged$Gene %in% marker_int
]

# Add an indicator to both dataframes - 'integrated vs unintegrated'
imba_dge_control_merged$int <- "No"
imba_dge_power_merged$int <- "Yes"

# Drop unecessary columns from the power control
imba_dge_power_merged <- imba_dge_power_merged[
  ,-c(
    "Total batches",
    "Batches downsampled",
    "Celltype intersection ratio",
    "Mean proportion cosine distance",
    "Length coeff var"
  )
]

# Put in placeholder 'method' column for the control experiment
imba_dge_control_merged$Method <- "Standard"

# Merge together both the dataframes
imba_dge_merged <- rbind(
  imba_dge_control_merged,
  imba_dge_power_merged
)

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
imba_dge_merged$type <- factor(imba_dge_merged$type)

# Convert int to factor
imba_dge_merged$int <- factor(imba_dge_merged$int)

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
## integrated/unintegrated, after factoring in marker gene, celltype 
## downsampled, and type (control/ds/ablate) - format and save the result

# Note that method is not factored in here because it is collinear with 
# 'int' 
dge_model_fit <- lm(
  as.formula(
    paste0(
      "max_rank", 
      "~",
      "gene+",
      "downsampled_celltypes+",
      "type+",
      "int"
    )
  ),
  data = imba_dge_merged
)
dge_anova_result <- anova(dge_model_fit, test = "F")
dge_anova_result_dt <- as.data.table(dge_anova_result, keep.rownames = TRUE)
colnames(dge_anova_result_dt)[1] <- "Covariate"
dge_anova_result_dt$dataset_name <- "PBMC 2 batch base balanced"
dge_anova_result_dt$metric <- "Max DGE rank marker gene"
dge_anova_result_dt$last_covariate <- "int (integrated or not)"
fwrite(
  dge_anova_result_dt,
  paste0(
    "outs/dge_control_exp/results/",
    "14_pbmc_base_control_exp_dge_rank_aov_results.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Plot the ANOVA F-values for a supplementary figure 
dge_anova_result_dt_nores <- dge_anova_result_dt[
  dge_anova_result_dt$Covariate != "Residuals"
]

f_vals <- dge_anova_result_dt_nores$`F value`
covars <- dge_anova_result_dt_nores$Covariate

dge_aov_comp_df <- data.frame(
  "Covariates" = covars,
  "F_values" = f_vals
)

dge_aov_comp_df$Covariates <- plyr::mapvalues(
  dge_aov_comp_df$Covariates,
  from = c(
    "gene",
    "downsampled_celltypes",
    "type",
    "int"
  ),
  to = c(
    "Marker gene",
    "Cell-type downsampled",
    "Balanced vs perturbed",
    "Integrated vs control"
  )
)
dge_aov_comp_df$Covariates <- factor(
  dge_aov_comp_df$Covariates,
  levels = rev(
      c(
      "Marker gene",
      "Cell-type downsampled",
      "Balanced vs perturbed",
      "Integrated vs control"
    )
  )
)

ggplot(data = dge_aov_comp_df, aes(Covariates, F_values)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(),
    fill = "#E41A1C"
  ) + 
  scale_fill_brewer(palette = "Set1") +
  theme_classic() + 
  coord_flip () +
  labs(x = "Covariate", y = "ANOVA F-statistic") +
  theme(axis.title.x = element_text(size = 16)) +
  theme(axis.title.y = element_text(size = 16)) +
  theme(strip.text.x = element_text(size = 16)) +
  theme(strip.text.y = element_text(size = 16)) +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 16)) +
  theme(legend.title = element_text(size = 16)) +
  theme(legend.text = element_text(size = 16)) +
  theme(aspect.ratio = 1)
ggsave(
  paste0(
    "outs/dge_control_exp/figures/",
    "14_pbmc_base_control_exp_dge_rank_aov_f_value_results.pdf"
  ),
  width = 12,
  height = 12,
  device = cairo_pdf
)
