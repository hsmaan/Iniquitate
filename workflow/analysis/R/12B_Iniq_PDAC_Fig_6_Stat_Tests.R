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

### Note that this analysis done without LIGER ### 
`%ni%` <- Negate(`%in%`)
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

# Load in and concatenate imbalance summary files 
setwd("../../../results/pdac_comp/imbalance_summaries/")
imba_files <- list.files()
imba_loaded <- lapply(imba_files, fread)
imba_concat <- Reduce(rbind, imba_loaded)
gc()

# Load in and concatenate celltype summary files
setwd("../celltype_imbalance_summaries")
cimba_files <- list.files()
cimba_loaded <- lapply(cimba_files, fread)
cimba_concat <- Reduce(rbind, cimba_loaded)
gc()

# Load in and concatenate the clustering summary results 
setwd("../clustering_summaries/")
clus_files <- list.files()
clus_loaded <- lapply(clus_files, fread)
clus_concat <- Reduce(rbind, clus_loaded)
clus_concat <- clus_concat[clus_concat$Method != "liger"]
gc()

# Load in and concatenate clustering concordance summaries 
setwd("../clustering_concord_summaries/")
clus_concord_files <- list.files()
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
dge_loaded <- lapply(dge_files, fread)
dge_concat <- Reduce(rbind, dge_loaded)
dge_concat <- dge_concat[dge_concat$`Method 1` != "liger"]
dge_concat <- dge_concat[dge_concat$`Method 2` != "liger"]
gc()

# Load in and concatenate dge ranking summaries, subset by marker genes
setwd("../dge_ranking_stats_marker_sub")
dge_rank_files <- list.files()
dge_rank_loaded <- lapply(dge_rank_files, fread)
dge_rank_concat <- Reduce(rbind, dge_rank_loaded)
dge_rank_concat <- dge_rank_concat[dge_rank_concat$Method != "liger"]
gc()

# Load in markers and corresponding celltypes 
setwd("../marker_results/")
base_marker_genes <- fread(
  "peng_pdac_tumor_annot_8_batch_preintegration_marker_selection.tsv"
) 

# Load in and concatenate knn classification summaries
setwd("../knn_classification_reports/")
knn_files <- list.files()
knn_loaded <- lapply(knn_files, fread)
knn_concat <- Reduce(rbind, knn_loaded)
knn_concat <- knn_concat[knn_concat$Method != "liger"]
gc()

# Change to top level dir 
setwd("../../..")

# Make pdac comp output dir if it doesn't exist
if (!dir.exists("outs/pdac_comp/figures")) {
  dir.create("outs/pdac_comp/figures", recursive = TRUE)
}
if (!dir.exists("outs/pdac_comp/results")) {
  dir.create("outs/pdac_comp/results", recursive = TRUE)
}

### Statistical test -  downsampling results on KNN classification scores
### of given subsets/compartments 
imba_knn_merged <- merge(
imba_concat,
knn_concat,
by = c(
  "Number of batches downsampled",
  "Number of celltypes downsampled",
  "Proportion downsampled",
  "Replicate"
)
)
imba_knn_merged <- distinct(imba_knn_merged)

# Subset for only cases where the celltype/compartment downsampled is equal to 
# the celltype being classified
imba_knn_merged_celltype <- imba_knn_merged[
  imba_knn_merged$Celltype == imba_knn_merged$`Downsampled celltypes` |
    imba_knn_merged$`Downsampled celltypes` %in% c("None")
]

# Indicate which panels are control and which ones are ablations or downsampling
imba_knn_merged_celltype$type <- ifelse(
  imba_knn_merged_celltype$`Number of batches downsampled` == 0,
  "Control",
  ifelse(
    imba_knn_merged_celltype$`Proportion downsampled` == 0,
    "Ablated",
    "Downsampled"
  )
) 

# Indicate the separate compartments
compartments <- unique(imba_knn_merged_celltype$Celltype)

# Create a function to do an ANOVA test for the F1 score based on each
# compartment utilized 
anova_compart_knn <- function(
    compartment, 
    dataset
){
  # Subset data for the given compartment
  dataset_sub <- dataset[dataset$Celltype == compartment]
  
  # Format the data columns for lm 
  colnames(dataset_sub) <- plyr::mapvalues(
    colnames(dataset_sub),
    from = c(
      "F1-score",
      "Method",
      "type"
    ),
    to = c(
      "f1_score",
      "method",
      "type"
    )
  )
  
  # Fit ANOVA model
  model_fit <- lm(
    as.formula(
      paste0(
        "f1_score", 
        "~",
        "method+",
        "type"
      )
    ),
    data = dataset_sub
  )
  anova_result <- anova(model_fit, test = "F")
  
  # Format results and return
  anova_result_dt <- as.data.table(anova_result, keep.rownames = TRUE)
  colnames(anova_result_dt)[1] <- "Covariate"
  anova_result_dt$compartment_name <- compartment
  anova_result_dt$metric <- "F1 score"
  anova_result_dt$last_covariate <- "type"
  return(anova_result_dt)
}

# Iterate over compartments and get the significance of ds/ablation 
knn_anova_comp_results <- mapply(
  anova_compart_knn,
  compartment = compartments,
  MoreArgs = list(
    dataset = imba_knn_merged_celltype
  ),
  SIMPLIFY = FALSE
)
knn_anova_comp_results

# Save the concatenated results and plot the ANOVA F-values for 
# a supplementary figure 
knn_anova_comp_results_concat <- Reduce(rbind, knn_anova_comp_results)
fwrite(
  knn_anova_comp_results_concat,
  "outs/pdac_comp/results/12B_comp_specific_ds_knn_f1_score_anovas.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

knn_anova_comp_results_concat_nores <- knn_anova_comp_results_concat[
  knn_anova_comp_results_concat$Covariate != "Residuals"
]

f_vals <- knn_anova_comp_results_concat_nores$`F value`
covars <- knn_anova_comp_results_concat_nores$Covariate
comps <- knn_anova_comp_results_concat_nores$compartment_name

knn_aov_comp_df <- data.frame(
  "Covariates" = covars,
  "F_values" = f_vals,
  "Compartment" = comps
)

knn_aov_comp_df_melted <- reshape2::melt(
  knn_aov_comp_df,
  id.vars = c("Compartment", "Covariates"),
  measure.vars = "F_values"
)
knn_aov_comp_df_melted$Covariates <- plyr::mapvalues(
  knn_aov_comp_df_melted$Covariates,
  from = c(
    "type",
    "method"
  ),
  to = c(
    "Unperturbed vs perturbed",
    "Integration method"
  )
) 

ggplot(data = knn_aov_comp_df_melted, aes(Covariates, value)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(),
    aes(
      fill = Compartment
    )
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
    "outs/pdac_comp/figures/",
    "12B_pdac_knn_aov_comp_ds_f_statistic.pdf"
  ),
  width = 12,
  height = 12,
  device = cairo_pdf
)
