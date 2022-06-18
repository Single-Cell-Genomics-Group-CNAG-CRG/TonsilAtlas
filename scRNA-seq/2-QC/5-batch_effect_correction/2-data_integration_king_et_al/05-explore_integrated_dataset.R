# This script explores the information that we can leverage from King et al.


# Load packages
library(Seurat)
library(tidyverse)


# Load data
path_to_data <- here::here("scRNA-seq/results/R_objects/seurat_merged_with_king_et_al_integrated.rds")
tonsil <- readRDS(path_to_data)


# Plot key markers
selected_markers <- c("TOP2A", "VBREP1", "MIR155HG", "HSP90AB1", "NCL", "PRDM1",
                      "FCRL4")
feature_plots <- purrr::map(selected_markers, FeaturePlot, pt.size = 0.1)


# Subset
tonsil_sub <- subset(tonsil, subset = hospital == "Royal London")
DimPlot(tonsil_sub, group.by = "cell_type")



# Save
ggsave



