# This script annotates NK and ILC cells


# Load packages
library(Seurat)
library(tidyverse)


# Read data
nk_ilc <- readRDS("~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/Cytotoxic/ILC_NK/ILC_NK_clustered_level_5.rds")


# Rename
nk_ilc$annotation_paper <- factor(nk_ilc$annotation_level_4)
levels(nk_ilc$annotation_paper) <- c(
  "NKp44+ ILC3",
  "ILC1",
  "CD16-CD56+ NK",
  "CD16-CD56- NK",
  "CD16+CD56- NK",
  "NKp44- ILC3"
)
ordered_levels <- c(
  "CD16-CD56+ NK",
  "CD16-CD56- NK",
  "CD16+CD56- NK",
  "ILC1",
  "NKp44+ ILC3",
  "NKp44- ILC3"
)
nk_ilc$annotation_paper <- factor(
  nk_ilc$annotation_paper,
  levels = ordered_levels
)
Idents(nk_ilc) <- "annotation_paper"


# Save
saveRDS(nk_ilc, "~/Desktop/CNAG/mnt_clust/tonsil_atlas/current/scRNA-seq/results/R_objects/level_5/Cytotoxic/ILC_NK/ILC_NK_annotated_level_5.rds")


