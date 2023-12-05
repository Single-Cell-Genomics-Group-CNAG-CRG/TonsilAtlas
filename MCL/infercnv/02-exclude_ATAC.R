# This script removes the ATAC slot, so we don't need to install signac
# in the infercnv conda environment


library(Seurat)
library(Signac)
library(here)
seurat_102 <- readRDS(here("MCL/results/R_objects/6.seurat_microenvironment_102.rds"))
seurat_102@assays$ATAC <- NULL
seurat_microenv@reductions$atacUMAP <- NULL
seurat_microenv@reductions$lsi <- NULL
saveRDS(seurat_102, here("MCL/results/R_objects/6.seurat_microenvironment_102_rna_only.rds"))


seurat_413 <- readRDS(here("MCL/results/R_objects/6.seurat_microenvironment_413.rds"))
seurat_413@assays$ATAC <- NULL
seurat_413@reductions$atacUMAP <- NULL
seurat_413@reductions$lsi <- NULL
saveRDS(seurat_413, here("MCL/results/R_objects/6.seurat_microenvironment_413_rna_only.rds"))