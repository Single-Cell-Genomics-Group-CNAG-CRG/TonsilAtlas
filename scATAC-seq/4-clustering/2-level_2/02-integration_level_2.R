# This script performs feature selection, dimensionality reduction and
# batch effect correction for a specific cell type (level 2)


# Load packages
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(harmony)
library(tidyverse)
set.seed(666)

integration_level2 <- function(cell_type){
  
  print(cell_type)
  
  # Define parameters
  path_to_dir <- here::here("scATAC-seq/results/R_objects/level_2/")
  path_to_obj <- str_c(
    path_to_dir,
    cell_type,
    "/",
    cell_type,
    "_subsetted_level_2.rds",
    sep = ""
  )
  
  path_to_save <- str_c(
    here::here("scATAC-seq/results/R_objects/level_2/"),
    cell_type,
    "/",
    cell_type,
    "_integrated_level_2.rds",
    sep = ""
  )

    # Load data
    seurat <- readRDS(path_to_obj)
  
    # Normalization, dimensionality reduction 
    seurat <-  RunTFIDF(seurat) %>% 
    FindTopFeatures(min.cutoff = 10) %>% RunSVD()  
    
    # DepthCor(seurat)
    # seurat <- RunUMAP(object = seurat, reduction = 'lsi', dims = 2:40)
    # DimPlot(seurat)

    # Batch correction
    seurat <- RunHarmony(
      object = seurat,
      dims = 2:40,
      group.by.vars = 'assay',
      reduction = 'lsi',
      assay.use = 'peaks_macs',
      project.dim = FALSE,
      max.iter.harmony = 20
    )
    
    # Non-linear dimension reduction and clustering
    seurat <- RunUMAP(seurat, dims = 2:40, reduction = 'harmony') %>%
    FindNeighbors(reduction = "harmony", dims = 2:40)
    #DimPlot(object = seurat)
    
    print(seurat)
    table(seurat$annotation_level_1)
  
  # Save
  saveRDS(seurat, path_to_save)
}

cell_types = c("CD4_T", "NBC_MBC", "myeloid", 
              "GCBC", "Cytotoxic", "PC", "PDC",
              "FDC")

# Run the function per each cell type
lapply(cell_types, integration_level2)

