# Pre-processing

## Load packages
library(Seurat)
library(here)
library(glue)


## Define parameters
source(here("scRNA-seq/bin/utils_final_clusters.R"))


## Load data
seurat <- readRDS(here("scRNA-seq/results/R_objects/final_clusters/20230124/20230124_GCBC_seurat_obj.rds"))


# Cell cycle regression
seurat <- seurat[, seurat$assay == "3P"]
gcbc_levels_old <- c("DZ early Sphase", "DZ late Sphase", "DZ early G2Mphase",
                     "DZ late G2Mphase", "DZ cell cycle exit", "DZ non proliferative",
                     "DZ_LZ transition", "LZ", "LZ_DZ reentry commitment", "LZ proliferative",
                     "LZ_DZ transition", "Precursor MBCs", "Reactivated proliferative MBCs",
                     "prePC")
gcbc_levels_new <- c("DZ early Sphase", "DZ late Sphase", "DZ early G2Mphase",
                     "DZ late G2Mphase", "DZ cell cycle exit", "DZ non proliferative",
                     "DZ-LZ transition", "LZ", "LZ-DZ reentry commitment", "LZ proliferative",
                     "LZ_DZ transition", "Precursor MBCs", "Reactivated proliferative MBCs",
                     "PC Commited Light Zone GCBC")
seurat$annotation_20230124 <- factor(
  seurat$names_level_5_clusters_eta,
  levels = gcbc_levels_old
)
levels(seurat$annotation_20230124) <- gcbc_levels_new
seurat <- FindVariableFeatures(seurat, nfeatures = 2500)
seurat <- ScaleData(
  seurat,
  vars.to.regress = c("S.Score", "G2M.Score"),
  features = rownames(seurat)
)


# Save
saveRDS(seurat, here("scRNA-seq/results/R_objects/revision/seurat_gcbc_cell_cycle_regressed.rds"))

