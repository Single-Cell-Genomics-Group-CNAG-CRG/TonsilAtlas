# This script plots the results obtained with inferCNV
# (https://github.com/broadinstitute/inferCNV/)


# Load packages
print("Loading packages...")
library(infercnv)
library(Seurat)
library(here)
library(glue)
set.seed(1234)


# Parse command line arguments
print("Parsing command-line arguments...")
args <- commandArgs(trailingOnly = TRUE)
donor_id <- args[[1]]


# Define paths
path_to_dir <- here()
path_to_infer_cnv_obj <- glue("{path_to_dir}/MCL/results/inferCNV/infercnv_obj_{donor_id}.rds")
path_to_tumoral <- glue("{path_to_dir}/MCL/results/R_objects/7.seurat_tumoral_{donor_id}_clustered.rds")
path_to_microenv <- glue("{path_to_dir}/MCL/results/R_objects/6.seurat_microenvironment_{donor_id}_rna_only.rds")
path_to_infercnv_out <- glue("{path_to_dir}/MCL/results/inferCNV/donor_{donor_id}")
path_to_save <-  glue("{path_to_dir}/MCL/results/inferCNV/{donor_id}")
# path_to_save <- "~/Desktop/diogenes/"


# Read data
infercnv_obj <- readRDS(path_to_infer_cnv_obj)
seurat_tumoral <- readRDS(path_to_tumoral)
seurat_microenv <- readRDS(path_to_microenv)


# Plot and save
plot_cnv(
  infercnv_obj = infercnv_obj,
  out_dir = path_to_save,
  title = "",
  obs_title = "",
  ref_title = "",
  cluster_references = FALSE,
  plot_chr_scale = TRUE,
  k_obs_groups = 1,
  custom_color_pal = color.palette(c("#f24153", "white", "#536db6"), between = c(2, 2)),
  color_safe_pal = FALSE,
  output_filename = glue("infercnv_custom_heatmap_{donor_id}"),
  output_format = "pdf",
  dynamic_resize = 0.25,
  write_expr_matrix = FALSE,
  useRaster = TRUE,
)


# Add to Seurat
seurat <- merge(x = seurat_tumoral, y = seurat_microenv)
seurat <- subset(seurat, cells = colnames(infercnv_obj@expr.data))
rm(seurat_microenv)
gc()
seurat <- infercnv::add_to_seurat(
  seurat_obj = seurat,
  infercnv_output_path = path_to_infercnv_out
)
metadata_df <- seurat@meta.data
rm(seurat)
gc()
seurat_tumoral$has_dupli_chrY <- NA
seurat_tumoral$has_loss_chrY <- NA
metadata_df <- metadata_df[rownames(metadata_df) %in% colnames(seurat_tumoral), ]
seurat_tumoral$has_dupli_chrY[rownames(metadata_df)] <- metadata_df$has_dupli_chrY
seurat_tumoral$has_loss_chrY[rownames(metadata_df)] <- metadata_df$has_loss_chrY
DimPlot(seurat_tumoral, group.by = "has_dupli_chrY", reduction = "umap")
DimPlot(seurat_tumoral, group.by = "has_loss_chrY", reduction = "umap")


# Save
saveRDS(
  seurat_tumoral,
  glue("{path_to_dir}/MCL/results/R_objects/8.seurat_tumoral_{donor_id}_clustered_infercnv.rds")
)


