# This script runs inferCNV (https://github.com/broadinstitute/inferCNV/)


# Load packages
print("Loading packages...")
library(Seurat)
library(here)
library(glue)
library(tidyverse)
library(infercnv)
set.seed(1234)


# Parse command line arguments
print("Parsing command-line arguments...")
args <- commandArgs(trailingOnly = TRUE)
donor_id <- args[[1]]


# Define parameters
# path_to_dir <- "~/Desktop/CNAG/mnt_clust/RICHTER/current"
print("Defining parameters...")
path_to_dir <- here()
path_to_utils <- glue("{path_to_dir}/scRNA-seq/bin/utils.R")
path_to_tumoral <- glue("{path_to_dir}/MCL/results/R_objects/7.seurat_tumoral_{donor_id}_clustered.rds")
path_to_microenv <- glue("{path_to_dir}/MCL/results/R_objects/6.seurat_microenvironment_{donor_id}_rna_only.rds")
path_to_gene_order_f <- glue("{path_to_dir}/MCL/5-infercnv/tmp/gencode_v21_gen_pos.complete.txt")
path_to_gene_order_donor <- glue("{path_to_dir}/MCL/5-infercnv/tmp/gencode_v21_gen_pos.complete.txt")
path_to_cell_annot <- glue("{path_to_dir}/MCL/5-infercnv/tmp/cell_annotation_{donor_id}.txt")
path_to_outdir <- glue("{path_to_dir}/MCL/results/inferCNV/donor_{donor_id}")
path_to_save <- glue("{path_to_dir}/MCL/results/inferCNV/infercnv_obj_{donor_id}.rds")
dir.create(path_to_outdir, recursive = TRUE)
cutoff_infercnv <- 0.1


# Source functions
print("Sourcing functions...")
source(path_to_utils)


# Load data
print("Loading data...")
seurat_tumoral <- readRDS(path_to_tumoral)
seurat_microenv <- readRDS(path_to_microenv)
gene_order_file <- read_tsv(
  path_to_gene_order_f,
  col_names = c("gene", "chromosome", "start", "end")
)


# Downsample objectsand merge
seurat_tumoral <- subset(seurat_tumoral, downsample = 350)
selected_cells <- sample(colnames(seurat_microenv), size = 3000, replace = FALSE)
seurat_microenv <- subset(seurat_microenv, cells = selected_cells)
seurat_tumoral$annotation_infercnv <- str_c(
  "malignant",
  seurat_tumoral$seurat_clusters,
  sep = "_"
)
seurat_microenv$annotation_infercnv <- "normal"
seurat <- merge(x = seurat_tumoral, y = seurat_microenv)
rm(seurat_microenv, seurat_tumoral)
gc()


# Write cell annotations file
print("Writing cell annotations file...")
seurat$cell_barcodes <- colnames(seurat)
cell_annotations <- seurat@meta.data[, c("cell_barcodes", "annotation_infercnv")]
write_tsv(cell_annotations, path_to_cell_annot, col_names = FALSE)


# Obtain gene ordering file
print("Preparing and saving gene annotation file...")
# gene_order_file <- as.data.frame(gene_order_file)
# gene_order_file$gene <- str_remove(gene_order_file$gene, "\\|ENSG.*$")
# gene_order_file <- gene_order_file[gene_order_file$gene %in% rownames(seurat), ]
# gene_order_file <- gene_order_file[!duplicated(gene_order_file$gene), ]
# seurat <- subset(seurat, features = gene_order_file$gene)
# write_tsv(gene_order_file, path_to_gene_order_donor, col_names = FALSE)
gene_order_file <- as.data.frame(gene_order_file)
gene_order_file$gene <- str_remove(gene_order_file$gene, "\\|ENSG.*$")
seurat <- subset(
  seurat,
  features = rownames(seurat)[rownames(seurat) %in% gene_order_file$gene]
)
gene_order_file <- gene_order_file[match(rownames(seurat), gene_order_file$gene), ]
write_tsv(gene_order_file, path_to_gene_order_donor, col_names = FALSE)


# Create inferCNV object
print("Creating inferCNV object...")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = seurat[["RNA"]]@counts,
  annotations_file = path_to_cell_annot,
  delim = "\t",
  gene_order_file = path_to_gene_order_donor,
  ref_group_names = "normal",
  chr_exclude = "chrM"
)


# Run inferCNV
print("Running inferCNV...")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = cutoff_infercnv,
  out_dir = path_to_outdir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 8,
  save_rds = FALSE,
  save_final_rds = TRUE,
  plot_steps = FALSE
)


# Save
print("Saving...")
saveRDS(infercnv_obj, path_to_save)
