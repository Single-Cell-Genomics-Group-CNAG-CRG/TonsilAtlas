# Load packages
library(Seurat)
library(Signac)
library(tidyverse)
library(glue)

# Set options
options(timeout = 10000)


# Download cellranger outs for scATAC-seq (version 3) and unzip
url_cellranger <- "https://zenodo.org/records/11355186/files/scATAC-seq.zip?download=1"
destfile_cellranger <- "scATAC-seq.zip"
download.file(url_cellranger, destfile_cellranger)
unzip(destfile_cellranger)

# Download Seurat objects for scATAC-seq and unzip
url_seurat <- "https://zenodo.org/records/8373756/files/TonsilAtlasSeuratATAC.tar.gz?download=1"
destfile_seurat <- "TonsilAtlasSeuratATAC.tar.gz"
download.file(url_cellranger, destfile_seurat)
unzip(destfile_seurat)

# Read object and metadata
metadata <- read_csv("~/Desktop/zenodo_20240525/scATAC-seq/metadata/cellranger_atac_metadata.csv")
se <- readRDS("~/Downloads/scATAC-seq 3/20230911_tonsil_atlas_atac_seurat_obj.rds")


# Read barcodes for each gem_id and check that they match
gem_ids <- unique(se$gem_id[se$assay == "scATAC"])
all_tests <- map_lgl(gem_ids, \(x) {
  print(x)
  subproject <- metadata$subproject[metadata$gem_id == x]
  barcodes <- read_tsv(
    glue("~/Desktop/zenodo_20240525/scATAC-seq/{subproject}/{x}/filtered_peak_bc_matrix/barcodes.tsv"),
    col_names = "barcode"
  )
  barcodes$barcode <- glue("{x}_{barcodes$barcode}")
  out <- all(se$barcode[se$gem_id == x] %in% barcodes$barcode)
  out
})

# Compare
if (all(all_tests) == TRUE) {
  print("Cell barcodes in accessibility matrix and Seurat objects match!")
} else {
  print("Cell barcodes in accessibility matrix and Seurat DO NOT objects match!")  
}

