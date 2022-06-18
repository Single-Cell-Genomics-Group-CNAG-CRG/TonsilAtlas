# Variables
ver <- "2020-09-22"

## PATHS and common variables
version_dir <- glue::glue("{ver}")

## Create directory for plots
plt_dir <- glue::glue("{ver}/plots_{ver}") 

## Create directory for RDS objects
robj_dir <- glue::glue("{ver}/R_objects_{ver}")

## Paths to all the folders
qc <- "02-QC"
clust <- "03-clustering"
anot <-  "04-annotation"
decon <- "05-sc_map"
regulon <- "06-regulon_integration"
spatialde <- "spatially_variable_genes"
misty <- "MISTy"
cd4 <- "CD4-Analysis"
cd8 <- "CD8-Analysis"
gct <- "GC-Analysis"
stereo <- "stereoseq"
myeloid <- "myeloid_integration"
epithelium <- "epithelium_integration"
plasma <- "plasma_integration"
fibro <- "fibroblasts"
chemo <- "chemokine-analysis"
cc_interact <- "cell-cell-interaction"
fig1 <- "Spatial-Figure-1"
fig2 <- "Spatial-Figure-2"
fig3 <- "Spatial-Figure-3"
fig4 <- "Spatial-Figure-4"
fig5 <- "Spatial-Figure-5"
fig_s <- "Supplementary-Figures"

## Paths to data
spaceranger <- "01-spaceranger/results/"
filtered_mtrx <- "outs/filtered_feature_bc_matrix/"

## QC dictionary
qc_dict <- list()
qc_dict[["c28w2r_7jne4i"]] <- list()
qc_dict[["c28w2r_7jne4i"]][["min_lib_size"]] <- 500
qc_dict[["c28w2r_7jne4i"]][["min_n_genes"]] <- 50
qc_dict[["c28w2r_7jne4i"]][["max_pct_mt"]] <- 50
qc_dict[["c28w2r_7jne4i"]][["min_cells"]] <- 8

qc_dict[["esvq52_nluss5"]] <- list()
qc_dict[["esvq52_nluss5"]][["min_lib_size"]] <- 500
qc_dict[["esvq52_nluss5"]][["min_n_genes"]] <- 50
qc_dict[["esvq52_nluss5"]][["max_pct_mt"]] <- 50
qc_dict[["esvq52_nluss5"]][["min_cells"]] <- 8

qc_dict[["p7hv1g_tjgmyj"]] <- list()
qc_dict[["p7hv1g_tjgmyj"]][["min_lib_size"]] <- 500
qc_dict[["p7hv1g_tjgmyj"]][["min_n_genes"]] <- 50
qc_dict[["p7hv1g_tjgmyj"]][["max_pct_mt"]] <- 50
qc_dict[["p7hv1g_tjgmyj"]][["min_cells"]] <- 8

qc_dict[["tarwe1_xott6q"]] <- list()
qc_dict[["tarwe1_xott6q"]][["min_lib_size"]] <- 500
qc_dict[["tarwe1_xott6q"]][["min_n_genes"]] <- 50
qc_dict[["tarwe1_xott6q"]][["max_pct_mt"]] <- 50
qc_dict[["tarwe1_xott6q"]][["min_cells"]] <- 8


## Sample ID data
id_sp_df <- readr::read_csv(file = here::here("01-spaceranger/data/sample_id.txt"))
