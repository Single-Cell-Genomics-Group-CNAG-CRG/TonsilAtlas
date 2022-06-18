# For each cluster, this script performs differential expression analysis (DEA)
# between the following conditions
# 1. Kid vs young adult
# 2. Kid vs old adult
# 3. Young adult vs old adult
# 4. Male vs female


# Load packages
library(Seurat)
library(tidyverse)


# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
cell_type <- args[[1]]


# Define paths
path_to_level_3 <- here::here("scRNA-seq/results/R_objects/level_3/")
path_to_level_3_cell_type <- str_c(path_to_level_3, cell_type, sep = "")
path_to_obj <- str_c(
  path_to_level_3_cell_type,
  "/",
  cell_type,
  "_clustered_level_3.rds",
  sep = ""
)
path_to_annotation <- str_c(here::here("annotation/level_3/"), cell_type, sep = "")
path_to_save_kid_vs_young <- str_c(path_to_annotation, "/", cell_type, "_DEA_kid_vs_young.xlsx")
path_to_save_kid_vs_old <- str_c(path_to_annotation, "/", cell_type, "_DEA_kid_vs_old.xlsx")
path_to_save_young_vs_old <- str_c(path_to_annotation, "/", cell_type, "_DEA_young_vs_old.xlsx")
path_to_save_male_vs_female <- str_c(path_to_annotation, "/", cell_type, "_DEA_male_vs_female.xlsx")


# Define function
perform_dea <- function(seurat_obj, idents, ident.1, ident.2) {
  x <- seurat_obj[[idents, drop = TRUE]]
  if (length(x[x == ident.1]) < 4 | length(x[x == ident.2]) < 4) {
    return(NULL)
  }
  Idents(seurat_obj) <- idents
  if (ident.1 %in% Idents(seurat_obj) & ident.2 %in% Idents(seurat_obj)) {
    df <- FindMarkers(
      seurat_obj,
      ident.1 = ident.1,
      ident.2 = ident.2,
      only.pos = FALSE,
      logfc.threshold = 0.3,
      verbose = TRUE
    )
    df <- df %>%
      rownames_to_column(var = "gene") %>%
      filter(p_val_adj < 0.001) %>%
      arrange(desc(avg_logFC))
    df 
  } else {
    return(NULL)
  }
}


# Load data
seurat <- readRDS(path_to_obj)


# AGE
## Kid vs young adult
seurat_list <- SplitObject(seurat, split.by = "seurat_clusters")
markers_kid_vs_young_adult <- purrr::map(
  seurat_list,
  perform_dea,
  idents = "age_group",
  ident.1 = "kid",
  ident.2 = "young_adult"
)
sorted_names <- as.character(sort(as.numeric(names(seurat_list))))
markers_kid_vs_young_adult <- markers_kid_vs_young_adult[sorted_names]


## Kid vs old adult
markers_kid_vs_old_adult <- purrr::map(
  seurat_list,
  perform_dea,
  idents = "age_group",
  ident.1 = "kid",
  ident.2 = "old_adult"
)
markers_kid_vs_old_adult <- markers_kid_vs_old_adult[sorted_names]


## Young adult vs old adult
markers_young_adult_vs_old_adult <- purrr::map(
  seurat_list,
  perform_dea,
  idents = "age_group",
  ident.1 = "young_adult",
  ident.2 = "old_adult"
)
markers_young_adult_vs_old_adult <- markers_young_adult_vs_old_adult[sorted_names]


# SEX
## Male vs female
markers_male_vs_female <- purrr::map(
  seurat_list,
  perform_dea,
  idents = "sex",
  ident.1 = "male",
  ident.2 = "female"
)
markers_male_vs_female <- markers_male_vs_female[sorted_names]


# Save
openxlsx::write.xlsx(markers_kid_vs_young_adult, path_to_save_kid_vs_young)
openxlsx::write.xlsx(markers_kid_vs_old_adult, path_to_save_kid_vs_old)
openxlsx::write.xlsx(markers_young_adult_vs_old_adult, path_to_save_young_vs_old)
openxlsx::write.xlsx(markers_male_vs_female, path_to_save_male_vs_female)