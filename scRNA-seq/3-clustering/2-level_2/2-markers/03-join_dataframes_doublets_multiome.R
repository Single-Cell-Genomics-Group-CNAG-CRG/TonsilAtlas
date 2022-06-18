# This script merges all the dataframes that contain the doublets that we exclude
# from the multiome dataset for each cell type (important to Pauli)


# Load packages
library(tidyverse)


# Load data
path_to_dfs <- here::here("scRNA-seq/3-clustering/2-level_2/tmp")
files_to_load <- list.files(
  path_to_dfs,
  pattern = "multiome_doublets_level_2.rds",
  full.names = TRUE
)
dfs <- purrr::map(files_to_load, readRDS)


# Merge
doublets_multiome_df <- bind_rows(dfs)
print(dim(doublets_multiome_df))
doublets_multiome_df <- doublets_multiome_df[doublets_multiome_df$assay == "multiome", ]
print(dim(doublets_multiome_df))


# Save
saveRDS(
  doublets_multiome_df,
  str_c(path_to_dfs, "/doublets_multiome_df_all.rds")
)
