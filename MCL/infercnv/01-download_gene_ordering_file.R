# This script downloads the gene ordering file to run inferCNV

library(tidyverse)
library(here)
library(glue)
path_to_dir <- here::here()
path_to_gene_order_f <- glue::glue("{path_to_dir}/MCL/5-infercnv/tmp/gencode_v21_gen_pos.complete.txt")
gene_order_file <- read_tsv(
  "https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v21_gen_pos.complete.txt",
  col_names = c("gene", "chromosome", "start", "end")
)
write_tsv(gene_order_file, path_to_gene_order_f, col_names = FALSE)
