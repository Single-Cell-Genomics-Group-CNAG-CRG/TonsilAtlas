#!/usr/bin/env bash

# Script to be run from where the Rmarkdown document is
set -euo pipefail

# Myeloid report scripts

R -e "rmarkdown::render('MAGIC_denoising_fibroblasts.Rmd',  output_file='MAGIC_denoising_fibroblasts.html')"

R -e "rmarkdown::render('MAGIC_visualization_fibroblasts.Rmd',  output_file='MAGIC_visualization_fibroblasts.html')"
