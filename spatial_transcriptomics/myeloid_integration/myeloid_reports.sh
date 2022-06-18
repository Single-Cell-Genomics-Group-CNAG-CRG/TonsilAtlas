#!/usr/bin/env bash

# Script to be run from where the Rmarkdown document is
set -euo pipefail

# Myeloid report scripts
R -e "rmarkdown::render('myeloid_signatures.Rmd', output_file='myeloid_signatures.html')"
                        
R -e "rmarkdown::render('myeloid_signature_projection.Rmd',  output_file='myeloid_signature_projection.html')"

R -e "rmarkdown::render('myeloid-deconvolution.Rmd',  output_file='myeloid-deconvolution.html')"

R -e "rmarkdown::render('myeloid-deconvolution-visualization.Rmd',  output_file='myeloid-deconvolution-visualization.html')"

R -e "rmarkdown::render('MAGIC_denoising-myeloid.Rmd',  output_file='MAGIC_denoising-myeloid.html')"

R -e "rmarkdown::render('MAGIC_visualization-myeloid.Rmd',  output_file='MAGIC_visualization-myeloid.html')"
