#!/usr/bin/env bash

# Script to be run from where the Rmarkdown document is
set -euo pipefail

# Myeloid report scripts
R -e "rmarkdown::render('epithelium_signatures.Rmd', output_file='epithelium_signatures.html')"
                        
R -e "rmarkdown::render('epithelium_signatures_projection.Rmd', output_file='epithelium_signatures_projection.html')"

R -e "rmarkdown::render('epithelium-deconvolution.Rmd', output_file='epithelium-deconvolution.html')"

R -e "rmarkdown::render('epithelium-deconvolution-viz.Rmd', output_file='epithelium-deconvolution-viz.html')"
