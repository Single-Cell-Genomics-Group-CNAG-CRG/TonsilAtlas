#!/usr/bin/env bash

# Script to be run from where the Rmarkdown document is
set -euo pipefail

# Myeloid report scripts
R -e "rmarkdown::render('MAGIC_denoising.Rmd',
                        output_file='MAGIC_denoising.html')"
                        
R -e "rmarkdown::render('MAGIC-trajectories-markers.Rmd',
                        output_file='MAGIC-trajectories-markers.html')"

R -e "rmarkdown::render('CD4-deconvolution.Rmd',
                        output_file='CD4-deconvolution.html')"

R -e "rmarkdown::render('CD4-deconvolution_assessment.Rmd',
                        output_file='CD4-deconvolution_assessment.html')"



