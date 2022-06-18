#!/bin/bash

#SBATCH --job-name="cd8_decon"
#SBATCH --workdir=.
#SBATCH --partition=genB,main
#SBATCH --mail-type=all
#SBATCH --mail-user=marc.elosua@cnag.crg.eu
#SBATCH --error=./cd8_deconv_%x_%J.err
#SBATCH --output=./cd8_deconv_%x_%J.out
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME
module purge
conda deactivate
conda activate spatial_r
R -e "rmarkdown::render('CD8-deconvolution.Rmd', output_file='CD8-deconvolution.html')"
echo [`date "+%Y-%m-%d %T"`] job finished
