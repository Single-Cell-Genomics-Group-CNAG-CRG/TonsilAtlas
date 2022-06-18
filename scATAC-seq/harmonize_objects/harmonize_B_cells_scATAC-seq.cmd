#!/bin/bash 
#SBATCH --job-name="harmonize_B"
#SBATCH --error=./log/harmonize_B.err
#SBATCH --output=./log/harmonize_B.log
#SBATCH --cpus-per-task=12
#SBATCH --time=03:00:00
#SBATCH --mail-user=ramonmassoni@gmail.com

module load gcc/3.6.0 gsl/2.4 gmp/6.1.2 hdf5/1.10.1 R/4.1.1 PANDOC
export R_LIBS_USER=~/R_libs_repo/4.1.1/R_libs/

echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

R -e "library(BiocStyle); path_to_knit <- '/scratch/devel/rmassoni/tonsil_atlas/current'; rmarkdown::render('harmonize_B_cells_scATAC-seq.Rmd', knit_root_dir = path_to_knit)"

echo [`date "+%Y-%m-%d %T"`] job finished
