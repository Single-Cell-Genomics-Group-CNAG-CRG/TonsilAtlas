#!/bin/bash
#SBATCH --job-name="ryh4el3i_biv0w7ca"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/ryh4el3i_biv0w7ca_%x_%J.err
#SBATCH --output=./log/ryh4el3i_biv0w7ca_%x_%J.out
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

export HDF5_USE_FILE_LOCKING="FALSE"

/scratch/groups/hheyn/software/cellranger-arc/1.0.0/cellranger-arc-1.0.0/cellranger-arc count --id ryh4el3i_biv0w7ca --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-arc-GRCh38-2020-A --libraries libraries.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished
