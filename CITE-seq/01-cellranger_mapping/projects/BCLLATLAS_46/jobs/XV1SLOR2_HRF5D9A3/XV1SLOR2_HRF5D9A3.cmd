#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=XV1SLOR2_HRF5D9A3
#SBATCH --workdir=.
#SBATCH --error=./log/XV1SLOR2_HRF5D9A3.err
#SBATCH --output=./log/XV1SLOR2_HRF5D9A3.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id XV1SLOR2_HRF5D9A3 --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
