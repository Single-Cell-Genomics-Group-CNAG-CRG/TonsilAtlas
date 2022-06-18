#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=rdFRFhrU_ZdYeOZlf
#SBATCH --workdir=.
#SBATCH --error=./log/rdFRFhrU_ZdYeOZlf.err
#SBATCH --output=./log/rdFRFhrU_ZdYeOZlf.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id rdFRFhrU_ZdYeOZlf --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
