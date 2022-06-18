#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=WToIzInl_LudU7hVX
#SBATCH --workdir=.
#SBATCH --error=./log/WToIzInl_LudU7hVX.err
#SBATCH --output=./log/WToIzInl_LudU7hVX.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id WToIzInl_LudU7hVX --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
