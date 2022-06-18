#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=BZEECBEG_GXkc6Q1y
#SBATCH --workdir=.
#SBATCH --error=./log/BZEECBEG_GXkc6Q1y.err
#SBATCH --output=./log/BZEECBEG_GXkc6Q1y.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id BZEECBEG_GXkc6Q1y --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
