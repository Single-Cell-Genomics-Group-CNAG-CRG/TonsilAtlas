#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=uqJAc4r9_BYScOzxA
#SBATCH --workdir=.
#SBATCH --error=./log/uqJAc4r9_BYScOzxA.err
#SBATCH --output=./log/uqJAc4r9_BYScOzxA.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id uqJAc4r9_BYScOzxA --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
