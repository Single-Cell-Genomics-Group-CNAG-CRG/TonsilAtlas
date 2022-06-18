#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=ifZOgenn_TpMNTvBa
#SBATCH --workdir=.
#SBATCH --error=./log/ifZOgenn_TpMNTvBa.err
#SBATCH --output=./log/ifZOgenn_TpMNTvBa.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id ifZOgenn_TpMNTvBa --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
