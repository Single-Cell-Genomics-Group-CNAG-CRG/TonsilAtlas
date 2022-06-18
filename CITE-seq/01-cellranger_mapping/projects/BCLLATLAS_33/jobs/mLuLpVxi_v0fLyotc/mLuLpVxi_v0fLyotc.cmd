#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=mLuLpVxi_v0fLyotc
#SBATCH --workdir=.
#SBATCH --error=./log/mLuLpVxi_v0fLyotc.err
#SBATCH --output=./log/mLuLpVxi_v0fLyotc.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id mLuLpVxi_v0fLyotc --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
