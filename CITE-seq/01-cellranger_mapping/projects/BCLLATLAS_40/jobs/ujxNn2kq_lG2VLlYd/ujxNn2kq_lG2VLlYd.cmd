#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=ujxNn2kq_lG2VLlYd
#SBATCH --workdir=.
#SBATCH --error=./log/ujxNn2kq_lG2VLlYd.err
#SBATCH --output=./log/ujxNn2kq_lG2VLlYd.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id ujxNn2kq_lG2VLlYd --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
