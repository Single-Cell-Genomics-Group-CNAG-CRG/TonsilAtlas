#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=SOJZt9kY_qpnv20QN
#SBATCH --workdir=.
#SBATCH --error=./log/SOJZt9kY_qpnv20QN.err
#SBATCH --output=./log/SOJZt9kY_qpnv20QN.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id SOJZt9kY_qpnv20QN --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
