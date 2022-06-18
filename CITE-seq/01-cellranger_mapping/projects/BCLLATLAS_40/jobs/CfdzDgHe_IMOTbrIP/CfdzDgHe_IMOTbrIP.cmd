#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=CfdzDgHe_IMOTbrIP
#SBATCH --workdir=.
#SBATCH --error=./log/CfdzDgHe_IMOTbrIP.err
#SBATCH --output=./log/CfdzDgHe_IMOTbrIP.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id CfdzDgHe_IMOTbrIP --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
