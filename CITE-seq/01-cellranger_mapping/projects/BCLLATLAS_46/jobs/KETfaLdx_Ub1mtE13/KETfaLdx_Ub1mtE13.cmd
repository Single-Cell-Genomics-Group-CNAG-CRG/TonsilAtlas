#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=KETfaLdx_Ub1mtE13
#SBATCH --workdir=.
#SBATCH --error=./log/KETfaLdx_Ub1mtE13.err
#SBATCH --output=./log/KETfaLdx_Ub1mtE13.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id KETfaLdx_Ub1mtE13 --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
