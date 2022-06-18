#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=pseMjZsU_qgZNOhOQ
#SBATCH --workdir=.
#SBATCH --error=./log/pseMjZsU_qgZNOhOQ.err
#SBATCH --output=./log/pseMjZsU_qgZNOhOQ.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id pseMjZsU_qgZNOhOQ --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
