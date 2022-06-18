#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=B20O1bh7_VmM99YZJ
#SBATCH --workdir=.
#SBATCH --error=./log/B20O1bh7_VmM99YZJ.err
#SBATCH --output=./log/B20O1bh7_VmM99YZJ.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id B20O1bh7_VmM99YZJ --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
