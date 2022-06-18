#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name=HjqdPU0E_aGDmEY5F
#SBATCH --workdir=.
#SBATCH --error=./log/HjqdPU0E_aGDmEY5F.err
#SBATCH --output=./log/HjqdPU0E_aGDmEY5F.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/6.0.1/cellranger multi --id HjqdPU0E_aGDmEY5F --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
