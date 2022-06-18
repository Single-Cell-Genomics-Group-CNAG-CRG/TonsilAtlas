#!/bin/bash
#SBATCH --job-name="co7dzuup_xuczw9vc"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/co7dzuup_xuczw9vc_%x_%J.err
#SBATCH --output=./log/co7dzuup_xuczw9vc_%x_%J.out
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger-arc/1.0.0/cellranger-arc-1.0.0/cellranger-arc count --id co7dzuup_xuczw9vc --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-arc-GRCh38-2020-A --libraries libraries.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished
