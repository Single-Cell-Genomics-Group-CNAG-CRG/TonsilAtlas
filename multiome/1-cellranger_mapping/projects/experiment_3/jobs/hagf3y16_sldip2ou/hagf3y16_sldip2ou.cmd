#!/bin/bash
#SBATCH --job-name="hagf3y16_sldip2ou"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/hagf3y16_sldip2ou_%J.err
#SBATCH --output=./log/hagf3y16_sldip2ou_%J.out
#SBATCH -t 23:59:00
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -q long
#SBATCH -p genD


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger-arc/1.0.0/cellranger-arc-1.0.0/cellranger-arc count --id hagf3y16_sldip2ou --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-arc-GRCh38-2020-A --libraries libraries.csv;

echo [`date "+%Y-%m-%d %T"`] job finished
