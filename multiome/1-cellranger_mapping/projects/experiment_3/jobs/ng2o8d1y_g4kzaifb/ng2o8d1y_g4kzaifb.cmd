#!/bin/bash
#SBATCH --job-name="ng2o8d1y_g4kzaifb"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/ng2o8d1y_g4kzaifb_%x_%J.err
#SBATCH --output=./log/ng2o8d1y_g4kzaifb_%x_%J.out
#SBATCH -t 23:59:00
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -q long
#SBATCH -p genD

echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger-arc/1.0.0/cellranger-arc-1.0.0/cellranger-arc count --id ng2o8d1y_g4kzaifb --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-arc-GRCh38-2020-A --libraries libraries.csv;

echo [`date "+%Y-%m-%d %T"`] job finished
