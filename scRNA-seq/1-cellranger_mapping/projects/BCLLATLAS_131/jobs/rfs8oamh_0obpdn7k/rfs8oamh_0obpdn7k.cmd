#!/bin/bash 


#SBATCH --job-name="rfs8oamh_0obpdn7k"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/%x_%J.err
#SBATCH --output=./log/%x_%J.out
#SBATCH -t 23:59:00
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -q long
#SBATCH -p genD
#SBATCH --mem=120G

echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id rfs8oamh_0obpdn7k --chemistry SC3Pv3 --expect-cells 20000 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;

echo [`date "+%Y-%m-%d %T"`] job finished
