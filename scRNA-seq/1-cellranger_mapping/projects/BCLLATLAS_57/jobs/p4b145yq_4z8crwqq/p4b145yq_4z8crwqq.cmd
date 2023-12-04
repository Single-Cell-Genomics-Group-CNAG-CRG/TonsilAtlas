#!/bin/bash 


#SBATCH --job-name="p4b145yq_4z8crwqq"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/p4b145yq_4z8crwqq_%J.err
#SBATCH --output=./log/p4b145yq_4z8crwqq_%J.out
#SBATCH -t 11:59:00
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -q normal
#SBATCH -p genD


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

export HDF5_USE_FILE_LOCKING=FALSE 

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id p4b145yq_4z8crwqq --chemistry SC3Pv3 --expect-cells 20000 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;

echo [`date "+%Y-%m-%d %T"`] job finished
