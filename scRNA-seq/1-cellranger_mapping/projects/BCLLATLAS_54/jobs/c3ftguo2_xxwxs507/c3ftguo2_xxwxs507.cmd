#!/bin/bash 


#SBATCH --job-name="c3ftguo2_xxwxs507"
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/c3ftguo2_xxwxs507_%x_%J.err
#SBATCH --output=./log/c3ftguo2_xxwxs507_%x_%J.out
#SBATCH -t 10:00:00
#SBATCH -n 1
#SBATCH -c 24
#SBATCH -q normal
#SBATCH -p genD



echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

export HDF5_USE_FILE_LOCKING=FALSE 

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id c3ftguo2_xxwxs507 --chemistry SC3Pv3 --expect-cells 20000 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A --jobmode /scratch/groups/singlecell/software/cellranger/6.1.1/external/martian/jobmanagers/new_cluster/slurm.template;

echo [`date "+%Y-%m-%d %T"`] job finished
