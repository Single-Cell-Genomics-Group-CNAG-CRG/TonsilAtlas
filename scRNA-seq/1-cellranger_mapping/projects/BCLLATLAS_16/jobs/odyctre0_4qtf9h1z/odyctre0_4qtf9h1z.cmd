#!/bin/bash
    
#SBATCH --job-name="odyctre0_4qtf9h1z"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/odyctre0_4qtf9h1z_%x_%J.err
#SBATCH --output=./log/odyctre0_4qtf9h1z_%x_%J.out
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_16/jobs/odyctre0_4qtf9h1z/fastq --id odyctre0_4qtf9h1z --chemistry SC3Pv3 --expect-cells 5000 --localcores 24 --localmem 64 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;


echo [`date "+%Y-%m-%d %T"`] job finished