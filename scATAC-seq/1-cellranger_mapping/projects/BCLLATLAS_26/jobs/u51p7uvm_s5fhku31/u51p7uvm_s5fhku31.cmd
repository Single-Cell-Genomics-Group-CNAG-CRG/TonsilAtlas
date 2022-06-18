#!/bin/bash
    

#SBATCH --job-name=u51p7uvm_s5fhku31
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/u51p7uvm_s5fhku31_%x_%J.err
#SBATCH --output=./log/u51p7uvm_s5fhku31_%x_%J.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_26/jobs/u51p7uvm_s5fhku31/fastq --id u51p7uvm_s5fhku31 --sample u51p7uvm_s5fhku31 --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished