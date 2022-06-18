#!/bin/bash
    

#SBATCH --job-name=lacrul5e_67e3wsyf
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/lacrul5e_67e3wsyf_%x_%J.err
#SBATCH --output=./log/lacrul5e_67e3wsyf_%x_%J.out
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

/scratch/groups/hheyn/software/cellranger-atac/1.2.0/cellranger-atac count --fastqs /scratch/devel/rmassoni/tonsil_atlas/current/scATAC_seq/1-cellranger_mapping/projects/BCLLATLAS_30/jobs/lacrul5e_67e3wsyf/fastq --id lacrul5e_67e3wsyf --sample lacrul5e_67e3wsyf --localcores 24 --localmem 64 --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-atac-GRCh38-1.2.0;

echo [`date "+%Y-%m-%d %T"`] job finished
