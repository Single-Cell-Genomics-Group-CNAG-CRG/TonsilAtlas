#!/bin/bash
    
#SBATCH --job-name="wf4su8ny_h4yj8bv7"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/velocyto_wf4su8ny_h4yj8bv7_%J.err
#SBATCH --output=./log/velocyto_wf4su8ny_h4yj8bv7_%J.out
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

export HDF5_USE_FILE_LOCKING="FALSE"

REPEAT_MASK_GTF="/scratch/groups/hheyn/data/reference/GRCh38_repeat_mask/GRCh38_repeat_mask.gtf"
CELLRANGER_OUTPUT_DIR="/scratch/devel/rmassoni/tonsil_atlas/current/scRNA-seq/1-cellranger_mapping/projects/BCLLATLAS_25/jobs/wf4su8ny_h4yj8bv7/wf4su8ny_h4yj8bv7"
REFERENCE_GTF="/scratch/devel/rmassoni/reference/human/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"

velocyto run10x -m $REPEAT_MASK_GTF $CELLRANGER_OUTPUT_DIR $REFERENCE_GTF

echo [`date "+%Y-%m-%d %T"`] job finished