#!/bin/bash 


#SBATCH --job-name="wf4su8ny_h4yj8bv7"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/wf4su8ny_h4yj8bv7_%x_%J.err
#SBATCH --output=./log/wf4su8ny_h4yj8bv7_%x_%J.out
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id wf4su8ny_h4yj8bv7 --chemistry SC3Pv3 --expect-cells 20000 --localcores 24 --localmem 62 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;

echo [`date "+%Y-%m-%d %T"`] job finished