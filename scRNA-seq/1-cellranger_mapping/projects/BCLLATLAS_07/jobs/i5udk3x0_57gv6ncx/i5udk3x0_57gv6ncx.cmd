#!/bin/bash 


#SBATCH --job-name="i5udk3x0_57gv6ncx"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/i5udk3x0_57gv6ncx_%x_%J.err
#SBATCH --output=./log/i5udk3x0_57gv6ncx_%x_%J.out
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger/4.0.0/cellranger count --libraries libraries.csv --feature-ref feature_reference.csv --id i5udk3x0_57gv6ncx --chemistry SC3Pv3 --expect-cells 20000 --localcores 24 --localmem 62 --transcriptome /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A;

echo [`date "+%Y-%m-%d %T"`] job finished