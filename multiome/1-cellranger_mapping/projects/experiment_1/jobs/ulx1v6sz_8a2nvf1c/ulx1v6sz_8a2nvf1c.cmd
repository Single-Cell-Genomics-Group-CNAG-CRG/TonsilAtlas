#!/bin/bash
#SBATCH --job-name="ulx1v6sz_8a2nvf1c"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/ulx1v6sz_8a2nvf1c_%x_%J.err
#SBATCH --output=./log/ulx1v6sz_8a2nvf1c_%x_%J.out
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

/scratch/groups/hheyn/software/cellranger-arc/1.0.0/cellranger-arc-1.0.0/cellranger-arc count --id ulx1v6sz_8a2nvf1c --reference /scratch/groups/hheyn/data/reference/refdata-cellranger-arc-GRCh38-2020-A --libraries libraries.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished
