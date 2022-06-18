# This script creates a GEM-specific job to run velocyto


# Import required packages
import numpy as np
import pandas as pd
import argparse
import os
import config_vars as cfg


# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--gem_id",
                    dest = "gem_id",
                    action = "store",
                    default = None,
                    help = "Gel Beads in Emulsion id")
options = parser.parse_args()
gem_id = options.gem_id


# Read metadata
metadata_df = pd.read_csv("data/tonsil_atlas_metadata.csv", sep = ",", header = 0)
mask = (metadata_df["gem_id"] == gem_id)
subproject = metadata_df.loc[mask, "subproject"]
subproject = subproject.values[0]


# Write velocyto jobscript
jobscript_path = "{}projects/{}/jobs/{}/velocyto_{}.cmd".format(cfg.project_path, subproject, gem_id, gem_id)
cellranger_output_dir = "{}projects/{}/jobs/{}/{}".format(cfg.project_path, subproject, gem_id, gem_id)
job_script_file = open(jobscript_path, "w")
job_script = """#!/bin/bash
    
#SBATCH --job-name="{}"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/velocyto_{}_%J.err
#SBATCH --output=./log/velocyto_{}_%J.out
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=20


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

export HDF5_USE_FILE_LOCKING="FALSE"

REPEAT_MASK_GTF="{}"
CELLRANGER_OUTPUT_DIR="{}"
REFERENCE_GTF="{}"

velocyto run10x -m $REPEAT_MASK_GTF $CELLRANGER_OUTPUT_DIR $REFERENCE_GTF

echo [`date "+%Y-%m-%d %T"`] job finished""".format(gem_id, gem_id, gem_id, cfg.REPEAT_MASK_GTF, cellranger_output_dir, cfg.REFERENCE_GTF)
job_script_file.write(job_script)
job_script_file.close()

