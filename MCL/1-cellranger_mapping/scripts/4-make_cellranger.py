# Import required packages
import numpy as np
import pandas as pd
import os
import argparse
import subprocess
import re
import sys
import config_vars as cfg
from utils import *


# Define command-line arguments
parser = argparse.ArgumentParser(description = "Create script to run cellranger")
parser.add_argument("--experiment",
                    dest = "experiment",
                    action = "store",
                    default = None,
                    help = "multiome experiment")
parser.add_argument("--gem_id",
		    dest = "gem_id",
		    action = "store",
		    default = None,
		    help = "Gel Beads in Emulsion id")
options = parser.parse_args()
experiment = options.experiment
gem_id = options.gem_id


# Read data
project_dir = "projects/{}".format(experiment)
fastq_path_rna = "{}/fastq_paths_rna.csv".format(project_dir)
fastq_path_atac = "{}/fastq_paths_atac.csv".format(project_dir)
fastq_path_rna_df = pd.read_csv(fastq_path_rna, sep = ",", header = 0)
fastq_path_atac_df = pd.read_csv(fastq_path_atac, sep = ",", header = 0)
fastq_path_df = fastq_path_rna_df.append(fastq_path_atac_df)
metadata_path = "data/sequencing_metadata.csv" 
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)



# Create directories
if not os.path.exists("{}/jobs".format(project_dir)):
	os.mkdir("{}/jobs".format(project_dir))
filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]
gem_id_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(gem_id_dir)
log_dir = "{}/log".format(gem_id_dir)
for direct in [gem_id_dir, fastq_dir, log_dir]:
	if not os.path.exists(direct):	
		os.mkdir(direct)
if not os.path.exists("{}/RNA".format(fastq_dir)):
	os.mkdir("{}/RNA".format(fastq_dir))
if not os.path.exists("{}/ATAC".format(fastq_dir)):
	os.mkdir("{}/ATAC".format(fastq_dir))


# Create symmlinks to fastq files
library_ids = metadata_df.loc[filt, "library_id"]
fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_ids), :]
for lib in library_ids:
	lib_type = metadata_df.loc[metadata_df["library_id"] == lib, "type"]
	lib_type = lib_type.values[0]
	if lib_type == "RNA":
		symlink_path = "{}/RNA".format(fastq_dir)
	elif lib_type == "ATAC":
		symlink_path = "{}/ATAC".format(fastq_dir)
	create_fastq_symlink(gem_id, lib, lib_type, fastq_sub_df, symlink_path)


# Create libraries.csv: file indicating fastq path and type of reads (RNA/ATAC)
write_libraries_csv(gem_id, os.path.abspath(gem_id_dir))


# Create cellranger script
make_cellranger_arc(gem_id, gem_id_dir)


