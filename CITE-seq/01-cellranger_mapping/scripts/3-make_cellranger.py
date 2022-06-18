# This script initializes the filesystem of this project:
# It creates a "jobs" folder which contains as many subdirectories as samples it has
# For each sample directory, it creates the following files/folders:
# 1. fastq: dir with the symlinks pointing to the fastq files
# 2. output: dir which contains standard error and output of cellranger 
# 3. (sample_id).cmd: job script to compute the features-barcode matrix using cellranger


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
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--subproject",
					dest = "subproject",
					action = "store",
					default = None,
					help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--gem_id",
					dest = "gem_id",
					action = "store",
					default = None,
					help = "Gel Beads in Emulsion id")
parser.add_argument("--verbose",
					dest = "verbose",
					action = "store_true",
					default = False,
					help = "Print log in standard error")
parser.add_argument("--metadata",
					dest = "metadata",
					action = "store",
					default = None,
					help = "Metadata csv file for the tonsil atlas project")
parser.add_argument("--feat_ref",
					dest = "feat_ref",
					action = "store",
					default = None,
					help = "Feature hashtag reference correspondence path")
parser.add_argument("--fastq_paths",
					dest = "fastq_paths",
					action = "store",
					default = None,
					help = "File that contains the paths of the fastqs for the subproject libraries")
options = parser.parse_args()
subproject = options.subproject
gem_id = options.gem_id
metadata_path = options.metadata
feat_ref = options.feat_ref
fastq_paths = options.fastq_paths


# Read data
project_dir = "projects/{}".format(subproject)
fastq_path_df = pd.read_csv(fastq_paths, sep = ",", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
feature_reference_df = pd.read_csv(feat_ref, sep = ",", header = 0)
if options.verbose:
	sys.stderr.write("Files read successfully!\n")


# Create directories
if not os.path.exists("{}/jobs".format(project_dir)):
	os.mkdir("{}/jobs".format(project_dir))
filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]
subproject_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(subproject_dir)
output_dir = "{}/log".format(subproject_dir)
for direct in [subproject_dir, fastq_dir, output_dir]:
	if not os.path.exists(direct):	
		os.mkdir(direct)

if not os.path.exists("{}/gex".format(fastq_dir)):
	os.mkdir("{}/gex".format(fastq_dir))
if not os.path.exists("{}/csp".format(fastq_dir)):
	os.mkdir("{}/csp".format(fastq_dir))
if not os.path.exists("{}/bcr".format(fastq_dir)):
	os.mkdir("{}/bcr".format(fastq_dir))
if not os.path.exists("{}/tcr".format(fastq_dir)):
	os.mkdir("{}/tcr".format(fastq_dir))


# Create symmlinks to fastq files
library_ids = metadata_df.loc[filt, "library_id"]
fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_ids), :]
type = metadata_df["type"]
type = type.values[0]
for lib in library_ids:
	lib_type = metadata_df.loc[metadata_df["library_id"] == lib, "type"]
	lib_type = lib_type.values[0]
	if lib_type == "GEX":
		symlink_path = "{}/gex".format(fastq_dir)
	elif lib_type == "CSP":
		symlink_path = "{}/csp".format(fastq_dir)
	elif lib_type == "BCR":
		symlink_path = "{}/bcr".format(fastq_dir)
	elif lib_type == "TCR":
		symlink_path = "{}/tcr".format(fastq_dir)	
	create_fastq_symlink_cite_seq(gem_id, lib, lib_type, fastq_sub_df, symlink_path)
	 

# Create libraries.csv: file indicating fastq path and type of reads (protein/cDNA)
write_multi_csv(gem_id, os.path.abspath(subproject_dir))


# Create cellranger script
make_cellranger_multi(gem_id, subproject_dir, 10000)

