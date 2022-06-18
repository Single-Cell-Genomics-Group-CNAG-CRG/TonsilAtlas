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
# import config_vars as cfg
from utils import *

# Define command-line arguments
parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts of this project")
parser.add_argument("--spaceranger",
                    dest = "spaceranger",
                    action = "store",
                    default = None,
                    help = "Path to spaceranger software - /scratch/groups/hheyn/software/spaceranger/1.1.0/")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--reference_path",
                    dest = "reference_path",
                    action = "store",
                    default = None,
                    help = "Path of folder containing 10x-compatible reference - /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A")
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
                    help = "Metadata csv file for the tonsil atlas project, visium part.")
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
# parser.add_argument("--slide",
#                     dest = "slide",
#                     action = "store",
#                     default = None,
#                     help = "Visium slide serial number, for example 'V10J25-015'")
# parser.add_argument("--area",
#                     dest = "area",
#                     action = "store",
#                     default = None,
#                     help = "Visium area identifier, for example 'A1'")
# parser.add_argument("--slidefile",
#                     dest = "slidefile",
#                     action = "store",
#                     default = None,
#                     help = "Spot design file for your slide, downloaded from 10x Genomics. NOTE: this is only required if your machine doesn't have internet access. You must still pass --slide and --area when using this argument")

options = parser.parse_args()
spaceranger = options.spaceranger
reference_path = options.reference_path
subproject = options.subproject
gem_id = options.gem_id
print(gem_id)
metadata_path = options.metadata
# feat_ref = options.feat_ref
# fastq_paths = options.fastq_paths
# slide = options.slide
# area = options.area
# slidefile = options.slidefile

# spaceranger = "/scratch/devel/melosua/phd/10x_software/spaceranger-1.1.0/spaceranger"
# reference_path = "/scratch/devel/melosua/phd/10x_software/refdata-gex-GRCh38-2020-A/"
# subproject = "BCLLATLAS_32"
# gem_id = "c28w2r_7jne4i"
# metadata_path = "../data/sample_id.txt"
# slide = "V19S23-039"
# area = "A1"
# slidefile = "data/spot_layout/V19S23-039.gpr"

# Read data
project_dir = "../projects/{}".format(subproject)
# fastq_path_df = pd.read_csv(fastq_paths, sep = ",", header = 0)
metadata_df = pd.read_csv(metadata_path, sep = ",", header = 0)
# feature_reference_df = pd.read_csv(feat_ref, sep = ",", header = 0)

if options.verbose:
    sys.stderr.write("Files read successfully!\n")

# For each sample, create directories and jobscript
if not os.path.exists("{}/jobs".format(project_dir)):
    os.mkdir("{}/jobs".format(project_dir))


filt = (metadata_df["gem_id"] == gem_id)
metadata_df = metadata_df.loc[filt]

# Create directories if they don't exist
subproject_dir = "{}/jobs/{}".format(project_dir, gem_id)
fastq_dir = "{}/fastq".format(subproject_dir)
output_dir = "{}/output".format(subproject_dir)

for direct in [subproject_dir, fastq_dir, output_dir]:
    if not os.path.exists(direct):    
        os.mkdir(direct)

# Define variables and subset dataframes
# library_id = metadata_df.loc[filt, "library_id"]
# fastq_sub_df = fastq_path_df.loc[fastq_path_df["library_id"].isin(library_id), :]
# type = metadata_df["type"]
# type = type.values[0]

# Create symmlinks to fastq files
# create_fastq_symlink_nh(gem_id, fastq_sub_df, fastq_dir)

# Create spaceranger script
if len(metadata_df) == 1:
    make_spaceranger_clustermode(subproject_dir, metadata_df, spaceranger, reference_path)
else:
    print("metadata_fd has number of rows != 1; {}".format(len(metadata_df)))
