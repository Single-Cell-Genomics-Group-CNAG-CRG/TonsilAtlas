# Writes fastq path by arranging proper flowcell, lane, index and read for a set of libraries

# Load packages
import numpy as np
import pandas as pd
import os
import argparse


# Define command-line arguments
parser = argparse.ArgumentParser(description = "Write paths to fastq files by combining flow cell, lane and index")
parser.add_argument("--experiment",
                    dest = "experiment",
                    action = "store",
                    default = None,
                    help = "Experiment performed")
parser.add_argument("--metadata",
                    dest = "metadata",
                    action = "store",
                    default = None,
                    help = "Metadata csv file for that project")
options = parser.parse_args()
experiment = options.experiment
metadata_path = options.metadata


# Read and subset files
metadata = pd.read_csv(metadata_path)
if experiment == "experiment_1":
    subproject = "BCLLATLAS_64"
#elif experiment == "experiment_2":
#    subproject = "BCLLATLAS_47"
mask = (metadata["subproject"] == subproject)
libraries = metadata.loc[mask, "library_id"]
libraries = list(libraries)
path_to_info = "projects/{}/info_{}.txt".format(experiment, subproject)
lims = pd.read_csv(path_to_info, sep = "\t", header = 0)


# Assemble fastq paths combining flowcell, lane and index
fastq_path = "/scratch/project/production/fastq"
fastq_dict = {"library_id":[], "fastq_path":[], "read":[], "pair_id":[]}
for idx in lims.index:
    fc = lims.loc[idx, "flowcell"]
    lane = lims.loc[idx, "lane"]
    index = lims.loc[idx, "index"]
    fastq_path_r1 = "{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_r2 = "{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_list = [fastq_path_r1, fastq_path_r2]
    library_id_list = [lims.loc[idx, "id"]] * 2
    read_list = ["R1", "R2"]
    pair_id_list = ["P" + str(idx + 1)] * 2
    fastq_dict["library_id"].extend(library_id_list)
    fastq_dict["fastq_path"].extend(fastq_path_list)
    fastq_dict["read"].extend(read_list)
    fastq_dict["pair_id"].extend(pair_id_list)

fastq_df = pd.DataFrame(fastq_dict)


# Write dataframe to csv
fastq_df.to_csv("projects/{}/fastq_paths_rna.csv".format(experiment), header = True, index = False)
