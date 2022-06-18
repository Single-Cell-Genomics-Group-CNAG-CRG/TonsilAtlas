# Writes fastq path by arranging proper flowcell, lane, index and read for a set of libraries

# Load packages
import numpy as np
import pandas as pd
import os
import argparse


# Define command-line arguments
parser = argparse.ArgumentParser(description = "Write fastq paths by using a unique combination of flow cell, lane and index")
parser.add_argument("--subproject",
                    dest = "subproject",
                    action = "store",
                    default = None,
                    help = "Subproject we are working on (i.e. BCLLATLAS_10)")
parser.add_argument("--info_file",
                    dest = "info_file",
                    action = "store",
                    default = None,
                    help = "Tab-delimited file with the information of Illumina sequence of libraries for that subproject")
parser.add_argument("--metadata",
                    dest = "metadata",
                    action = "store",
                    default = None,
                    help = "Metadata csv file for the tonsil atlas project")
options = parser.parse_args()
subproject = options.subproject
info_file = options.info_file

# Read and subset files
metadata_path = options.metadata
tonsil_metadata = pd.read_csv(metadata_path)
mask = (tonsil_metadata["subproject"] == subproject)
libraries = tonsil_metadata.loc[mask, "library_id"]
libraries = list(libraries)
lims = pd.read_csv(info_file, sep = "\t", header = 0)
lims = lims.loc[lims.id.isin(libraries)]


# Assemble fastq paths combining flowcell, lane and index
fastq_path = "/scratch/project/production/fastq"
fastq_dict = {"library_id":[], "fastq_path":[], "read":[], "pair_id":[]}
for idx in lims.index:
    fc = lims.loc[idx, "flowcell"]
    lane = lims.loc[idx, "lane"]
    index = lims.loc[idx, "index"]
    fastq_path_r1 = "{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_r2 = "{}/{}/{}/fastq/{}_{}_{}_UMI2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_r3 = "{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, index)
    fastq_path_list = [fastq_path_r1, fastq_path_r2, fastq_path_r3]
    library_id_list = [lims.loc[idx, "id"]] * 3
    read_list = ["R1", "R2", "R3"]
    pair_id_list = ["P" + str(idx + 1)] * 3
    fastq_dict["library_id"].extend(library_id_list)
    fastq_dict["fastq_path"].extend(fastq_path_list)
    fastq_dict["read"].extend(read_list)
    fastq_dict["pair_id"].extend(pair_id_list)

fastq_df = pd.DataFrame(fastq_dict)


# Write dataframe to csv
fastq_df.to_csv("projects/{}/fastq_paths.csv".format(subproject), header = True, index = False)
