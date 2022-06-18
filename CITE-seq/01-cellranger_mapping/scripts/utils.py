# Python functions used in different scripts

import numpy as np
import pandas as pd
import subprocess
import os
import config_vars as cfg
import re


def create_fastq_symlink_cite_seq(gem_id, library, lib_type, fastq_path_df, symlink_path):
	"""Creates a symbolic link pointing to a fastq file using cellranger notation for cell-hashed samples.
	
	Args:
	  gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
	  library: library id
	  lib_type: type of library (cDNA or protein)
	  fastq_path_df: pandas dataframe with the fastq paths for that gem_id
	  symlink_path: string specifying where to create the symlinks (fastq/protein or fastq/cDNA)
	
	Returns:
	  None 
	"""
	fastq_path_sub = fastq_path_df.loc[fastq_path_df["library_id"] == library, :]
	pair_ids = np.unique(fastq_path_sub["pair_id"])
	for i in range(len(pair_ids)):
		filt = (fastq_path_df["pair_id"] == pair_ids[i])
		pair_df = fastq_path_df.loc[filt, :]
		for j in pair_df.index:
			lane = str(i + 1)
			symlink_path_lane = "{}/lane{}".format(symlink_path, lane)
			if not os.path.exists(symlink_path_lane):
				os.mkdir(symlink_path_lane)
			fastq_path = pair_df.loc[j, "fastq_path"]
			read = pair_df.loc[j, "read"]
			read = read.replace("R", "")
			if lib_type == "CSP":
				gem_id_sp = "{}_csp".format(gem_id)
			if lib_type == "GEX":
				gem_id_sp = "{}_gex".format(gem_id)
			if lib_type == "TCR":
				gem_id_sp = "{}_tcr".format(gem_id)
			if lib_type == "BCR":
				gem_id_sp = "{}_bcr".format(gem_id)
			subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path_lane, gem_id_sp, lane, read)])



def write_multi_csv(gem_id, gem_id_path):
	"""Creates the file "libraries.csv" which is required by cellranger in feature-barcoding analysis.
	
	Args:
	 gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
	 fastq_path: path to GEM-specific directory
	 gem_id_path: absolute path to gem_id-specific directory.
	
	Returns:
	  None 
	"""
	multi_csv = open("{}/multi.csv".format(gem_id_path), "w")
	multi_csv_header="""[gene-expression],,,,
reference,/scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A,,,
chemistry,SC5P-R2,,,
no-bam,TRUE,,,
,,,,
[feature],,,,
reference,/home/devel/srashmi/tonsil_atlas/1-cellranger_mapping/data/tonsil_atlas_cite_seq_reference.csv,,,
,,,,
[vdj],,,,
reference,/scratch/groups/hheyn/data/reference/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0,,,
,,,,
[libraries],,,,
fastq_id,fastqs,lanes,feature_types,subsample_rate
""".format(cfg.gex_reference,cfg.csp_reference,cfg.vdj_reference)

	multi_csv.write(multi_csv_header)
	fastq_dirs = os.listdir("{}/fastq".format(gem_id_path))
	for d in fastq_dirs:
		fastq_sub_dirs = os.listdir("{}/fastq/{}".format(gem_id_path, d))
		for sub_d in fastq_sub_dirs:
				if d == "csp":
					fastq_id = "{}_csp".format(gem_id)
					feature_type = "Antibody Capture"
				elif d == "gex":
					fastq_id = "{}_gex".format(gem_id)
					feature_type = "Gene Expression"
				elif d == "tcr":
					fastq_id = "{}_tcr".format(gem_id)
					feature_type = "vdj-t"
				elif d == "bcr":
					fastq_id = "{}_bcr".format(gem_id)
					feature_type = "vdj-b"
				fastqs = "{}/fastq/{}/{}".format(gem_id_path, d, sub_d)
				lane=re.findall('[0-9]+', sub_d)[0]	
				multi_csv.write(fastq_id+","+fastqs+","+lane+","+feature_type+"\n")
	multi_csv.close()
	
	
	
def make_cellranger_multi(gem_id, jobscript_path, expected_cells):
	"""Creates a cellranger script for a multi experiment (GEX, CITE-seq/CSP, VDJ) 
	
	Args:
	  gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
	  jobscript_path: path to save the jobscript
	  expected_cells: expected number of high-quality cells in this experiment

	Returns:
	  None 
	"""
	job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
	job_script = """#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH--job-name={}
#SBATCH --workdir=.
#SBATCH --error=./log/{}.err
#SBATCH --output=./log/{}.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

{} multi --id {} --csv multi.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished!!
""".format(gem_id,gem_id,gem_id, cfg.cellranger_5_path, gem_id)
	job_script_file.write(job_script)
	job_script_file.close()
	
