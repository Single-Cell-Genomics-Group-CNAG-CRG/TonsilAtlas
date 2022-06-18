# Python functions used in different scripts

import numpy as np
import pandas as pd
import subprocess
import os
import config_vars as cfg


def create_fastq_symlink(gem_id, lib, lib_type, fastq_path_df, symlink_path):
	"""Creates a symbolic link pointing to a fastq file using cellranger notation.
	
	Args:
		gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
		library: library id
		lib_type: type of library (RNA or ATAC)
		fastq_path_df: pandas dataframe with the fastq paths for that gem_id
		symlink_path: string specifying where to create the symlinks (fastq/RNA or fastq/ATAC)
	
	Returns:
		None 
	"""
	fastq_path_sub = fastq_path_df.loc[fastq_path_df["library_id"] == lib, :]
	pair_ids = np.unique(fastq_path_sub["pair_id"])
	for i in range(len(pair_ids)):
		filt = (fastq_path_sub["pair_id"] == pair_ids[i])
		pair_df = fastq_path_sub.loc[filt, :]
		for j in pair_df.index:
			lane = str(i + 1)
			symlink_path_lane = "{}/lane{}".format(symlink_path, lane)
			if not os.path.exists(symlink_path_lane):
				os.mkdir(symlink_path_lane)
			fastq_path = pair_df.loc[j, "fastq_path"]
			read = pair_df.loc[j, "read"]
			read = read.replace("R", "")
			if lib_type == "RNA":
				gem_id_sp = "{}_rna".format(gem_id)
			if lib_type == "ATAC":
				gem_id_sp = "{}_atac".format(gem_id)
			subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path_lane, gem_id_sp, lane, read)])



def write_libraries_csv(gem_id, gem_id_path):
	"""Creates the file "libraries.csv" which is required by cellranger-arc.
	
	Args:
	 gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
	 gem_id_path: absolute path to gem_id-specific directory.
	
	Returns:
	  None 
	"""
	lib_csv = open("{}/libraries.csv".format(gem_id_path), "w")
	lib_csv.write("fastqs,sample,library_type")
	fastq_dirs = os.listdir("{}/fastq".format(gem_id_path))
	for d in fastq_dirs:
		if d == "RNA":
			gem_id_sp = "{}_rna".format(gem_id)
			lib_type = "Gene Expression"
		elif d == "ATAC":
			gem_id_sp = "{}_atac".format(gem_id)
			lib_type = "Chromatin Accessibility"
		fastq_sub_dirs = os.listdir("{}/fastq/{}".format(gem_id_path, d))
		for sub_d in fastq_sub_dirs:
			sub_d_abs_path = "{}/fastq/{}/{}".format(gem_id_path, d, sub_d)
			output_line = "\n{},{},{}".format(sub_d_abs_path, gem_id_sp, lib_type)    
			lib_csv.write(output_line)
	lib_csv.close()



def make_cellranger_arc(gem_id, jobscript_path):
	"""Creates a cellranger script for a multiome experiment (RNA + ATAC)
	
	Args:
	  gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
	  jobscript_path: path to save the jobscript
	Returns:
	  None 
	"""
	job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
	job_script = """#!/bin/bash
#SBATCH --job-name="{}"
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/{}_%x_%J.err
#SBATCH --output=./log/{}_%x_%J.out
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=12


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME

{} count --id {} --reference {} --libraries libraries.csv --localcores 24 --localmem 64;

echo [`date "+%Y-%m-%d %T"`] job finished""".format(gem_id, gem_id, gem_id, cfg.cellranger_path, gem_id, cfg.reference_path)
	job_script_file.write(job_script)
	job_script_file.close()
