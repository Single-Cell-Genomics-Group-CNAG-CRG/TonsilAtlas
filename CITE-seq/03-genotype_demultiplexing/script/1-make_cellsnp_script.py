"""
Author: "Sonal Rashmi"
Date: 2021-05-18
Description: This script creates the directory and executable script structure for cellSNP (customised for cellranger 6, tonsil atlas project and to run on cluster)
Input : path to cellranger working directory and CellSNP working directory
"""

import sys
import pandas as pd 
import os
import argparse

def make_cellsnp_script(cellranger_run_path, subproject, gem_id, jobscript_path, min_maf = 0.1, min_count = 20):

	BAM=os.path.join(cellranger_run_path,"projects",subproject,"jobs",gem_id,gem_id,"outs/per_sample_outs",gem_id,"count/sample_alignments.bam")
	BARCODE=os.path.join(cellranger_run_path,"projects",subproject,"jobs",gem_id,gem_id,"outs/per_sample_outs",gem_id,"count/sample_feature_bc_matrix/barcodes.tsv.gz")
	REGION_VCF="/home/devel/srashmi/reference/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
	OUT_DIR=os.path.join(jobscript_path,gem_id)
	MIN_MAF = float(min_maf)
	MIN_COUNT = int(min_count)
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	log_dir = os.path.join(OUT_DIR, "log")
	if not os.path.exists(log_dir):
		os.makedirs(log_dir)
	job_script_file = open("{}/{}.cmd".format(OUT_DIR, gem_id), "w")	
	
	job_script = """#!/usr/bin/env bash
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --job-name={}
#SBATCH --workdir=.
#SBATCH --error=./log/{}.err
#SBATCH --output=./log/{}.out
echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME
source /scratch/groups/hheyn/software/anaconda3/bin/activate cellsnp
cellsnp-lite -s {} -b {} -O {} -R {} -p 20 --minMAF {} --minCOUNT {} --gzip

echo [`date "+%Y-%m-%d %T"`] job finished!!
""".format(gem_id,gem_id,gem_id,BAM,BARCODE,OUT_DIR,REGION_VCF, MIN_MAF, MIN_COUNT)
	job_script_file.write(job_script)
	job_script_file.close()
 

def main():
	parser = argparse.ArgumentParser(description = "options to initialize the filesystem and scripts for cellsnp run")
	parser.add_argument("--cellranger_output_path",
						dest = "cellranger_output_path",
						action = "store",
						default = None,
						help = "Successful run cellranger output")
	parser.add_argument("--cellsnp_output_dir",
						dest = "cellsnp_output_dir",
						action = "store",
						default = None,
						help = "Cellsnp Output Dir")
	parser.add_argument("--min_minor_allele_frequency",
						dest = "min_maf",
						action = "store",
						default = 0.1,
						type = float,
						help = "threshold for the minimum minor allele frequency required for a SNP to be considered")
	parser.add_argument("--min_reads",
						dest = "min_count",
						action = "store",
						default = 20,
						type = int,
						help = "minimum number of supporting reads or UMIs required for a SNP to be considered")
	if len(sys.argv) < 5:
		parser.print_help()
		sys.exit(0)
	options = parser.parse_args()
	cellranger_run_path = options.cellranger_output_path
	jobscript_path = options.cellsnp_output_dir	
	min_maf = options.min_maf
	min_count = options.min_count
	if not os.path.exists(cellranger_run_path):
		print("Cellranger output directory not found. Program terminated!")
		sys.exit(0)
	else:
		if not os.path.exists(jobscript_path):
			os.makedirs(jobscript_path) 
		metadata = pd.read_csv(os.path.join(cellranger_run_path,"data/tonsil_atlas_metadata.csv"))
		metadata_multiplexed = metadata[metadata['donor_id'] == "multiplexed"]
		iterable = list(sorted(set(zip(metadata_multiplexed["subproject"], metadata_multiplexed["gem_id"]))))
		for subproject, gem_id in iterable:
			make_cellsnp_script(cellranger_run_path, subproject, gem_id, jobscript_path, min_maf, min_count)

if __name__ =="__main__":
		main()			
			
