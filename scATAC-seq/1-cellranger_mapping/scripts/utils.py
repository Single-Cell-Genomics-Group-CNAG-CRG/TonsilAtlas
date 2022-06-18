# Python functions used in different scripts

import numpy as np
import pandas as pd
import subprocess
import os
import config_vars as cfg


def create_fastq_symlink(gem_id, fastq_path_df, symlink_path):
    """Creates a symbolic link pointing to a fastq file using cellranger notation
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      fastq_path_df: pandas dataframe with the fastq paths for that gem_id
      symlink_path: string specifying where to create the symlinks
    
    Returns:
      None 
    """
    pair_ids = np.unique(fastq_path_df["pair_id"])
    for i in range(len(pair_ids)):
        filt = (fastq_path_df["pair_id"] == pair_ids[i])
        pair_df = fastq_path_df.loc[filt, :]
        for j in pair_df.index:
            fastq_path = pair_df.loc[j, "fastq_path"]
            lane = str(i + 1)
            read = pair_df.loc[j, "read"]
            read = read.replace("R", "")
            subprocess.run(["ln", "-s", fastq_path, "{}/{}_S1_L00{}_R{}_001.fastq.gz".format(symlink_path, gem_id, lane, read)])



def make_cellranger(gem_id, jobscript_path, fastq_path):
    """Creates a cellranger-atac script for a GEM well 
    
    Args:
      gem_id: identifier of the Gelbeads-in-Emulsion (GEM) well that will be used as prefix in the symlink
      jobscript_path: path to save the jobscript
      fastq_path: path to the fastq files

    Returns:
      None 
    """
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    job_script = """#!/bin/bash
    

#SBATCH --job-name={}
#SBATCH --workdir=.
#SBATCH --mail-type=all
#SBATCH --mail-user=ramon.massoni@cnag.crg.eu
#SBATCH --error=./log/{}_%x_%J.err
#SBATCH --output=./log/{}_%x_%J.out
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=14


echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME 

{} count --fastqs {} --id {} --sample {} --localcores 24 --localmem 64 --reference {};

echo [`date "+%Y-%m-%d %T"`] job finished""".format(gem_id, gem_id, gem_id, cfg.cellranger_path, fastq_path, gem_id, gem_id, cfg.reference_path)
    job_script_file.write(job_script)
    job_script_file.close()


