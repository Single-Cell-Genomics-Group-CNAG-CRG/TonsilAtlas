# Python functions used in different scripts

import numpy as np
import pandas as pd
import subprocess
import os
# import config_vars as cfg


def make_spaceranger(jobscript_path, metadata_df, spaceranger, reference_path):
    """Creates a cellranger script for a hashed GEM well 
    
    Args:
      jobscript_path: path to save the jobscript
      metadata_df: pandas object with 1 row containing the info for that gem_id of interest
      spaceranger: Path to spaceranger software - /scratch/groups/hheyn/software/spaceranger/1.1.0/
      reference_path: Path to reference transcriptome on which to map the reads - /scratch/groups/hheyn/data/reference/refdata-gex-GRCh38-2020-A

    Returns:
      None 
    """
    
    # Extract variables
    subproject = metadata_df["subproject"].values[0]
    gem_id = metadata_df["gem_id"].values[0]
    slide = metadata_df["slide"].values[0]
    area = metadata_df["area"].values[0]
    
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    print("Writing spaceranger script for {} in {}".format(gem_id, jobscript_path))
    job_script = """#!/bin/bash
    
    

module load PYTHON/2.7.5 
module load lims/1.2 

{} count --id={} \
         --transcriptome={} \
         --fastqs=../projects/{}/fastq \
         --sample={} \
         --image=../img/{}_{}_{}.jpg \
         --slide={} \
         --area={} \
         --localcores=8 \
         --localmem=64 \
         --slidefile=../data/spot_layout/{}.gpr

""".format(spaceranger, gem_id, reference_path, subproject, gem_id, gem_id, slide, area, slide, area, slide)
    # Write jobscript to file
    job_script_file.write(job_script)
    job_script_file.close()


def make_spaceranger_clustermode(jobscript_path, metadata_df, spaceranger, reference_path):
    """Creates a cellranger script for a hashed GEM well 
    
    Args:
      jobscript_path: path to save the jobscript
      metadata_df: pandas object with 1 row containing the info for that gem_id of interest
      spaceranger: Path to spaceranger software - /scratch/devel/melosua/phd/10x_software/spaceranger-1.1.0
      reference_path: Path to reference transcriptome on which to map the reads

    Returns:
      None 
    """
    
    # Extract variables
    subproject = metadata_df["subproject"].values[0]
    gem_id = metadata_df["gem_id"].values[0]
    slide = metadata_df["slide"].values[0]
    area = metadata_df["area"].values[0]
    
    job_script_file = open("{}/{}.cmd".format(jobscript_path, gem_id), "w")
    print("Writing spaceranger script for {} in {}".format(gem_id, jobscript_path))
    job_script = """#!/bin/bash
    
    

module load PYTHON/2.7.5 
module load lims/1.2 

{}/spaceranger count --id={} \
         --transcriptome={} \
         --fastqs=../projects/{}/fastq \
         --sample={} \
         --image=../img/{}_{}_{}.jpg \
         --slide={} \
         --area={} \
         --localcores=2 \
         --slidefile=../data/spot_layout/{}.gpr \
	 --jobmode={}/external/martian/jobmanagers/slurm.template

""".format(spaceranger, gem_id, reference_path, subproject, gem_id, gem_id, slide, area, slide, area, slide, spaceranger)
    # Write jobscript to file
    job_script_file.write(job_script)
    job_script_file.close()
