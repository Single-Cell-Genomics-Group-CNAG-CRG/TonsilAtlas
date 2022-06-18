#!/bin/bash

# Author Juan Luis Trincado
module load PYTHON/2.7.5
module load lims/1.2

# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project
/scratch/production/DAT/apps/LIMSQ/limsq -sp $1 | sed 's/;/\t/g' > info.txt
