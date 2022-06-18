#!/bin/bash

# Load required modules
module load PYTHON/2.7.5


# Get information for each library (flow cell, lane, sample id, etc.)
# $1  needs to be the name of the project
/scratch/project/production/DAT/apps/LIMSQ/limsq -sp $1 | sed 's/;/\t/g' > "projects/${1}/info.txt"
