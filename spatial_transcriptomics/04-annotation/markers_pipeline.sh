#!/usr/bin/env bash

# Script to be run from where the Rmarkdown document is
# set -euo pipefail
# 
# declare -A id_dict=( \
#   ["esvq52_nluss5"]="BCLL-10-T" \
#   ["c28w2r_7jne4i"]="BCLL-8-T" \
#   ["p7hv1g_tjgmyj"]="BCLL-12-T" \
#   ["tarwe1_xott6q"]="BCLL-2-T")
# 
# for id in esvq52_nluss5 c28w2r_7jne4i p7hv1g_tjgmyj tarwe1_xott6q
# do
# # Get the donor id in human form
# donor_id=${id_dict[${id}]}
# R -e "rmarkdown::render('markers_resolutions.Rmd',
#                         params = list( 
#                           sample_id = '$id',
#                           donor_id = '$donor_id'),
#                         output_file='markers_resolutions_${donor_id}.html')"
# done


R -e "rmarkdown::render('04-DE_markers_integrated.Rmd',
                        output_file='04-DE_markers_integrated.html')"
