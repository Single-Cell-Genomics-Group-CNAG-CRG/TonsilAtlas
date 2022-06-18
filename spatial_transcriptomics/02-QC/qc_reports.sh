#!/usr/bin/env bash
# Script to be run from where the Rmarkdown document is
set -euo pipefail

R -e "rmarkdown::render('01-metric_summary.Rmd',
                        output_file='01-metric_summary.html')"

R -e "rmarkdown::render('02-QC_common.Rmd',
                        output_file='02-QC_common.html')"


# declare -A id_dict=( \
#   ["tarwe1_xott6q"]="BCLL-2-T" \
#   ["c28w2r_7jne4i"]="BCLL-8-T" \
#   ["esvq52_nluss5"]="BCLL-10-T" \
#   ["p7hv1g_tjgmyj"]="BCLL-12-T" \
#   ["gcyl7c_cec61b"]="BCLL_13-T" \
#   ["zrt7gl_lhyyar"]="BCLL_14-T" \
#   ["qvwc8t_2vsr67"]="BCLL_9-T" \
#   ["exvyh1_66caqq"]="BCLL_11-T")

## SPOTlight visualization
# for id in ge14g1_8z5p2c gcvpj3_12dt1y bwdgo4_8matxx gqdsmq_vsrozg myhqei_7edxbq dlmtnt_p8fghv mw4d5x_pbn1ug jvajpo_5aqtp9
# do
# Get the donor id in human form
# donor_id=${id_dict[${id}]}
# R -e "rmarkdown::render('02-QC_common.Rmd',
#                         output_file='02-QC_common.html')"
# done