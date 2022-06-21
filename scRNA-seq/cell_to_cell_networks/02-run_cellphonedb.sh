#run cellphoneDB v3.0
META=path/to/CPDB/MDcelllabels.txt
COUNTS=/path/to/CPDB/input/data
DEG=path/to/CPDB/input/data/TableDEG_allvsRest_Wilcox.txt
OUTPUT_PATH=/path/to/CPDB/output
cellphonedb \
method degs_analysis \
$META $COUNTS $DEG \
--output-path=$OUTPUT_PATH \
--counts-data=hgnc_symbol