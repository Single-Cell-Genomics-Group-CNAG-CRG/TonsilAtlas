PATH_OBJECTS=../../results/R_objects/level_2/
cell_types=$(ls $PATH_OBJECTS)
for cell_type in $cell_types
do
Rscript -e "rmarkdown::render('03-visualize_integration.Rmd',params=list(cell_type='$cell_type'), output_file = paste('reports/03-','$cell_type', '_visualize_integration.html', sep = ''))"
done