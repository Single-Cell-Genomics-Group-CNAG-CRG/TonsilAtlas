R -e "rmarkdown::render('MAGIC_denoising-plasma.Rmd',  output_file='MAGIC_denoising-plasma.html')"

R -e "rmarkdown::render('MAGIC_visualization-plasma.Rmd',  output_file='MAGIC_visualization-plasma.html')"

R -e "rmarkdown::render('plasma-deconvolution.Rmd',  output_file='plasma-deconvolution.html')"

R -e "rmarkdown::render('plasma-deconvolution-visualization.Rmd',  output_file='plasma-deconvolution-visualization.html')"

R -e "rmarkdown::render('PC-figure-plots.Rmd',  output_file='PC-figure-plots.html')"
