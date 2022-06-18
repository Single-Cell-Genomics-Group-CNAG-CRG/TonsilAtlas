R -e "rmarkdown::render('MAGIC-denoising-GC.Rmd',
                        output_file='MAGIC-denoising-GC.html')"

R -e "rmarkdown::render('MAGIC-visualization-GC.Rmd',
                        output_file='MAGIC-visualization-GC.html')"

R -e "rmarkdown::render('GC-deconvolution.Rmd',
                        output_file='GC-deconvolution.html')"

R -e "rmarkdown::render('GC-deconvolution-visualization.Rmd',
                        output_file='GC-deconvolution-visualization.html')"
