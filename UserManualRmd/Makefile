bookdown-demo.html: *Rmd
	Rscript -e "bookdown::render_book('index.Rmd', output_format='bookdown::pdf_book')"

%.Rmd: %.Rnw
	Rscript -e "source('latex2Rmd.R'); latex2Rmd(\"$(basename $(@))\")"

clean:
	rm bookdown-demo.pdf


