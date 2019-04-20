NimbleUserManual.pdf: *Rmd
	Rscript -e "bookdown::render_book('index.Rmd', output_format='bookdown::pdf_book')"

NimbleUserManual.html: *Rmd
	Rscript -e "bookdown::render_book('index.Rmd', output_format='bookdown::gitbook')"

%.Rmd: %.Rnw
	Rscript -e "source('latex2Rmd.R'); latex2Rmd(\"$(basename $(@))\")"

clean:
	rm ../bookdown-demo.pdf


