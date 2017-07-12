NimbleUserManual.pdf :	NimbleUserManual.Rnw chapter_*.Rnw *Table*tex
	MAKEFLAGS=-j1 Rscript -e "library(nimble); library(knitr);  knit2pdf('NimbleUserManual.Rnw')"
# should we remove the compiled .o/.so in model1 too?
# loading library(nimble) externally to knitr prevents library-loading messages from getting into the final pdf, I think

includesMCMConly.pdf :	includesMCMConly.Rnw chapter_MCMC.Rnw
	Rscript -e "library(knitr);  knit2pdf('includesMCMConly.Rnw')"
	open -a "Google Chrome" includesMCMConly.pdf

NimbleUserManualNoEval.pdf :	NimbleUserManual.Rnw
	Rscript makeNimbleUserManualNoEval.R 
# makeNimbleUserManualNoEval.R doesn't seem to exist
# one option here is to programmatically (via sed, e.g)
# change all chunks to have eval=FALSE and then run knit/knit2pdf

clean :
	rm -f *.aux *.log *.toc *.out chapter_*.tex NimbleUserManual.tex

