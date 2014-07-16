all: tex copy

tex:
	R CMD Sweave BiostringsTools.Rnw
	texi2dvi --pdf BiostringsTools.tex

copy:
	cp BiostringsTools.pdf ../../vignettes/BTools_vignette.pdf
