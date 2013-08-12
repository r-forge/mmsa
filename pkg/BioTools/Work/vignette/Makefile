all: tex copy

tex:
	R CMD Sweave BioTools.Rnw
	texi2dvi --pdf BioTools.tex

copy:
	cp BioTools.pdf ../../vignettes
