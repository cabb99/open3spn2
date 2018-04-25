BASE=  erratum
TEX=   $(BASE).tex
#BIB=   bib/$(BASE).bib
PDF=   $(BASE).pdf

default: pdf
remake:  clean pdf

clean:
	rm -v ${BASE}.ps ${BASE}.pdf ${BASE}.dvi ${BASE}.aux ${BASE}.dvi ${BASE}.log ${BASE}.bbl ${BASE}.blg

pdf:    $(TEX) $(BIB) $(FIG)
	latex  $(TEX)      > texout
	bibtex    $(BASE).aux     > bibout
	latex  $(TEX)     >> texout
	latex  $(TEX)     >> texout
	dvips  $(BASE).dvi
	ps2pdf $(BASE).ps

svn:
	svn commit -m "updating $(BASE)"
	svn update
