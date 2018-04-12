all: gaia_like.pdf

.SUFFIXES: .tex .dvi .ps .pdf .mp .eps

%.pdf: %.tex
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

clean:
	\rm -vf *.log *.aux
