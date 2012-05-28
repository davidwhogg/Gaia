all: attitude_model.pdf

.SUFFIXES: .tex .dvi .ps .pdf .mp .eps

%.pdf: %.tex
	pdflatex $<
	pdflatex $<
	pdflatex $<

clean:
	\rm -vf *.log *.aux
