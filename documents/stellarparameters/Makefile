LATEX       = pdflatex
CHECK_RERUN = grep "Rerun to get" $*.log

all: parallax.pdf

parallax.pdf: parameters02.png

%.pdf: %.tex
	${LATEX} $<
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."
	( ${CHECK_RERUN} && ${LATEX} $< ) || echo "Done."

parameters02.png: parameters_pgm.py
	python3 parameters_pgm.py
