paper.pdf: paper.tex MyBib.bib pic/*
	pdflatex -shell-escape -jobname paper paper.tex

.PHONY: quick
quick: .quick
.quick: paper.tex MyBib.bib pic/*
	touch .quick
	pdflatex -shell-escape -jobname paper paper.tex -o paper.pdf

.PHONY: 4
4: .4
.4: paper.tex MyBib.bib pic/*
	touch .4
	pdflatex -shell-escape -jobname paper paper.tex -o paper.pdf
	bibtex *.aux
	pdflatex -shell-escape -jobname paper paper.tex
	pdflatex -shell-escape -jobname paper paper.tex

.PHONY: dvi
dvi: .dvi
.dvi: paper.tex MyBib.bib
	touch .dvi
	latex -jobname paper paper.tex
	bibtex *.aux
	latex -jobname paper paper.tex
	latex -jobname paper paper.tex

clean:
	rm -f *.aux *.bbl *.blg *.log *.pdf *.spl *.dvi *.out .4 .dvi tmpbib.tex

clean-all:
	rm -f *.aux *.bbl *.blg *.log *.pdf *.spl *.dvi *.out .4 .dvi tmpbib.tex pic/*.pdf
