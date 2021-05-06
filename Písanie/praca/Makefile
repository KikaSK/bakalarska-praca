all: main.pdf

main.pdf: main.tex *.tex *.bib images/*
	pdflatex main
	bibtex main
	pdflatex main
	pdflatex main

clear:
	rm -f main.pdf main.aux