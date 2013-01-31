TARGETS=paragraphs.pdf slides.pdf submission.pdf
all: $(TARGETS)

paragraphs.pdf: paragraphs.tex refs.bib
	pdflatex -interaction batchmode paragraphs.tex
	bibtex paragraphs.aux
	pdflatex -interaction batchmode paragraphs.tex
	pdflatex -interaction batchmode paragraphs.tex

submission.pdf: submission.tex refs.bib
	- pdflatex -interaction batchmode submission.tex
	bibtex submission.aux
	- pdflatex -interaction batchmode submission.tex
	- pdflatex -interaction batchmode submission.tex

slides.pdf: slides.tex refs.bib
	pdflatex -interaction batchmode slides.tex
	bibtex slides.aux
	pdflatex -interaction batchmode slides.tex
	pdflatex -interaction batchmode slides.tex

.PHONY: clean

clean:
	rm -rf *.log *.aux *.out *.bbl *.nav *.snm *.toc *.blg $(TARGETS)
