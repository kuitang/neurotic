TARGETS=paragraphs.pdf slides.pdf submission.pdf
all: $(TARGETS)

paragraphs.pdf: paragraphs.tex refs.bib
	pdflatex paragraphs.tex
	bibtex paragraphs.aux
	pdflatex paragraphs.tex
	pdflatex paragraphs.tex

submission.pdf: submission.tex refs.bib
	- pdflatex submission.tex
	bibtex submission.aux
	- pdflatex submission.tex
	- pdflatex submission.tex

slides.pdf: slides.tex refs.bib
	pdflatex slides.tex
	bibtex slides.aux
	pdflatex slides.tex
	pdflatex slides.tex

.PHONY: clean

clean:
	rm -rf *.log *.aux *.out *.bbl *.nav *.snm *.toc *.blg $(TARGETS)
