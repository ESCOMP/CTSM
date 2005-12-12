.SUFFIXES : .tex .dvi .ps .pdf

all : CodeReference

#-----------------------------------------------------------------------
# Time manager design

CodeReference : CodeReference.ps CodeReference.pdf CodeReference/CodeReference.html

CodeReference/CodeReference.html : CodeReference.dvi
	latex2html  -split +1 -show_section_numbers CodeReference.tex
	cp clm_offline_tree.html CodeReference/.

CodeReference.dvi : CodeReference.tex
	latex $<
	latex $<

#-----------------------------------------------------------------------
.tex.dvi :
	latex $<
	latex $<

.dvi.ps :
	dvips -o $@ $<

.dvi.pdf :
	dvipdf $< $@

clean:
	rm -f *.aux *.bbl *.blg *.dvi *.log *.toc

realclean: clean
	rm -rf CodeReference  *.pdf *.ps
