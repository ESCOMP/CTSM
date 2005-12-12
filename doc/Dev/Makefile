.SUFFIXES : .tex .dvi .ps .pdf

all : time-manager

#-----------------------------------------------------------------------
# Time manager design

time-manager : time-manager.ps time-manager.pdf time-manager/time-manager.html

time-manager/time-manager.html : time-manager.dvi
	latex2html  -split +1 -show_section_numbers time-manager.tex

time-manager.dvi : time-manager.tex
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
	rm -rf time-manager *.pdf *.ps
