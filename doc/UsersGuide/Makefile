#
# Makefile to convert DocBook CLM Users-Guide into html and/or pdf
# (rtf, txt, ps, tex, man, dvi, and texi are also valid docbook formats)
#
.SUFFIXES: .xml .pdf .html

VPATH := . ../../bld

PDFUG   := index.pdf
HTMLUG  := book1.html
DOCBKUG := index.xml
CFGLOG  := config_help.log
BNMLOG  := buildnml_help.log

all: $(HTMLUG) $(PDFUG)

$(HTMLUG): $(DOCBKUG)
	docbook2html $<

$(PDFUG): $(DOCBKUG)
	docbook2pdf $<

$(BNMLOG): build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(BNMLOG) is good and redo your make"
	../../bld/build-namelist -help >& $(BNMLOG)

$(CFGLOG): configure
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(CFGLOG) is good and redo your make"
	../../bld/configure -help >& $(CFGLOG)
clean:
	rm -f $(HTMLUG) $(PDFUG)
realclean: clean
	rm -f f*.html c*.html x*.html
