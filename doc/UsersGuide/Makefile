#
# Makefile to convert DocBook CLM Users-Guide into html and/or pdf
# (rtf, txt, ps, tex, man, dvi, and texi are also valid docbook formats)
#
VPATH := . .. ../../bld ../../tools/ncl_scripts

PDFUG     := clm_ug.pdf
HTMLUG    := book1.html
DOCBKUG   := clm_ug.xml
CFGLOG    := config_help
BNMLOG    := buildnml_help
RESLOG    := buildnml_resolutions
USCLOG    := buildnml_usecases
QCKLOG    := quickstart_guide
USRLOG    := quickstart_usrdat
GETREG    := getregional_datasets
COMPLIST  := compsets_list_ofIcases.xml
SOURCES   := $(DOCBKUG) $(COMPLIST) $(CFGLOG).xml $(BNMLOG).xml \
             $(RESLOG).xml $(USCLOG).xml $(QCKLOG).xml \
             $(USRLOG).xml $(GETREG).xml preface.xml custom.xml special_cases.xml \
             tools.xml adding_files.xml single_point.xml

CONVAMP    := sed 's/</\&amp;/g'
CONVSIGNS  := sed 's/>/\&gt;/g' | sed 's/&/\&lt;/g'

all: $(HTMLUG) $(PDFUG)

.SUFFIXES:
.SUFFIXES: .xml .pdf .html .log

$(CFGLOG).xml: $(CFGLOG).log
$(BNMLOG).xml: $(BNMLOG).log
$(RESLOG).xml: $(RESLOG).log
$(USCLOG).xml: $(USCLOG).log
$(QCKLOG).xml: $(QCKLOG).log
$(USRLOG).xml: $(USRLOG).log
$(GETREG).xml: $(GETREG).log

clm_ug.xml: preface.xml custom.xml special_cases.xml \
            tools.xml adding_files.xml single_point.xml

custom.xml: $(CFGLOG).xml $(RESLOG).xml $(BNMLOG).xml $(USCLOG).xml $(COMPLIST)
preface.xml: $(QCKLOG).xml
single_point.xml: $(USRLOG).xml $(GETREG).xml

.log.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > $@

debug:
	@echo "SOURCES:   $(SOURCES)"
	@echo "CONVAMP:   $(CONVAMP)"
	@echo "CONVSIGNS: $(CONVSIGNS)"

$(COMPLIST): 
	get_Icaselist.pl > $@

$(HTMLUG): $(SOURCES)
	docbook2html $<

$(PDFUG): $(SOURCES)
	docbook2pdf $<

$(BNMLOG).log: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -help >& $@

$(RESLOG).log: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -res list >& $@

$(USCLOG).log: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -use_case list >& $@

$(CFGLOG).log: configure
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(CFGLOG) is good and redo your make"
	../../bld/configure -help >& $@

$(QCKLOG).log: Quickstart.GUIDE
	cp $<  $@

$(USRLOG).log: Quickstart.userdatasets
	cp $<  $@

$(GETREG).log:  getregional_datasets.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(GETREG) is good and redo your make"
	../../tools/ncl_scripts/getregional_datasets.pl -help >& $@

clean:
	rm -f $(HTMLUG) $(PDFUG) *.log

realclean: clean
	rm -f f*.html c*.html x*.html $(COMPLIST) $(CFGLOG).xml $(BNMLOG).xml \
	$(RESLOG).xml $(USCLOG).xml $(USRLOG).xml $(GETREG).xml $(QCKLOG).xml
