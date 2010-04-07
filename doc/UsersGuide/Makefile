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
CO2DIF    := addco2_datm.buildnml
COMPLIST  := compsets_list_ofIcases.xml
SOURCES   := $(DOCBKUG) $(COMPLIST) $(CFGLOG).xml $(BNMLOG).xml \
             $(RESLOG).xml $(USCLOG).xml $(QCKLOG).xml \
             $(USRLOG).xml $(GETREG).xml preface.xml custom.xml special_cases.xml \
             tools.xml adding_files.xml single_point.xml addco2_datm.buildnml.xml

CONVAMP    := sed 's/\&/\&amp;/g'
CONVSIGNS  := sed 's/>/\&gt;/g' | sed 's/</\&lt;/g'

all: $(HTMLUG) $(PDFUG)

.SUFFIXES:
.SUFFIXES: .xml .pdf .html .tlog .diff

$(CFGLOG).xml: $(CFGLOG).tlog
$(BNMLOG).xml: $(BNMLOG).tlog
$(RESLOG).xml: $(RESLOG).tlog
$(USCLOG).xml: $(USCLOG).tlog
$(QCKLOG).xml: $(QCKLOG).tlog
$(USRLOG).xml: $(USRLOG).tlog
$(GETREG).xml: $(GETREG).tlog
$(CO2DIF).xml: $(CO2DIF).diff

.tlog.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > $@

.diff.xml:
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

$(BNMLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -help >& $@

$(RESLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -res list >& $@

$(USCLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -use_case list >& $@

$(CFGLOG).tlog: configure
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(CFGLOG) is good and redo your make"
	../../bld/configure -help >& $@

$(QCKLOG).tlog: Quickstart.GUIDE
	cp $<  $@

$(USRLOG).tlog: Quickstart.userdatasets
	cp $<  $@

$(GETREG).tlog:  getregional_datasets.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(GETREG) is good and redo your make"
	../../tools/ncl_scripts/getregional_datasets.pl -help >& $@

clean:
	rm -f $(HTMLUG) $(PDFUG)

realclean: clean
	rm -f f*.html c*.html x*.html a*.html i*.html $(COMPLIST) $(CFGLOG).xml $(BNMLOG).xml \
	$(RESLOG).xml $(USCLOG).xml $(USRLOG).xml $(GETREG).xml $(QCKLOG).xml \
	$(CO2DIF).xml *.tlog
