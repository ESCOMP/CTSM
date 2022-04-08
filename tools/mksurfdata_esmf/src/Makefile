# GNU Makefile template for user ESMF application

################################################################################
################################################################################
## This Makefile must be able to find the "esmf.mk" Makefile fragment in the  ##
## 'include' line below. Following the ESMF User's Guide, a complete ESMF     ##
## installation should ensure that a single environment variable "ESMFMKFILE" ##
## is made available on the system. This variable should point to the         ##
## "esmf.mk" file.                                                            ##
##                                                                            ##
## This example Makefile uses the "ESMFMKFILE" environment variable.          ##
##                                                                            ##
## If you notice that this Makefile cannot find variable ESMFMKFILE then      ##
## please contact the person responsible for the ESMF installation on your    ##
## system.                                                                    ##
## As a work-around you can simply hardcode the path to "esmf.mk" in the      ##
## include line below. However, doing so will render this Makefile a lot less ##
## flexible and non-portable.                                                 ##
################################################################################


COMPILEPATH = -I$(MKSURFDATA_BLD) -I/glade/u/apps/ch/opt/pio/2.5.5/mpt/2.22/intel/19.1.1/include
LINKPATH = -L$(MKSURFDATA_BLD)
PIO_LINKPATH = -L/glade/work/jedwards/tools/pio/2.5.5/intel/19.1.1/mpt/2.22/lib -L/glade/u/apps/ch/opt/pnetcdf/1.12.2/mpt/2.22/intel/19.1.1//lib -L/glade/u/apps/ch/opt/netcdf/4.7.3/intel/19.0.5/lib

ESMFMKFILE=/glade/work/oehmke/ESMF/scalable_mesh_from_file/lib/libg/Linux.intel.64.mpt.default/esmf.mk
include $(ESMFMKFILE)

################################################################################
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $(COMPILEPATH)  $<

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(COMPILEPATH) $<

%.o : %.c
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

%.o : %.C
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

# -----------------------------------------------------------------------------
mksurfdata: mksurfdata.o
	 $(ESMF_F90LINKER) -g $(LINKPATH) $(PIO_LINKPATH) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) -l mksurfdata -l piof -l pnetcdf -lnetcdf

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean info edit
dust:
	rm -f PET*.ESMF_LogFile *.nc *.stdout
clean:
	rm -f esmApp *.o *.mod
distclean: dust clean

info:
	@echo ==================================================================
	@echo ESMFMKFILE=$(ESMFMKFILE)
	@echo ==================================================================
	@cat $(ESMFMKFILE)
	@echo ==================================================================

run:
	$(ESMF_INTERNAL_MPIRUN) -np 4 ./esmApp
