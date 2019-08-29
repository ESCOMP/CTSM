
#================================================================================
# Makefile to compile the lilac  program
#================================================================================
## This is temporary Makefile for building lilac against CTSM pre-compiled library



#================================================================================
### Finding and including esmf.mk
#================================================================================

# Note: This fully portable Makefile template depends on finding environment
# #       variable "ESMFMKFILE" set to point to the appropriate "esmf.mk" file,
# #       as is discussed in the User's Guide.
# #       However, you can still use this Makefile template even if the person
# #       that installed ESMF on your system did not provide for a mechanism to
# #       automatically set the environment variable "ESMFMKFILE". In this case
# #       either manually set "ESMFMKFILE" in your environment or hard code the
# #       location of "esmf.mk" into the include statement below.
# #       Notice that the latter approach has negative impact on flexibility and 
# #       portability.


ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

#================================================================================
### Define directory paths
#================================================================================
# Temporarily hard-coded
# TODO: Please fix this part.
CTSM_BLD_DIR      =      /glade/scratch/negins/lilac_ctsm/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf
CTSM_INC          =    -I$(CTSM_BLD_DIR)/include
CTSM_LIB          =    -L$(CTSM_BLD_DIR)/lib -lclm
TRACEBACK_FLAGS   =    -g -traceback -debug all -check all
# -----------------------------------------------------------------------------
#EXTRA_LIBS       =    $(EXTRA_LIBS) -I/glade/scratch/negins/lilac_ctsm/bld/intel/mpt/nodebug/nothreads/nuopc/pio/pio2
EXTRA_LIBS        =   -L/glade/scratch/negins/lilac_ctsm/bld/intel/mpt/nodebug/nothreads/nuopc/nuopc/esmf/c1a1l1/lib -lcsm_share -L/glade/scratch/negins/lilac_ctsm/bld/intel/mpt/nodebug/nothreads/nuopc/lib -lpiof -lpioc -lgptl -lmct -lmpeu -mkl=cluster  -L/glade/u/apps/ch/opt/pnetcdf/1.11.0/mpt/2.19/intel/19.0.2//lib -lpnetcdf -L/glade/u/home/dunlap/ESMF-INSTALL/intel19/8.0.0bs32/lib/libO/Linux.intel.64.mpt.default -L/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib -L/glade/u/home/dunlap/YAML-INSTALL/lib  -Wl,-rpath,/glade/u/home/dunlap/ESMF-INSTALL/intel19/8.0.0bs32/lib/libO/Linux.intel.64.mpt.default  -Wl,-rpath,/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib  -Wl,-rpath,/glade/u/home/dunlap/YAML-INSTALL/lib -lesmf  -cxxlib -lrt -ldl -lnetcdff -lnetcdf -lyaml-cpp  -cxxlib
MORE_LIBSi        =   -I/glade/scratch/negins/lilac_ctsm/bld/intel/mpt/nodebug/nuopc/nuopc/esmf/c1a1l1/csm_share/ -L/glade/scratch/negins/lilac_ctsm/bld/intel/mpt/nodebug/nuopc/nuopc/esmf/c1a1l1/csm_share/
# -----------------------------------------------------------------------------


#================================================================================
### Compiler and linker rules using ESMF_ variables supplied by esmf.mk
#================================================================================

.SUFFIXES: .f90 .F90 .c .C

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \ 
		$(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
		$(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) \
		$(CTSM_INC) $(CTSM_LIB) $(TRACEBACK_FLAGS)           \
		$(EXTRA_LIBS) $(MORE_LIBS) $<

%.o : %.c
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) \ 
		$(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

% : %.C
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) \ 
		$(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

demo_driver: demo_driver.o atmos_cap.o lilac_mod.o lilac_utils.o cpl_mod.o lnd_cap.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(CTSM_INC) $(CTSM_LIB) $(EXTRA_LIBS) $(TRACEBACK_FLAGS) $(MORE_LIBS)
	mv demo_driver demo_driver.exe
	rm *.o *.mod

# module dependencies:
#demo_driver.o: lilac_mod.o atmos_cap.o lilac_utils.o cpl_mod.o shr_string_mod.o shr_kind_mod.o shr_sys_mod.o shr_log_mod.o shr_timer_mod.o
demo_driver.o: lilac_mod.o atmos_cap.o lilac_utils.o cpl_mod.o
lilac_mod.o: atmos_cap.o lilac_utils.o cpl_mod.o lnd_cap.o #shr_pio_mod.o
atmos_cap.o: lilac_utils.o

# -----------------------------------------------------------------------------

.PHONY: clean berzerk remake
clean:
	rm -f *.exe
berzerk:
	rm -f PET*.ESMF_LogFile job_name* *.o *.mod *.exe
remake:
	rm lilac_mod.o demo_driver.o demo_driver.exe & make
# -----------------------------------------------------------------------------
