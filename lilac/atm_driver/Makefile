
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


#ifneq ($(origin ESMFMKFILE), environment)
#$(error Environment variable ESMFMKFILE was not set.)
#endif

ESMFMKFILE =  /glade/work/dunlap/ESMF-INSTALL/8.0.0bs38/lib/libO/Linux.intel.64.mpt.default/esmf.mk
include $(ESMFMKFILE)

#================================================================================
### Define directory paths
#================================================================================
# Temporarily hard-coded
# TODO: Please fix this part.
CASE_NAME        =  why01-g
CTSM_BLD_DIR     =  /glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/nothreads/mct/mct/esmf
CTSM_INC         =  -I$(CTSM_BLD_DIR)/include
CTSM_LIB         =  -L$(CTSM_BLD_DIR)/lib -lclm
#TRACEBACK_FLAGS =  -g -traceback -debug all -check all -O2 -r8
TRACEBACK_FLAGS  = -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -qopt-report -xCORE_AVX2 -no-fma  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created -DLINUX  -DFORTRANUNDERSCORE -DCPRINTEL -DDEBUG -DHAVE_MPI -DPIO1 -DHAVE_SLASHPROC -D_PNETCDF -DESMF_VERSION_MAJOR=7 -DESMF_VERSION_MINOR=1 -free

# -----------------------------------------------------------------------------
EXTRA_LIBS  =   -L/glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/nothreads/mct/mct/esmf/c1a1l1i1o1r1g1w1i1e1/lib -lcsm_share -L/glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/nothreads/mct/lib -lpio -lgptl -lmct -lmpeu -mkl=cluster  -L/glade/u/apps/ch/opt/pnetcdf/1.11.0/mpt/2.19/intel/19.0.2//lib -lpnetcdf -L/glade/u/home/dunlap/ESMF-INSTALL/intel19/8.0.0bs32/lib/libO/Linux.intel.64.mpt.default -L/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib -L/glade/u/home/dunlap/YAML-INSTALL/lib  -Wl,-rpath,/glade/u/home/dunlap/ESMF-INSTALL/intel19/8.0.0bs32/lib/libO/Linux.intel.64.mpt.default  -Wl,-rpath,/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib  -Wl,-rpath,/glade/u/home/dunlap/YAML-INSTALL/lib -lesmf  -cxxlib -lrt -ldl -lnetcdff -lnetcdf -lyaml-cpp  -cxxlib
MORE_LIBS  =   -I/glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/mct/mct/esmf/c1a1l1i1o1r1g1w1i1e1/csm_share/ -L/glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/mct/mct/esmf/c1a1l1i1o1r1g1w1i1e1/csm_share/     -I/glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/nothreads/mct/mct/esmf/clm/obj/ -I//glade/scratch/negins/$(CASE_NAME)/bld/intel/mpt/debug/nothreads/mct/mct/esmf/c1a1l1i1o1r1g1w1i1e1/csm_share/
# -----------------------------------------------------------------------------

DRIVER_DIR = $(CURDIR)
LILAC_DIR = $(DRIVER_DIR)/../lilac

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

lilac_atmcap.o : $(LILAC_DIR)/lilac_atmcap.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
		$(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(ESMF_F90COMPILEFREENOCPP) \
		$(CTSM_INC) $(CTSM_LIB) $(TRACEBACK_FLAGS)           \
		$(EXTRA_LIBS) $(MORE_LIBS) $<

lilac_cpl.o : $(LILAC_DIR)/lilac_cpl.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
		$(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(ESMF_F90COMPILEFREENOCPP) \
		$(CTSM_INC) $(CTSM_LIB) $(TRACEBACK_FLAGS)           \
		$(EXTRA_LIBS) $(MORE_LIBS) $<

lilac_utils.o : $(LILAC_DIR)/lilac_utils.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
		$(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(ESMF_F90COMPILEFREENOCPP) \
		$(CTSM_INC) $(CTSM_LIB) $(TRACEBACK_FLAGS)           \
		$(EXTRA_LIBS) $(MORE_LIBS) $<

lilac_mod.o : $(LILAC_DIR)/lilac_mod.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
		$(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(ESMF_F90COMPILEFREENOCPP) \
		$(CTSM_INC) $(CTSM_LIB) $(TRACEBACK_FLAGS)           \
		$(EXTRA_LIBS) $(MORE_LIBS) $<

atm_driver.o : $(DRIVER_DIR)/atm_driver.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) \
		$(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $(ESMF_F90COMPILEFREENOCPP) \
		$(CTSM_INC) $(CTSM_LIB) $(TRACEBACK_FLAGS)           \
		$(EXTRA_LIBS) $(MORE_LIBS) $<

atm_driver: atm_driver.o lilac_atmcap.o lilac_mod.o lilac_utils.o lilac_cpl.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(CTSM_INC) $(CTSM_LIB) $(EXTRA_LIBS) $(TRACEBACK_FLAGS) $(MORE_LIBS)
	mv atm_driver atm_driver.exe
	rm *.o *.mod

# module dependencies:
#atm_driver.o: lilac_mod.o lilac_atmcap.o lilac_utils.o lilac_cpl.o shr_string_mod.o shr_kind_mod.o shr_sys_mod.o shr_log_mod.o shr_timer_mod.o
atm_driver.o: lilac_mod.o lilac_atmcap.o lilac_utils.o lilac_cpl.o
lilac_mod.o: lilac_atmcap.o lilac_utils.o lilac_cpl.o  #shr_pio_mod.o
lilac_atmcap.o: lilac_utils.o

# -----------------------------------------------------------------------------

.PHONY: clean berzerk remake
clean:
	rm -f *.exe *.o
berzerk:
	rm -f PET*.ESMF_LogFile job_name* *.o *.mod *.exe
remake:
	rm lilac_mod.o atm_driver.o atm_driver.exe & make
# -----------------------------------------------------------------------------
