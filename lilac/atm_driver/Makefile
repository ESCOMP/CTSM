#================================================================================
# Makefile to compile atm_driver on cheyenne
#================================================================================

#================================================================================
# Define directory paths
#================================================================================

ESMFMKFILE   = /glade/work/turuncu/ESMF/8.0.0b50/lib/libg/Linux.intel.64.mpt.default/esmf.mk
ESMF_LIB_DIR = /glade/work/turuncu/ESMF/8.0.0b50/lib/libg/Linux.intel.64.mpt.default
include $(ESMFMKFILE)

BLD_DIR       =  /glade/scratch/mvertens/test_lilac/bld

CTSM_BLD_DIR  = $(BLD_DIR)/intel/mpt/debug/nothreads/mct/mct/esmf
MCT_LIB       = $(BLD_DIR)/intel/mpt/debug/nothreads/mct/lib
SHR_LIB       = $(BLD_DIR)/intel/mpt/debug/nothreads/mct/mct/esmf/c1a1l1i1o1r1g1w1i1e1/lib
SHR_INC       = $(BLD_DIR)/intel/mpt/debug/mct/mct/esmf/c1a1l1i1o1r1g1w1i1e1/csm_share/
CTSM_INC      = $(CTSM_BLD_DIR)/clm/obj

FFLAGS  = -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -qopt-report -xCORE_AVX2 -no-fma  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created -DLINUX  -DFORTRANUNDERSCORE -DCPRINTEL -DDEBUG -DHAVE_MPI -DPIO1 -DHAVE_SLASHPROC -D_PNETCDF -DESMF_VERSION_MAJOR=7 -DESMF_VERSION_MINOR=1 -free

LIBS = -L$(CTSM_BLD_DIR)/lib -lclm -L$(SHR_LIB) -lcsm_share -L$(MCT_LIB) -lpio -lgptl -lmct -lmpeu -mkl=cluster  -L/glade/u/apps/ch/opt/pnetcdf/1.11.0/mpt/2.19/intel/19.0.2//lib -lpnetcdf -L$(ESMF_LIB_DIR) -L/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib  -Wl,-rpath,$(ESMF_LIB_DIR) -Wl,-rpath,/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib  -lesmf  -cxxlib -lrt -ldl -lnetcdff -lnetcdf -cxxlib

INCLUDES  = -I$(CTSM_BLD_DIR)/include -I$(SHR_INC) -I$(CTSM_INC)

#================================================================================
# Compiler and linker rules using ESMF_ variables supplied by esmf.mk
#================================================================================

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(INCLUDES) $(FFLAGS) $<

atm_driver.o : $(CURDIR)/atm_driver.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(INCLUDES) $(FFLAGS) $<

atm_driver: atm_driver.o 
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS) $(LIBS)
	mv atm_driver atm_driver.exe

# module dependencies:
atm_driver.o: 

.PHONY: clean berzerk remake
clean:
	rm -f *.exe *.o *.mod *.optr*
berzerk:
	rm -f PET*.ESMF_LogFile job_name* *.o *.mod *.exe
remake:
	rm lilac_mod.o atm_driver.o atm_driver.exe & make

