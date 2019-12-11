#================================================================================
# Makefile to compile atm_driver on cheyenne
#================================================================================

#================================================================================
# Note: You must set the environment variable CTSM_MKFILE before running this - e.g.
# export CTSM_MKFILE=/glade/scratch/sacks/test_lilac_1205a/bld/ctsm.mk
#
# ESMFMKFILE must also be set in the environment
#================================================================================

include $(CTSM_MKFILE)

# Most atmosphere model builds shouldn't need this directly, but we use
# it here in order to easily get a f90 compiler and f90 compile opts for
# building atm_driver.o. (This is a bit of a kludge that we should
# change later.)
include $(ESMFMKFILE)

FFLAGS  = -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source -qopt-report -xCORE_AVX2 -no-fma  -O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created -DLINUX  -DFORTRANUNDERSCORE -DCPRINTEL -DDEBUG -DHAVE_MPI -DPIO1 -DHAVE_SLASHPROC -D_PNETCDF -DESMF_VERSION_MAJOR=7 -DESMF_VERSION_MINOR=1 -free

#================================================================================
# Compiler and linker rules using ESMF_ variables supplied by esmf.mk
#================================================================================

.SUFFIXES: .F90

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(CTSM_INCLUDES) $(FFLAGS) $<

atm_driver.o : $(CURDIR)/atm_driver.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(CTSM_INCLUDES) $(FFLAGS) $<

atm_driver: atm_driver.o 
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) -o $@ $^ $(CTSM_LIBS)
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

