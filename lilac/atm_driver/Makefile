#================================================================================
# Makefile to compile atm_driver on cheyenne
#================================================================================

#================================================================================
# NOTE: Before running this, you must:
#
# (1) Run cime's configure tool in order to generate a Macros.make file
#
# (2) Source the .env_mach_specific.sh file created by the configure
#     tool in order to set up the environment correctly. Among other
#     things, this should set the environment variable ESMFMKFILE. (See
#     notes below about the need for this.)
#
# (3) Set the environment variable CTSM_MKFILE - e.g.
#
#     export CTSM_MKFILE=/glade/scratch/sacks/test_lilac_1205a/bld/ctsm.mk
#
#================================================================================

include Macros.make

include $(CTSM_MKFILE)

# Most atmosphere model builds shouldn't need this directly, but we use
# it here in order to easily get a f90 compiler and f90 compile opts for
# building atm_driver.o. (This is a bit of a kludge that we should
# change later.)
include $(ESMFMKFILE)

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

