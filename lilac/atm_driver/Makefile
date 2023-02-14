#================================================================================
# Makefile to compile atm_driver on cheyenne
#================================================================================

#================================================================================
# NOTE: Before running this, you must:
#
# (1) Run cime's configure tool in order to generate a Macros.make file
#
# (2) Source the .env_mach_specific.sh file created by the configure
#     tool in order to set up the environment correctly.
#
# (3) Set the environment variable CTSM_MKFILE - e.g.
#
#     export CTSM_MKFILE=/glade/scratch/sacks/test_lilac_1205a/bld/ctsm.mk
#
#================================================================================

include Macros.make

include $(CTSM_MKFILE)

.SUFFIXES: .F90

%.o : %.F90
	$(MPIFC) -c $(CTSM_INCLUDES) $(FFLAGS) $<

atm_driver.o : $(CURDIR)/atm_driver.F90
	$(MPIFC) -c $(CTSM_INCLUDES) $(FFLAGS) $<

atm_driver: atm_driver.o 
	$(MPIFC) -o $@ $^ $(LDFLAGS) $(CTSM_LIBS)
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

