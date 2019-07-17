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

ifneq ($(origin ESMFMKFILE), environment)
$(error Environment variable ESMFMKFILE was not set.)
endif

include $(ESMFMKFILE)

################################################################################
################################################################################

.SUFFIXES: .f90 .F90 .c .C

%.o : %.f90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREENOCPP) $<

%.o : %.F90
	$(ESMF_F90COMPILER) -c $(ESMF_F90COMPILEOPTS) $(ESMF_F90COMPILEPATHS) $(ESMF_F90COMPILEFREECPP) $(ESMF_F90COMPILECPPFLAGS) $<
        
%.o : %.c
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<

% : %.C
	$(ESMF_CXXCOMPILER) -c $(ESMF_CXXCOMPILEOPTS) $(ESMF_CXXCOMPILEPATHSLOCAL) $(ESMF_CXXCOMPILEPATHS) $(ESMF_CXXCOMPILECPPFLAGS) $<


# -----------------------------------------------------------------------------
demo_driver: demo_driver.o atmos_cap.o lilac_mod.o lilac_utils.o cpl_mod.o lnd_cap.o lnd_comp_esmf.o
	$(ESMF_F90LINKER) $(ESMF_F90LINKOPTS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) -o $@ $^ $(ESMF_F90ESMFLINKLIBS)
	mv demo_driver demo_driver.exe
	#rm *.o *.mod

# module dependencies:
demo_driver.o: lilac_mod.o atmos_cap.o lilac_utils.o cpl_mod.o shr_string_mod.o shr_kind_mod.o shr_sys_mod.o shr_log_mod.o shr_timer_mod.o
#demo_driver.o: lilac_mod.o atmos_cap.o lilac_utils.o cpl_mod.o shr_string_mod.o
lilac_mod.o: atmos_cap.o lilac_utils.o cpl_mod.o lnd_cap.o lnd_comp_esmf.o
#lilac_mod.o: atmos_cap.o lilac_utils.o cpl_mod.o lnd_cap.o
shr_string_mod.o: shr_kind_mod.o shr_sys_mod.o shr_log_mod.o shr_timer_mod.o
shr_log_mod.o:
atmos_cap.o: lilac_utils.o
#lnd_comp_esmf.o: shr_string_mod.o abort_utils.o
shr_sys_mod.o: shr_log_mod.o
abort_utils.o: 
clm_varctl.o:
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
.PHONY: dust clean distclean berzerk
dust:
	rm -f PET*.ESMF_LogFile DE.nc FLAG.nc OMEGA.nc SIE.nc U_velocity.nc V_velocity.nc
clean:
	rm -f *.exe *.o *.mod
distclean: dust clean
berzerk:
	rm -f PET*.ESMF_LogFile job_name* *.o *.mod *.exe
