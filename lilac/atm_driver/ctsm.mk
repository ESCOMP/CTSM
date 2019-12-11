# ======================================================================
# Include this file to get makefile variables needed to include / link
# LILAC/CTSM in an atmosphere model's build
#
# Variables of interest are:
# - CTSM_INCLUDES: add this to the compilation line
# - CTSM_LIBS: add this to the link line
#
# Note: You must set the environment BLDDIR before running this - e.g.
# export BLDDIR=/glade/scratch/sacks/test_lilac_1205a/bld
#
# ESMFMKFILE must also be set in the environment
# ======================================================================

include $(ESMFMKFILE)

SHARED_BLD_DIR = $(BLDDIR)/intel/mpt/debug/nothreads/nuopc
CTSM_BLD_DIR   = $(SHARED_BLD_DIR)/nuopc/esmf
DEPENDS_LIB    = $(SHARED_BLD_DIR)/lib
SHR_LIB        = $(SHARED_BLD_DIR)/nuopc/esmf/c1a1l1/lib
CTSM_INC       = $(CTSM_BLD_DIR)/clm/obj

LIBS = -L$(CTSM_BLD_DIR)/lib -lclm -L$(SHR_LIB) -lcsm_share -L$(DEPENDS_LIB) -lpiof -lpioc -lgptl -lmct -lmpeu -mkl=cluster  -L/glade/u/apps/ch/opt/pnetcdf/1.11.0/mpt/2.19/intel/19.0.2//lib -lpnetcdf -L/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib -Wl,-rpath,/glade/u/apps/ch/opt/netcdf-mpi/4.6.1/mpt/2.19/intel/19.0.2/lib -cxxlib -lrt -ldl -lnetcdff -lnetcdf -cxxlib

CTSM_INCLUDES = $(ESMF_F90COMPILEPATHS) -I$(CTSM_INC)

CTSM_LIBS =  $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS) $(LIBS)
