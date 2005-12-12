#-----------------------------------------------------------------------
# This Makefile is for building any of the executables in this directory.
# USAGE: eg. gmake EXENAME=convert_lanwat
#------------------------------------------------------------------------

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

# The -Otask0 flag disables autotasking.  This is important due to the
# brain-damaged nature of local memory allocation on Cray PVP machines.

ifeq ($(UNAMEM),CRAY)
LIB_NETCDF := /usr/local/lib
INC_NETCDF := /usr/local/include
FPPFLAGS := -DCRAY -N -P $(cpp_path)
FFLAGS   := -DCRAY -c -G1 -Otask0 
LDFLAGS  := -L $(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
LIB_NETCDF := /usr/local/lib64/r4i4
INC_NETCDF := /usr/local/include
FFLAGS     := -DSGI -64 -r8 -i4 -O2 -c -I$(INC_NETCDF)
LDFLAGS    := -64 -L$(LIB_NETCDF) -lnetcdf 
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
LIB_NETCDF := /contrib/lib
INC_NETCDF := /contrib/include
FFLAGS     := -DSUN -c -O2 -stackvar -I$(INC_NETCDF)
LDFLAGS    := -L $(LIB_NETCDF) -lnetcdf -lc
endif

#------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------

FC := f90

.SUFFIXES:
.SUFFIXES: .F90 .o 

.F90.o:
	$(FC) $(FFLAGS) $<

$(EXENAME): $(EXENAME).o
	$(FC) -o $@ $(EXENAME).o $(LDFLAGS)

clean:
	$(RM) -f *.o $(EXENAME)










