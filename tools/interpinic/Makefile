# Makefile to build interpinic on various platforms
# Note: If netcdf library is not built in the standard location, 
# you must set the environment variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = interpinic
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .o

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
   LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
   INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -c -trapuv -I$(INC_NETCDF) -g -C -DEBUG:trap_uninitialized=ON
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF) -g
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
LIB_NETCDF := /usr/local/apps/netcdf-3.5/lib32/r4i4
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF) -qsuffix=f=f90 -g -qfullpath
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC = pgf90
FFLAGS =  -c -I$(INC_NETCDF) -fast
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := shr_kind_mod.o interpinic.o fmain.o addglobal.o wrap_nf.o 

.f90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

addglobal.o: addglobal.f90
interpinic.o : shr_kind_mod.o interpinic.f90
fmain.o : interpinic.o fmain.f90
shr_kind_mod.o: shr_kind_mod.f90
wrap_nf.o: shr_kind_mod.o wrap_nf.f90 
