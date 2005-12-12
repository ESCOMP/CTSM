#-----------------------------------------------------------------------
# This Makefile is for building cprlndnc
# on a CRAY, SGI or SUN platform
#
#------------------------------------------------------------------------
# Set up special characters
null  :=


#
# Default variables
#

# Check for the NetCDF library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

EXENAME := cprlndnc
RM      := rm

UNAMES  := $(shell uname -s)
# If path to source code not given
ifeq ($(VPATH),$(null))
  VPATH:= .
endif
dirs := $(subst :,$(space),$(VPATH))

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.c $(dir)/*.F $(dir)/*.F90)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

# If executable directory not given
ifeq ($(EXEDIR),$(null))
  EXEDIR  := $(HOME)/bin
  ifeq ($(UNAMES),IRIX64)
    EXEDIR  := $(HOME)/bin/IRIX64
  endif
endif

.SUFFIXES:
.SUFFIXES: .F .o

$(EXEDIR)/$(EXENAME):	$(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------
ifeq ($(UNAMES),IRIX64)
FC      := f90

CPP     := /lib/cpp
FFLAGS  := -64 -DSGI -I$(INC_NETCDF) -mips4 -g -c -C -freeform
LDFLAGS := -64 -g -L$(LIB_NETCDF) -lnetcdf -lfpe
MODULES = *.mod

.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------
ifeq ($(UNAMES),OSF1)
FC      := f90

FFLAGS  := -I$(INC_NETCDF) -g -c -free
LDFLAGS := -L$(LIB_NETCDF) -lnetcdf
MODULES = *.mod

.F.o:
	$(FC) $(FFLAGS) $<

endif

#-----------------------------------------------------------------------
# SUN
#-----------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC      := f90
CPP     := cpp
CFLAGS  := -DSUN -I.. -g
FFLAGS  := -I$(INC_NETCDF) -DSUN -c -O2 -stackvar -xpp=cpp -free
LDFLAGS := -L$(LIB_NETCDF) -L$(HOME)/lib -lnetcdf
MODULES = *.M *.mod

.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# LINUX
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC    := pgf90

FFLAGS := -DLINUX -I$(INC_NETCDF) -c -Mfreeform -Mrecursive -Mdalign -Mfreeform
ifeq ($(DEBUG),TRUE)
  FFLAGS += -g -Mbounds -Ktrap=fp -DDEBUG
else
  FFLAGS += -fast 
endif

LDFLAGS := -g -L$(LIB_NETCDF) -lnetcdf

.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
FC    := xlf90

FFLAGS  := -I$(INC_NETCDF) -c -qarch=auto -qstrict -O3 -qdpc=e -q64
LDFLAGS := -bmaxdata:0x20000000 -q64 -L$(LIB_NETCDF) -lnetcdf


.F.o:
	$(FC) $(FFLAGS) $<

endif

#------------------------------------------------------------------------
# Cray X1
#------------------------------------------------------------------------
ifeq ($(UNAMES),UNICOS/mp)
FC      := ftn

FFLAGS   := -DUNICOSMP -DVECTOR -f free -x omp -Otask0 -em -Ovector3,stream3,scalar3,modinline -dp -s real64 -F -dy -I$(INC_NETCDF) -I. -p.
LDFLAGS  := -L$(LIB_NETCDF) -lnetcdf
MODULES   = *.mod

.F.o:
	$(FC) -c $(FFLAGS) $<

endif


#
# Clean rule
#
clean:
	$(RM) -f $(OBJS) $(MODULES) $(EXEDIR)/$(EXENAME)

precision.o: precision.F
header.o: header.F precision.o
cpr.o: cpr.F header.o stats.o precision.o nldat.o
cprtps.o: cprtps.F header.o nldat.o
header.o: header.F precision.o
initstats.o: initstats.F stats.o precision.o
ismax.o: ismax.F
lenchr.o: lenchr.F
nldat.o: nldat.F
prhddiff.o: prhddiff.F precision.o header.o
printstats.o: printstats.F stats.o precision.o
stats.o: stats.F precision.o
wrap_nf.o: precision.o wrap_nf.F
