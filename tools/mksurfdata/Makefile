#-----------------------------------------------------------------------
# This Makefile is for building cprnc on AIX, Compaq, Linux (with pgf90 or 
# lf95 compiler), IRIX or SUN platforms. The code produces output nearly 
# identical to that of tools/cprnc. The latter code can be jettisoned when
# test-model goes away. Since test-model relies on cprnc, and demands certain
# quirky features of it (e.g. no .f90 files), the old cprnc code cannot be
# abandoned until test-model is also abandoned.
#
# These macros can be changed by setting environment variables:
#
# LIB_NETCDF --- Library directory location of netcdf. (defaults to /usr/local/lib)
# INC_NETCDF --- Include directory location of netcdf. (defaults to /usr/local/include)
#
#------------------------------------------------------------------------

# Set up special characters
null  :=

EXENAME = mksrfdat
RM = rm

# Check for the netcdf library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Determine platform 
UNAMES := $(shell uname -s)

# Architecture-specific flags and rules
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC      = f90
FFLAGS  = -64 -c -r8 -i4 -I$(INC_NETCDF) -O2 -DIRIX64
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
LIB_NETCDF := /usr/local/lib64/r4i4
FC      = xlf90
FFLAGS  = -c -I$(INC_NETCDF) -q64 -qsuffix=f=f90 -O2 -qmaxmem=-1 -qsuffix=f=f90:cpp=F90 -WF,-DAIX -g
LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF) -DOSF1
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#-----------------------------------------------------------------------
# SUN
#-----------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF) -DSUNOS -stackvar
LDFLAGS = -L$(LIB_NETCDF) -L$(HOME)/lib -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),$(null))
    FC := pgf90
  else
    FC := $(USER_FC)
  endif

  ifeq ($(FC),pgf90)
    ifeq ($(DEBUG),TRUE)
      FFLAGS = -c -I$(INC_NETCDF) -DLINUX -g -Ktrap=fp -Mbounds
    else
      FFLAGS = -c -I$(INC_NETCDF) -DLINUX -fast
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
  endif

  ifeq ($(FC),lf95)
    ifeq ($(DEBUG),TRUE)
      FFLAGS =  -c -I$(INC_NETCDF) -DLINUX -g --chk a,e,s,u
    else
      FFLAGS =  -c -I$(INC_NETCDF) -DLINUX -O
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
  endif
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := fileutils.o shr_kind_mod.o  shr_const_mod.o shr_sys_mod.o ncdio.o domainMod.o creategridMod.o \
        areaMod.o mkvarctl.o mkvarpar.o mkrank.o mkvarsur.o mklaiMod.o mkglacier.o shr_timer_mod.o \
        mklanwat.o mkurban.o mksoitex.o mksoicol.o mkpftMod.o mkfileMod.o mksrfdat.o nanMod.o 

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) $(FFLAGS) $<

$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXENAME)

fileutils.o:
shr_kind_mod.o: 
shr_const_mod.o: shr_kind_mod.o
shr_sys_mod.o:   shr_kind_mod.o
shr_timer_mod.o: shr_kind_mod.o
nanMod.o:        shr_kind_mod.o
domainMod.o:     shr_kind_mod.o nanMod.o
areaMod.o:       shr_kind_mod.o shr_const_mod.o 
creategridMod.o: shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o domainMod.o mkfileMod.o
ncdio.o:         shr_kind_mod.o shr_sys_mod.o mkvarsur.o
mkvarctl.o:      shr_kind_mod.o domainMod.o
mkvarpar.o:      shr_kind_mod.o
mkrank.o:        shr_kind_mod.o 
mkvarsur.o:      shr_kind_mod.o mkvarpar.o domainMod.o
mkfileMod.o:     shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarctl.o ncdio.o 
mklaiMod.o:      shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mkglacier.o:     shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mklanwat.o:      shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mkurban.o:       shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mksoitex.o:      shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mksoicol.o:      shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mkpftMod.o:      shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o
mksrfdat.o:      shr_kind_mod.o shr_sys_mod.o fileutils.o mkvarpar.o mkvarsur.o mkvarctl.o areaMod.o ncdio.o mkfileMod.o mkpftMod.o mklaiMod.o creategridMod.o shr_timer_mod.o

