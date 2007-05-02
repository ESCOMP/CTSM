#-----------------------------------------------------------------------
# This Makefile is for building clm tools on AIX, Compaq, Linux (with pgf90 or
# lf95 compiler), IRIX or SUN platforms.
#
# These macros can be changed by setting environment variables:
#
# LIB_NETCDF --- Library directory location of netcdf. (defaults to /usr/local/lib)
# INC_NETCDF --- Include directory location of netcdf. (defaults to /usr/local/include)
# USER_FC ------ Allow user to override the default Fortran compiler specified in Makefile.
# USER_CC ------ Allow user to override the default C compiler specified in Makefile (linux only).
# USER_LINKER -- Allow user to override the default linker specified in Makefile.
# USER_CPPDEFS - Additional CPP defines.
# USER_CFLAGS -- Additional C compiler flags that the user wishes to set.
# USER_FFLAGS -- Additional Fortran compiler flags that the user wishes to set.
# USER_LDLAGS -- Additional load flags that the user wishes to set.
#
#------------------------------------------------------------------------

# Set up special characters
null  :=

EXENAME = mksurfdata
RM = rm

# Check for the netcdf library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

# Set user specified Fortran compiler
ifneq ($(strip $(USER_FC)),)
  FC := $(USER_FC)
endif
# Set user specified C compiler
ifneq ($(strip $(USER_CC)),)
  CC := $(USER_CC)
endif


CPPDEF := -DNO_SHR_VMATH -DOFFLINE $(USER_CPPDEFS)

# Determine platform 
UNAMES := $(shell uname -s)

# Load dependency search path.
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(INC_MPI)
cpp_path := $(foreach dir,$(cpp_dirs),-I$(dir)) # format for command line

# Expand any tildes in directory names. Change spaces to colons.
VPATH    := $(foreach dir,$(cpp_dirs),$(wildcard $(dir))) 
VPATH    := $(subst $(space),:,$(VPATH))               

#Primary Target: build the tool
all: $(EXENAME)

# Get list of files and build dependency file for all .o files
#   using perl scripts mkSrcfiles and mkDepends

SOURCES   := $(shell cat Srcfiles)

OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

# Newer makes set the CURDIR variable.
CURDIR := $(shell pwd)

ifeq ($(CLM_ROOT),$(null))
   ROOTDIR := $(CURDIR)/../../
else
   ROOTDIR := $(CLM_ROOT)
endif

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Filepath
	$(ROOTDIR)/bld/mkDepends Filepath Srcfiles > $@


# Architecture-specific flags and rules
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
CPPDEF  += -DIRIX64
FC      = f90
FFLAGS  = -64 -c -r8 -i4 -I$(INC_NETCDF) -O2 $(CPPDEF)
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf
CFLAGS  := -64 -g -O2 $(CPPDEF) -DFORTRANUNDERSCORE
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
CPPDEF += -DAIX -DFORTRAN_SAME
cpre = $(null)-WF,-D$(null)
FPPFLAGS   := $(patsubst -D%,$(cpre)%,$(CPPDEF))
LIB_NETCDF := /usr/local/lib64/r4i4
FC      = xlf90
FFLAGS  = -c -I$(INC_NETCDF) -q64 -qsuffix=f=f90 -O2 -qmaxmem=-1 -qsuffix=f=f90:cpp=F90 \
	   $(FPPFLAGS)  -g
CFLAGS  := -q64 -g $(CPPDEF) -O2
LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)
endif

#------------------------------------------------------------------------
# Darwin
#------------------------------------------------------------------------

ifeq ($(UNAMES),Darwin)
CPPDEF += -DAIX -DDarwin -DFORTRAN_SAME
cpre = $(null)-WF,-D$(null)
FPPFLAGS   := $(patsubst -D%,$(cpre)%,$(CPPDEF))
LIB_NETCDF := /usr/local/lib
CFLAGS  := -I/usr/include -I/usr/include/malloc -g -O2 $(CPPDEF)
FC      = xlf90
FFLAGS  = -c -I$(INC_NETCDF) -qsuffix=f=f90 -O2 -qmaxmem=-1 -qsuffix=f=f90:cpp=F90 \
	  -WF,-DAIX -g -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -C -qfullpath \
	  -qhalt=e $(FPPFLAGS) -I.
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lSystemStubs -lSystemStubs_profile
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
CPPDEF += -DOSF1
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF) $(CPPDEF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
CFLAGS := $(CPPDEF)
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)
endif

#-----------------------------------------------------------------------
# SUN
#-----------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
CPPDEF += -DSUNOS
FC      = f90
FFLAGS  = -c -I$(INC_NETCDF) $(CPPDEF) -stackvar
LDFLAGS = -L$(LIB_NETCDF) -L$(HOME)/lib -lnetcdf
CFLAGS := $(CPPDEF)
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),$(null))
    FC := pgf90
  endif
  CPPDEF += -DLINUX -DFORTRANUNDERSCORE
  CFLAGS := $(CPPDEF)

  ifeq ($(FC),pgf90)
    FFLAGS := -c -I$(INC_NETCDF) $(CPPDEF)
    CC     := pgcc
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g -Ktrap=fp -Mbounds
    else
      FFLAGS += -fast
      CFLAGS += -fast
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
  endif

  ifeq ($(FC),lf95)
    FFLAGS =  -c -I$(INC_NETCDF) $(CPPDEF)
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g --chk a,e,s,u
    else
      FFLAGS += -O
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
  endif
  ifeq ($(FC),pathf90)
    FFLAGS =  -c -I$(INC_NETCDF) $(CPPDEF) -extend_source -ftpp -fno-second-underscore
    ifeq ($(DEBUG),TRUE)
      FFLAGS += -g
    else
      FFLAGS += -O
    endif
    LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
  endif
  FFLAGS += $(cpp_path)
  CFLAGS += $(cpp_path)
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .F90 .c .o

# Append user defined compiler and load flags to Makefile defaults
CFLAGS   += $(USER_CFLAGS)
FFLAGS   += $(USER_FFLAGS)
LDFLAGS  += $(USER_LDFLAGS)

# Set user specified linker
ifneq ($(strip $(USER_LINKER)),)
  LINKER := $(USER_LINKER)
else
  LINKER := $(FC)
endif

.F90.o:
	$(FC) $(FFLAGS) $<

.c.o:
	$(CC) -c $(CFLAGS) $<


$(EXENAME): $(OBJS)
	$(LINKER) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod Depends

include $(CURDIR)/Depends
