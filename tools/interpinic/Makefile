#-----------------------------------------------------------------------
# This Makefile is for building clm tools on AIX, Linux (with pgf90 or 
# lf95 compiler), Darwin or IRIX platforms.
#
# These macros can be changed by setting environment variables:
#
# LIB_NETCDF --- Library directory location of netcdf. (defaults to /usr/local/lib)
# INC_NETCDF --- Include directory location of netcdf. (defaults to /usr/local/include)
# MOD_NETCDF --- Module directory location of netcdf.  (defaults to $LIB_NETCDF)
# USER_FC ------ Allow user to override the default Fortran compiler specified in Makefile.
# USER_CC ------ Allow user to override the default C compiler specified in Makefile (linux only).
# USER_LINKER -- Allow user to override the default linker specified in Makefile.
# USER_CPPDEFS - Additional CPP defines.
# USER_CFLAGS -- Additional C compiler flags that the user wishes to set.
# USER_FFLAGS -- Additional Fortran compiler flags that the user wishes to set.
# USER_LDLAGS -- Additional load flags that the user wishes to set.
# SMP ---------- Shared memory Multi-processing (TRUE or FALSE) [default is FALSE]
# OPT ---------- Use optimized options.
#
#------------------------------------------------------------------------

# Set up special characters
null  :=

EXENAME = interpinic
RM = rm

# Check for the netcdf library and include directories 
ifeq ($(LIB_NETCDF),$(null))
LIB_NETCDF := /usr/local/lib
endif

ifeq ($(INC_NETCDF),$(null))
INC_NETCDF := /usr/local/include
endif

ifeq ($(MOD_NETCDF),$(null))
MOD_NETCDF := $(LIB_NETCDF)
endif

# Set user specified Fortran compiler
ifneq ($(strip $(USER_FC)),)
  FC := $(USER_FC)
endif
# Set user specified C compiler
ifneq ($(strip $(USER_CC)),)
  CC := $(USER_CC)
endif

# Set if Shared memory multi-processing will be used
ifneq ($(strip $(SMP)),)
  SMP := FALSE
endif


CPPDEF := -DNO_SHR_VMATH -DOFFLINE $(USER_CPPDEFS)

# Determine platform 
UNAMES := $(shell uname -s)

# Load dependency search path.
dirs := . $(shell cat Filepath)

# Set cpp search path, include netcdf
cpp_dirs := $(dirs) $(INC_NETCDF) $(INC_MPI) $(MOD_NETCDF)
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
   ROOTDIR := $(shell ls -1d $(CLM_ROOT)/models/lnd/clm*)
endif

$(CURDIR)/Depends: $(CURDIR)/Srcfiles $(CURDIR)/Filepath
	$(ROOTDIR)/../../../scripts/ccsm_utils/Build/mkDepends Filepath Srcfiles > $@


# Architecture-specific flags and rules
#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
CPPDEF += -DAIX -DFORTRAN_SAME
cpre = $(null)-WF,-D$(null)
FPPFLAGS   := $(patsubst -D%,$(cpre)%,$(CPPDEF))
LIB_NETCDF := /usr/local/lib64/r4i4
FC      = xlf90_r
FFLAGS  = -c -I$(INC_NETCDF) -q64 -qsuffix=f=f90 -O2 -qsuffix=f=f90:cpp=F90 \
	   $(FPPFLAGS)  -g -qfullpath -qarch=auto -qspillsize=2500  \
           -qflttrap=ov:zero:inv:en
ifneq ($(OPT),TRUE)
   FFLAGS += -qmaxmem=-1 -qinitauto=FF911299 -C
endif
CFLAGS  := -q64 -g $(CPPDEF) -O2 -C
LDFLAGS = -L$(LIB_NETCDF) -q64 -lnetcdf
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)

ifeq ($(SMP),TRUE)
  FFLAGS  += -qsmp=omp
  LDFLAGS += -qsmp=omp
endif

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
FFLAGS  = -c -I$(INC_NETCDF) -qsuffix=f=f90 -O2 -qsuffix=f=f90:cpp=F90 \
	  -WF,-DAIX -g -qfullpath -qhalt=e $(FPPFLAGS) -I.
ifneq ($(OPT),TRUE)
   FFLAGS  += -qmaxmem=-1 -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -C
else
   FFLAGS  += -qmaxmem=5000
endif
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf -lSystemStubs -lSystemStubs_profile
FFLAGS += $(cpp_path)
CFLAGS += $(cpp_path)
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
  ifeq ($(USER_FC),$(null))
    FC    := pgf90
    PGF90 := TRUE
  endif
  ifeq ($(USER_FC),ftn)
    PGF90 := TRUE
  endif
  CPPDEF += -DLINUX -DFORTRANUNDERSCORE
  CFLAGS := $(CPPDEF)
  LDFLAGS = -L$(LIB_NETCDF) -lnetcdf

  ifeq ($(PGF90),TRUE)
    CC     := pgcc
    ifneq ($(OPT),TRUE)
      FFLAGS += -g -Ktrap=fp -Mbounds -Kieee
    else
      FFLAGS += -fast -Kieee
      CFLAGS += -fast
    endif

    ifeq ($(SMP),TRUE)
      FFLAGS   += -mp
      LDFLAGS  += -mp
    endif

  endif

  ifeq ($(FC),lf95)
    ifneq ($(OPT),TRUE)
      FFLAGS += -g --chk a,e,s,u
    else
      FFLAGS += -O
    endif
    # Threading only works by putting thread memory on the heap rather than the stack
    # (--threadheap).
    # As of lf95 version 6.2 the thread stacksize limits are (still) too small to run
    # even small
    # resolution problems (FV at 10x15 res fails).
    ifeq ($(SMP),TRUE)
      FFLAGS  += --openmp --threadheap 4096
      LDFLAGS += --openmp --threadheap 4096
    endif
  endif
  ifeq ($(FC),pathf90)
    FFLAGS = -extend_source -ftpp -fno-second-underscore
    ifneq ($(OPT),TRUE)
      FFLAGS += -g
    else
      FFLAGS += -O
    endif
    ifeq ($(SMP),TRUE)
      FFLAGS   += -mp
      LDFLAGS  += -mp
    endif
  endif
  FFLAGS += -c -I$(INC_NETCDF) $(CPPDEF) $(cpp_path)
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
