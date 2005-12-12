#-----------------------------------------------------------------------
# This Makefile is for building dependency generator "makdep" on a CRAY, 
# SGI or SUN.  The executable will be put in $(HOME)/bin.  If this directory
# does not exist or is not in $PATH, problems will occur when the CCM is built.
#------------------------------------------------------------------------

# Set up special characters
null  :=

# If executable directory not given
ifeq ($(EXEDIR),$(null))
  EXEDIR  := $(HOME)/bin
endif
EXENAME := makdep
# If path to source code not given
ifeq ($(VPATH),$(null))
  VPATH:= .
endif
dirs := $(subst :,$(space),$(VPATH))

# Get list of files and determine objects and dependency files
FIND_FILES = $(wildcard $(dir)/*.c)
FILES      = $(foreach dir, $(dirs),$(FIND_FILES))
SOURCES   := $(sort $(notdir $(FILES)))
OBJS      := $(addsuffix .o, $(basename $(SOURCES)))

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

#
# Make main program
#
$(EXEDIR)/$(EXENAME): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

#
# Debug
#
debug:
	echo "OBJS: $(OBJS)"
	echo "SOURCES: $(SOURCES)"
	echo "FILES: $(FILES)"
	echo "dirs: $(dirs)"
# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
CFLAGS := 
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
CFLAGS := -O -64
LDFLAGS := -64
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
CFLAGS := -O
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
CFLAGS := -O
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
CFLAGS := -O
endif

#------------------------------------------------------------------------
# Default rules
#------------------------------------------------------------------------
.c.o:
	$(CC) -c $(CFLAGS) $<


RM := rm

clean:
	$(RM) -f *.o $(EXEDIR)/$(EXENAME)
