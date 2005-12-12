#-----------------------------------------------------------------------
# This Makefile is for doing csm_share unit testing
#------------------------------------------------------------------------

VPATH := .:..
.SUFFIXES: .F90 .o
OBJS     := test_shr_sys.o shr_sys_mod.o shr_kind_mod.o
OBJS_ORB := test_shr_orb.o shr_sys_mod.o shr_kind_mod.o shr_orb_mod.o

# Determine platform 
UNAMES := $(shell uname -s)

.F90.o:
	$(FC) -c $(FFLAGS) $<
#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)

FC          := f90
FFLAGS      = -64 -r8 -i4 -cpp -extend_source -DIRIX64 -DHIDE_MPI \
		-g -DEBUG:trap_uninitialized=ON -C
LDFLAGS     = -64

endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC      := f90

FFLAGS  := -xs -stackvar -Qoption f90comp -r8const -e -fast \
        -DSunOS -DHIDE_MPI -g -C -DHIDE_SHR_MSG
LDFLAGS := -fast

endif
#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)

FC = pgf90

CPPFLAGS     = -DHIDE_SHR_MSG -DLinux -DHIDE_MPI
F90FLAGS     = -Mfree
FFLAGS       = $(CPPFLAGS) $(F90FLAGS) -r8 -i4 -Mrecursive -Mdalign \
		-Msecond_underscore -Mextend \
		-g -Ktrap=fp -Mbounds

#.F90.o:
#	$(FC) $(CPPFLAGS) $<
#	$(FC) $(F90FLAGS) $*.f

endif
#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------
ifeq ($(UNAMES),AIX)

FC := xlf90_r
FPPFLAGS := -WF,-DHIDE_SHR_MSG,-DHIDE_MPI,-P,-DAIX
FREEFLAGS  := -qsuffix=f=f90:cpp=F90
FFLAGS   := $(FPPFLAGS) $(FREEFLAGS) -qarch=auto -qrealsize=8  \
	-qstrict -qdpc=e -qsmp=noauto -qspillsize=2500
LDFLAGS    := -qsmp=noauto -bmaxdata:0x80000000

endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)

FC          := f90
FFLAGS_DOTF90 := -DHIDE_SHR_MSG -DOSF1 -free -fpe3
FFLAGS      := -r8 -i4 -omp -automatic -fpe3 -g3 -C  $(FFLAGS_DOTF90)
LDFLAGS     := -omp

endif

test_shr_sys: $(OBJS)
	$(FC) -o test_shr_sys $(LDFLAGS) $(OBJS)

test_shr_sys.o: shr_sys_mod.o shr_sys_mod.o test_shr_sys.F90
	$(FC) -c $(FFLAGS) test_shr_sys.F90

test_shr_orb: $(OBJS_ORB)
	$(FC) -o test_shr_orb $(LDFLAGS) $(OBJS_ORB)

test_shr_sys.o: shr_sys_mod.o test_shr_sys.F90
	$(FC) -c $(FFLAGS) test_shr_sys.F90

test_shr_orb.o: shr_orb_mod.o test_shr_orb.F90
	$(FC) -c $(FFLAGS) test_shr_orb.F90

shr_sys_mod.o: shr_kind_mod.o ../shr_sys_mod.F90
	$(FC) -c $(FFLAGS) ../shr_sys_mod.F90

shr_orb_mod.o: shr_kind_mod.o shr_sys_mod.o ../shr_orb_mod.F90
	$(FC) -c $(FFLAGS) ../shr_orb_mod.F90

clean:
	rm *.mod *.o *.f
