# CIME cmake macro drop-in for the CTSM CI container. Carries the
# container-specific settings CIME needs and ccs_config does not provide:
# where pFUnit is, the gfortran flags ccs_config's own version guard cannot
# apply during a unit-test build, and -- for a serial build -- which ESMF
# flavor, which netCDF, and where the image's mpi-serial and PIO are (see the
# sections below).
#
# WHY THIS FILE EXISTS
# cime/scripts/fortran_unit_testing/run_tests.py locates pFUnit in
# find_pfunit(), which scrapes PFUNIT_PATH out of the CMake macros. It reads
# neither the environment nor --cmake-args:
# ccs_config/machines/cmake_macros/Macros.cmake snapshots the CMake variable
# set BEFORE including the macro files and emits only variables that appeared
# afterwards, so anything passed as -DPFUNIT_PATH=... is filtered right back
# out. PFUNIT_PATH has to be set by a macro *file*.
#
# HOW IT GETS PICKED UP
# CIME's copy_local_macros_to_dir() (cime/CIME/utils.py) copies
# $HOME/.cime/*.cmake into the case's cmake_macros/ directory. Macros.cmake
# includes ${COMPILER}_${MACH}.cmake last, so gnu_container.cmake takes
# precedence over everything ccs_config ships -- and ccs_config has no
# gnu_container.cmake, so nothing is being shadowed.
#
# This keeps the container self-contained: no ccs_config or cime edit needed.
# The eventual upstream fix is a one-line set(PFUNIT_PATH ...) in
# ccs_config/machines/container/container.cmake.
#
# NOTE ON $HOME: GitHub Actions container jobs override HOME (to /github/home),
# so the copy baked into /root/.cime is not found in CI. The workflow step has
# to place this file itself:
#   mkdir -p "$HOME/.cime"
#   cp /opt/ctsm-container/cime-macros/gnu_container.cmake "$HOME/.cime/"
#
# WHY THE SHARED PREFIX AND NOT THE PFUNIT-4.8 SUBDIRECTORY
# run_tests.py turns this value into -DCMAKE_PREFIX_PATH for the test build.
# pFUnit installs four sibling packages under one prefix -- PFUNIT-4.8,
# GFTL-1.11, GFTL_SHARED-1.7, FARGPARSE-1.6 -- and PFUNITConfig.cmake locates
# its own three dependencies ONLY by set()-ing GFTL_ROOT / GFTL_SHARED_ROOT /
# FARGPARSE_ROOT before find_dependency(). Those <pkg>_ROOT variables are
# honored only under policy CMP0074, and CTSM's src/CMakeLists.txt opens with
# cmake_minimum_required(VERSION 2.8), which leaves CMP0074 unset -- so CMake
# ignores them and reports "PFUNIT could not be found because dependency GFTL
# could not be found", after which add_pfunit_ctest is never defined.
#
# Naming the shared prefix sidesteps the policy entirely: find_package locates
# all four as siblings. Do not "tighten" this to the PFUNIT-4.8 subdirectory.
#
# The env-var branch lets an alternate pFUnit build be tested without
# rebuilding the image. An environment variable is not a CMake variable, so
# set()-ing it here is precisely what makes it visible to find_pfunit().
if (DEFINED ENV{PFUNIT_PATH})
  set(PFUNIT_PATH "$ENV{PFUNIT_PATH}")
else()
  set(PFUNIT_PATH "/usr/local/pfunit-4.8.0")
endif()

# ---------------------------------------------------------------------------
# gfortran >= 10 argument-mismatch flags
# ---------------------------------------------------------------------------
# CTSM's mpi-serial code calls mpi_bcast through an implicit interface with
# different actual argument types -- e.g. src/biogeochem/ch4varcon.F90 passes
# a LOGICAL at line 153 and an INTEGER at line 157 -- which gfortran 10+
# rejects outright:
#
#   Error: Type mismatch between actual argument at (1) and actual argument
#          at (2) (INTEGER(4)/LOGICAL(4)).
#
# ccs_config's cmake_macros/gnu.cmake adds -fallow-argument-mismatch for
# exactly this, but guards it on CMAKE_Fortran_COMPILER_VERSION >= 10. That
# guard cannot fire here: CTSM's src/CMakeLists.txt includes
# CIME_initial_setup (which pulls in Macros.cmake, and so this file) at line
# 4, but does not call project() until line 10 -- so when the guard is
# evaluated CMake has not yet probed the compiler and the version variable is
# empty. gnu.cmake's own "Fortran compiler version is" message prints blank in
# the build log, which is the visible symptom.
#
# The container pins gfortran to a single known version in the Dockerfile
# (GCC_VERSION), so the version test has a known answer and the flags can be
# set unconditionally. That stays true across a GCC bump: every gfortran the
# image could pin is >= 10, which is the only thing the guard tests for. This
# file is included last (${COMPILER}_${MACH}.cmake), after gnu.cmake, so
# appending to FFLAGS here reaches CIME_utils.cmake's
#     set(CMAKE_Fortran_FLAGS "${CPPDEFS} ${FFLAGS}")
#
# Not a container-only problem: any machine running CTSM's unit tests with
# gnu hits it. It goes unnoticed upstream because ccs_config defines
# PFUNIT_PATH only for the intel builds on derecho, casper and izumi, so the
# gnu unit-test path is effectively untravelled.
string(APPEND FFLAGS " -fallow-argument-mismatch -fallow-invalid-boz")

# ---------------------------------------------------------------------------
# ESMF flavor for serial (unit-test) builds
# ---------------------------------------------------------------------------
# run_tests.py forces mpilib="mpi-serial", so test executables link with plain
# gfortran and no MPI library -- but CTSM's src/CMakeLists.txt calls
# find_package(ESMF REQUIRED) and link_libraries(esmf) unconditionally. Linking
# an ESMF_COMM=mpich ESMF into a serial executable fails with "libesmf.so:
# undefined reference to symbol 'MPI_Bcast' ... DSO missing from command line",
# because nothing on the link line provides libmpi. The image used to carry two
# such flavors; since item 4b it carries none, and the ESMF_COMM=mpiuni build
# named below is the only one installed.
#
# derecho has a separate ESMF for this case, and its config_machines.xml
# selects it on exactly this condition: mpilib="mpi-serial" loads the plain
# esmf module rather than the -debug variant it loads for MPI builds, even when
# DEBUG="TRUE" (config_machines.xml:90 and :97, both esmf/8.9.1 today). That
# install reports ESMF_COMM=mpiuni and carries no -lmpi. Mirror that here.
#
# FindESMF.cmake (share/cmake) honors an already-defined ESMFMKFILE variable
# before consulting the environment, and this file is included from
# CIME_initial_setup at src/CMakeLists.txt line 4 -- well before the
# find_package(ESMF) at line 24 -- so setting it here wins over the image's
# baked-in ENV ESMFMKFILE. Since item 4b the two name the same file, so this
# is no longer an override; it keeps the macro self-contained and keeps
# working against an image whose default ESMFMKFILE ever moves again.
# EVERYTHING in an mpi-serial executable must speak ONE MPI.
#
# CTSM's mpi-serial build statically links mpi-serial, which defines every
# MPI_* symbol the executable needs (measured: `nm -D --undefined-only
# cesm.exe` reports ZERO undefined MPI_* symbols). A real-MPI ESMF cannot be
# mixed in: libesmf.so carries its own DT_NEEDED on libmpi.so.12, so ESMF's
# MPI calls bind to MPICH while CTSM's bind to mpi-serial. Only mpi-serial is
# then initialized, and the run dies immediately with a 59-byte cesm.log:
#     Attempting to use an MPI routine before initializing MPICH
# Nor can it be fixed by exporting the executable's symbols: ESMF is compiled
# against MPICH headers, so handing it mpi-serial's implementation would mean
# MPICH's MPI_COMM_WORLD constant reaching mpi-serial's functions.
#
# So: the mpiuni ESMF, for unit tests AND cases. It needs PIO built in for the
# CDEPS mesh reads (see the Dockerfile); ESMF disables PIO for mpiuni unless
# ESMF_PIO is set in the environment at build time.
#
# NETCDF_PATH is set for the same one-MPI reason, even though since item 4b it
# names the same prefix ccs_config/machines/container/container.cmake:4 already
# does. The image used to carry a second netCDF at /usr/local built
# --enable-parallel against MPICH (`ldd libnetcdf.so` showed libmpi.so.12), and
# linking that would have put a second, uninitialized MPI in the executable;
# the serial stack lived at /usr/local/serial and this line was what reached
# it. 4b deleted the parallel stack and promoted the serial one to /usr/local,
# so the hazard is gone and this line is now a restatement rather than a
# redirection -- kept because it says, at the point of use, which netCDF an
# mpi-serial executable is allowed to link. derecho splits the same way,
# loading the serial netcdf module (netcdf/4.9.3 today,
# config_machines.xml:87) for mpilib="mpi-serial" rather than netcdf-mpi.
# cime/CIME/Tools/Makefile turns NETCDF_PATH into INC_NETCDF/LIB_NETCDF and
# drops PNETCDF_PATH outright for mpi-serial (lines 214-217), which is what
# makes the container machine's now-unbacked PNETCDF_PATH=/usr/local harmless.
#
# That prefix is static (see "WHY /usr/local, AND WHY STATIC" in the
# Dockerfile), so there is no .so for the loader to pick at all; the image no
# longer ships an /etc/ld.so.conf.d/ctsm.conf either, because nothing shared
# installs there.
if (MPILIB STREQUAL "mpi-serial")
  set(ESMFMKFILE "/usr/local/esmf-8.9.1-mpiuni/lib/esmf.mk")
  set(NETCDF_PATH "/usr/local")

  # ---------------------------------------------------------------------------
  # The image's mpi-serial and PIO
  # ---------------------------------------------------------------------------
  # Without these a container case build compiles CTSM's own
  # libraries/mpi-serial and libraries/parallelio submodules, so the two sides
  # link different mpi-serials -- 2.5.4 here against derecho's 2.5.3 -- while
  # every version check reports a match. This mirrors derecho, where
  # ccs_config/machines/derecho/derecho.cmake sets MPI_SERIAL_PATH from
  # $ENV{NCAR_ROOT_MPI_SERIAL} (line 4) and PIO_LIBDIR/PIO_INCDIR from
  # $ENV{PIO} (lines 8-10), both exported by the modules it loads. The
  # container has no modules, so they are hardcoded, as PFUNIT_PATH and
  # ESMFMKFILE already are.
  #
  # THESE PATHS REQUIRE THE ITEM-4b IMAGE. They name the serial stack at
  # /usr/local with a Fortran-enabled PIO in it, which is what the 4b Dockerfile
  # builds. Any earlier image -- the published :20261007 and before -- put the
  # serial stack under /usr/local/serial and built PIO with
  # -DPIO_ENABLE_FORTRAN=OFF, since PIO was only ever meant to satisfy ESMF,
  # which links -lpioc and not -lpiof. A container case build against one of
  # those fails, deliberately and loudly, rather than silently reverting to the
  # submodules.
  #
  # THE TWO SETTINGS ARE ONE SWITCH, which was measured rather than reasoned.
  # With only MPI_SERIAL_PATH set, Tools/Makefile:912-918 points $(MPISERIAL)
  # at the image's archive; that file already exists, so make skips the
  # submodule build and leaves $(INSTALL_SHAREDPATH) holding only libgptl.a --
  # while line 734 still hands PIO's from-source build
  # -D MPISERIAL_PATH=$(INSTALL_SHAREDPATH), ignoring MPI_SERIAL_PATH
  # entirely. PIO's configure then dies with "Could NOT find MPISERIAL_C". No
  # case-level override exists: CMAKE_OPTS is appended again at 738, 915 and
  # 922, and none of those sets MPISERIAL_PATH.
  #
  # WHAT THIS FILE CAN AND CANNOT DECIDE. It can only set make variables, and
  # that is enough for the CTSM component build: Tools/Makefile:463-466 turns
  # PIO_INCDIR and PIO_LIBDIR into -I and -L. It cannot decide whether PIO is
  # built from source at all. That is cime/CIME/build_scripts/buildlib.pio:82-91,
  # which reads PIO_VERSION_MAJOR from the ENVIRONMENT -- os.environ, which no
  # macro file reaches -- and returns before running make. Setting PIO_LIBDIR
  # here does not substitute: buildlib.pio:115 passes its own PIO_LIBDIR on the
  # gmake command line, which beats anything Macros.make carries.
  #
  # SO THE IMAGE MUST EXPORT PIO_VERSION_MAJOR=2 AND PIO_TYPENAME_VALID_VALUES.
  # buildlib.pio:76-81 states the contract. Miss PIO_TYPENAME_VALID_VALUES and
  # line 90 hands None to _set_pio_valid_values, which does
  # valid_values += ",nothing" at line 224 and raises TypeError. Miss
  # PIO_VERSION_MAJOR and PIO is built from source while MPI_SERIAL_PATH has
  # already emptied the directory that build needs -- the "Could NOT find
  # MPISERIAL_C" failure above.
  #
  # MPI_SERIAL_PATH is /usr/local/mpi-serial, NOT /usr/local. The Dockerfile
  # copies the archive into both -- the mpiuni ESMF's link line wants it beside
  # libpioc -- but only the dedicated prefix carries mpi.h and mpif.h, kept out
  # of /usr/local/include because that directory is PIO_INCDIR and therefore on
  # the -I path of every case compile. The Dockerfile asserts this path at
  # image-build time, beside its PFUNIT_PATH and ESMFMKFILE assertions.
  #
  # ACCEPTANCE TEST, met on 2026-10-09 by a throwaway probe image (since folded
  # into the Dockerfile and deleted; see NEXT_STEPS.md item 4b): a
  # freshly built cesm.exe contains neither "/ctsm/libraries/mpi-serial" nor
  # "/ctsm/libraries/parallelio" among its compiled-source strings, and does
  # contain the image's own prefixes (strings -a cesm.exe). A passing version
  # check without that is the old vacuous behaviour unchanged.
  set(MPI_SERIAL_PATH "/usr/local/mpi-serial")
  set(PIO_LIBDIR "/usr/local/lib")
  set(PIO_INCDIR "/usr/local/include")
endif()
