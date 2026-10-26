# Container image for building CTSM in GitHub Actions CI.
#
# Built from scratch on AlmaLinux 9 so every component can match what
# CIME loads for a gnu + mpi-serial build on derecho (per
# ccs_config/machines/derecho/config_machines.xml) exactly where an open
# equivalent exists; see README.md for the version-match table and the two
# accepted deviations (cray-libsci -> reference LAPACK/BLAS, and the distro
# runtime).
#
# ONE SOFTWARE STACK, AND IT IS SERIAL. There is no MPICH here, no parallel
# HDF5/netCDF, no PnetCDF and no ESMF_COMM=mpich flavor: a gnu + mpi-serial run
# on derecho loads none of them, and a second stack only ever created ways for
# a build to link the wrong half. Everything installs under /usr/local, which
# is where CIME's "container" machine hard-codes NETCDF_PATH
# (ccs_config/machines/container/container.cmake:4, and the matching <env> in
# that directory's config_machines.xml). The same machine hard-codes
# PNETCDF_PATH=/usr/local, which is now a dead letter rather than a lie:
# cime/CIME/Tools/Makefile:214-217 undefines PNETCDF_PATH outright for
# MPILIB=mpi-serial, so nothing goes looking for the PnetCDF that is no longer
# installed. The toolchain lives under /opt. ESMF is located via ESMFMKFILE
# instead, so its one build flavor (ESMF_COMM=mpiuni) lives in its own prefix.
#
# The image serves two CIME workflows, both of them mpi-serial: create_test
# builds, and CTSM's Fortran unit tests
# (cime/scripts/fortran_unit_testing/run_tests.py, which forces
# mpilib=mpi-serial). A create_test case has to ask for mpi-serial explicitly
# -- "_Mmpi-serial" in the test name, or --mpilib mpi-serial -- because the
# container machine's <MPILIBS> still says mpich, so a test that names no
# modifier takes that default and finds no MPI in this image.
#
# A case build here links the IMAGE's mpi-serial and PIO rather than compiling
# CTSM's libraries/mpi-serial and libraries/parallelio submodules, which is
# what derecho gets from its mpi-serial and parallelio-serial modules. That
# takes two things together: cime-macros/gnu_container.cmake (MPI_SERIAL_PATH,
# PIO_LIBDIR, PIO_INCDIR) and the PIO environment variables in the ENV block
# near the bottom of this file. Neither works without the other; see those
# sections.
#
# All environment needed at runtime is baked in with ENV: GitHub
# Actions runs steps in non-login shells, so profile.d scripts never
# fire in CI.

FROM docker.io/almalinux:9

ARG GCC_VERSION=14.3.0
ARG HDF5_VERSION=1.14.6
ARG NETCDF_C_VERSION=4.9.3
ARG NETCDF_FORTRAN_VERSION=4.6.2
ARG ESMF_VERSION=8.9.1
ARG PFUNIT_VERSION=4.8.0
# These two are load-bearing for a CTSM case build in the container, not just
# for the image's internals. cime-macros/gnu_container.cmake sets
# MPI_SERIAL_PATH, PIO_LIBDIR and PIO_INCDIR, and the ENV block near the bottom
# of this file exports PIO/PIO_VERSION_MAJOR/PIO_TYPENAME_VALID_VALUES, so CIME
# links these installs instead of compiling CTSM's own libraries/mpi-serial and
# libraries/parallelio submodules -- mirroring derecho, where
# ccs_config/machines/derecho/derecho.cmake:4,8-10 reads the same three
# settings out of what the mpi-serial and parallelio-serial modules export. The
# serial PIO is ALSO the mpiuni ESMF's external PIO (see "ESMF, the only
# flavor" below), so one install serves both; that is why it is built first.
#
# They track derecho like every other version ARG here except one: GIT_VERSION
# (git is not part of derecho's stack at all; see the git section below). The
# image mirrors what you get from CTSM's current ccs_config plus derecho's
# default software stack, so check-derecho-versions.py compares both "direct"
# against the mpilib="mpi-serial" <modules> blocks of derecho's
# config_machines.xml (mpi-serial and parallelio-serial).
ARG MPI_SERIAL_VERSION=2.5.3
ARG PIO_VERSION=2.6.8

# Parallelism for the source builds below. Bounded (not $(nproc)) because
# podman inside a PBS job can report the whole node's CPU count rather than
# the cores actually allocated, which oversubscribes the build. Default 16
# matches the recommended ncpus=16 build allocation; override with
# --build-arg MAKE_JOBS=N.
ARG MAKE_JOBS=16

SHELL ["/bin/bash", "-c"]

# System packages. lapack-devel/blas-devel (stand-ins for derecho's
# proprietary cray-libsci; the CIME "container" machine links
# -llapack -lblas) live in the CRB repo. The perl modules are needed by
# CTSM's build-namelist. The system gcc/g++ only bootstraps the
# GCC ${GCC_VERSION} build below.
#
# git is NOT dnf-installed: it pulls in openssh, whose rpm chowns a setuid
# helper (ssh-keysign) to a non-root id, which fails under NCAR's rootless
# single-UID podman ("cpio: chown failed"). git IS needed at runtime --
# CIME's cprnc/PIO builds clone genf90 via git -- so it is built from
# source in a later layer instead (see the git section near the bottom).
RUN dnf -y install dnf-plugins-core \
    && dnf config-manager --set-enabled crb \
    && dnf -y install \
        gcc gcc-c++ \
        make cmake m4 patch file diffutils hostname which findutils \
        wget bzip2 xz \
        python3 \
        perl-core perl-XML-LibXML \
        zlib-devel libcurl-devel \
        lapack-devel blas-devel \
    && dnf clean all

# Extract source tarballs as the build user, not the uid/gid recorded in
# the archive. Upstream tarballs (GCC's, etc.) are owned by the packager's
# uid/gid; under rootless single-UID podman tar cannot chown to an unmapped
# id ("Cannot change ownership ... Invalid argument"). TAR_OPTIONS is
# honored by every tar invocation, including the nested ones inside GCC's
# contrib/download_prerequisites.
ENV TAR_OPTIONS=--no-same-owner

# Compute-node network can be flaky; have wget retry connection refusals
# (it already retries most other transient errors by default).
RUN printf '%s\n' 'retry_connrefused = on' 'tries = 10' 'waitretry = 15' >> /etc/wgetrc

# ---- GCC (matches derecho's gcc/14.3.0 module, config_machines.xml:59) ----
# The tarball is tried against several mirrors in turn because ftp.gnu.org (and
# ftpmirror.gnu.org, which resolves to the same host) is unreachable from NCAR
# compute nodes -- it times out rather than refusing, so wget burns its full
# retry budget and the step exits 4 before anything compiles. DNS resolves
# fine; this is routing. The per-attempt --tries/--timeout override the
# deliberately patient global wgetrc below so a dead mirror is abandoned
# quickly instead of eating the job's walltime.
#
# gmp/mpfr/mpc come from AlmaLinux (dnf) rather than GCC's
# contrib/download_prerequisites, which fetches them from gcc.gnu.org --
# unreliable from the compute node. The distro ships the versions GCC 12.2
# expects (gmp 6.2, mpfr 4.1.0, libmpc 1.2.1), so the compiler is
# equivalent and skipping their in-tree build is faster. configure finds
# them in the default /usr paths (no --with-* needed). Installed in this
# layer, not the base dnf layer, to keep that layer's cache valid.
# --disable-bootstrap: single-stage build; same compiler, just not
# self-compiled three times.
RUN dnf -y install gmp-devel mpfr-devel libmpc-devel && dnf clean all \
    && cd /tmp \
    && for base in "https://gcc.gnu.org/pub/gcc/releases" \
                   "https://mirrors.kernel.org/gnu/gcc" \
                   "https://ftp.gnu.org/gnu/gcc"; do \
           echo "fetching gcc from ${base}"; \
           wget -nv --tries=3 --timeout=30 -O "gcc-${GCC_VERSION}.tar.xz" \
               "${base}/gcc-${GCC_VERSION}/gcc-${GCC_VERSION}.tar.xz" && break; \
       done \
    && tar xf "gcc-${GCC_VERSION}.tar.xz" \
    && cd "gcc-${GCC_VERSION}" \
    && mkdir build && cd build \
    && ../configure --prefix="/opt/gcc/${GCC_VERSION}" \
        --enable-languages=c,c++,fortran \
        --disable-multilib --disable-bootstrap --disable-nls \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && printf '%s\n' "/opt/gcc/${GCC_VERSION}/lib64" > /etc/ld.so.conf.d/gcc.conf \
    && ldconfig \
    && rm -rf /tmp/gcc-*

ENV PATH=/usr/local/bin:/opt/gcc/${GCC_VERSION}/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin \
    LD_LIBRARY_PATH=/usr/local/lib:/opt/gcc/${GCC_VERSION}/lib64

# WHY /usr/local, AND WHY STATIC
# The serial stack is the only stack, so it installs at the prefix CIME's
# "container" machine hard-codes NETCDF_PATH to
# (ccs_config/machines/container/container.cmake:4). SERIAL_PREFIX stays a
# variable because every step below names it, and keeping it keeps this layout
# legible against the two-stack one it replaces.
#
# Static (--disable-shared) is no longer about shadowing. That WAS the reason:
# a parallel libnetcdf.so also lived in /usr/local/lib, so a shared serial
# build would have linked correctly and then silently loaded the parallel one.
# There is no parallel twin any more. Static stays for two reasons that remain
# true:
#   - PIO has to be static. cime/CIME/Tools/Makefile:667 names
#     $(PIO_LIBDIR)/libpiof.a and libpioc.a BY PATH as prerequisites of both
#     $(CSMSHARELIB) (line 931) and the executable (line 960), and the only
#     rule for them (line 929) carries no recipe -- shared PIO libraries stop a
#     case build outright with "No rule to make target".
#   - PIO, the mpiuni ESMF and every case executable link this netCDF/HDF5
#     chain, so it is built the same way. Static also puts each library's
#     version banner inside cesm.exe, which is where DRIFT_CHECK_DESIGN.md's
#     instrument map reads the container's netCDF-C and HDF5 versions from
#     (`strings`); a shared stack would move them into the .so instead.
# -fPIC because ESMF links these into position-independent objects.
#
# NOTHING SHARED LANDS IN /usr/local/lib, which is why there is no
# /etc/ld.so.conf.d/ctsm.conf any more. That file was written by the MPICH
# block, listed /usr/local/lib plus the mpich lib dir, and went with it.
# Traced over what is left: HDF5, netCDF-C and netCDF-Fortran are
# --disable-shared --enable-static; PIO is BUILD_SHARED_LIBS=OFF; mpi-serial is
# a hand-copied .a; git installs to bin/, libexec/ and share/; pFUnit and ESMF
# have their own prefixes, which ctsm.conf never listed either. ldconfig has
# nothing to cache here. ENV LD_LIBRARY_PATH above still lists /usr/local/lib,
# so anything shared that did appear would still be found at run time.
ENV SERIAL_PREFIX=/usr/local

# ---- HDF5 (serial) ----
# Note the compilers: named by absolute path into /opt. They used to be an
# override of the ENV CC=mpicc/FC=mpifort the MPICH block baked in; those are
# gone, but the explicit paths stay, because AlmaLinux's bootstrap gcc (11.5.0,
# the one that built GCC ${GCC_VERSION}) is still installed in /usr/bin and
# nothing in this stack should be able to fall back to it.
RUN cd /tmp \
    && wget -nv "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_${HDF5_VERSION}.tar.gz" \
    && tar xzf "hdf5_${HDF5_VERSION}.tar.gz" \
    && cd "hdf5-hdf5_${HDF5_VERSION}" \
    && CC="/opt/gcc/${GCC_VERSION}/bin/gcc" \
       CXX="/opt/gcc/${GCC_VERSION}/bin/g++" \
       FC="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
       CFLAGS=-fPIC CXXFLAGS=-fPIC FCFLAGS=-fPIC \
       ./configure --prefix="${SERIAL_PREFIX}" --enable-hl --with-zlib \
        --disable-shared --enable-static \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/hdf5-*

# ---- netCDF-C (serial) ----
# libxml2-devel is installed here rather than in the base dnf layer to keep the
# (expensive) GCC layer cache valid; fold it into the base layer whenever a
# clean rebuild happens anyway. It used to be installed by the MPI-linked
# netCDF-C build that ran before this one; that build is gone, and this is now
# the first step that needs it. netCDF-C links libxml2 -- it is in the
# -l list `nc-config --libs --static` returns, which the netCDF-Fortran step
# below depends on.
RUN dnf -y install libxml2-devel && dnf clean all \
    && cd /tmp \
    && wget -nv "https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF_C_VERSION}/netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
    && tar xzf "netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
    && cd "netcdf-c-${NETCDF_C_VERSION}" \
    && CC="/opt/gcc/${GCC_VERSION}/bin/gcc" \
       CFLAGS=-fPIC \
       CPPFLAGS="-I${SERIAL_PREFIX}/include" \
       LDFLAGS="-L${SERIAL_PREFIX}/lib" \
       ./configure --prefix="${SERIAL_PREFIX}" --enable-netcdf-4 \
        --disable-shared --enable-static --disable-dap --disable-byterange \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/netcdf-c-*

# ---- netCDF-Fortran (serial) ----
# --libs --static, not --libs. netCDF-C 4.9.3 dropped $LIBS from the NC_LIBS it
# substitutes into nc-config (configure.ac: NC_LIBS="$LDFLAGS $NC_LIBS" where
# 4.9.2 had "$LDFLAGS $NC_LIBS $LIBS"), so plain --libs now yields a bare
# -lnetcdf. That is fine against a shared netCDF, where DT_NEEDED pulls the
# rest in, but this prefix is static-only: without the transitive
# -lhdf5_hl -lhdf5 -lm -lz -lxml2 -lcurl the configure link test for nc_open
# fails and the step dies with "Could not link to netcdf C library".
RUN cd /tmp \
    && wget -nv "https://downloads.unidata.ucar.edu/netcdf-fortran/${NETCDF_FORTRAN_VERSION}/netcdf-fortran-${NETCDF_FORTRAN_VERSION}.tar.gz" \
    && tar xzf "netcdf-fortran-${NETCDF_FORTRAN_VERSION}.tar.gz" \
    && cd "netcdf-fortran-${NETCDF_FORTRAN_VERSION}" \
    && CC="/opt/gcc/${GCC_VERSION}/bin/gcc" \
       FC="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
       F77="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
       CFLAGS=-fPIC FCFLAGS=-fPIC FFLAGS=-fPIC \
       CPPFLAGS="-I${SERIAL_PREFIX}/include" \
       LDFLAGS="-L${SERIAL_PREFIX}/lib" \
       LIBS="$(${SERIAL_PREFIX}/bin/nc-config --libs --static)" \
       ./configure --prefix="${SERIAL_PREFIX}" --disable-shared --enable-static \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/netcdf-fortran-*

# ---- mpi-serial: PIO builds against it, and so does every case ----
# cime-macros/gnu_container.cmake points MPI_SERIAL_PATH at ${MPISERIAL_PREFIX}
# below, so cime/CIME/Tools/Makefile:912-918 links THIS archive and skips
# building CTSM's libraries/mpi-serial submodule. The copy into
# ${SERIAL_PREFIX}/lib is what the mpiuni ESMF's link line asks for
# (ESMF_PIO_LIBS="-lpioc -lmpi-serial" against ESMF_PIO_LIBPATH, below).
#
# The HEADERS go only into the dedicated prefix, never ${SERIAL_PREFIX}: that
# is now /usr/local, whose include/ is PIO_INCDIR and therefore on the -I path
# of every case compile (Tools/Makefile:465). Keeping mpi.h/mpif.h out of it
# means exactly one mpi.h is reachable, the one MPI_SERIAL_PATH names.
#
# The files are copied by hand rather than via `make install`: mpi-serial's
# install rule expands $(INSTALL) and $(MKINSTALLDIRS), and neither is defined
# anywhere in the release, so the rule runs "libmpi-serial.a -m 644 ..." and
# dies with "libmpi-serial.a: No such file or directory". CIME's own
# buildlib.mpi-serial copies the artifacts directly for the same reason.
ENV MPISERIAL_PREFIX=/usr/local/mpi-serial
RUN cd /tmp \
    && wget -nv "https://github.com/ESMCI/mpi-serial/archive/refs/tags/MPIserial_${MPI_SERIAL_VERSION}.tar.gz" \
    && tar xzf "MPIserial_${MPI_SERIAL_VERSION}.tar.gz" \
    && cd "mpi-serial-MPIserial_${MPI_SERIAL_VERSION}" \
    && CC="/opt/gcc/${GCC_VERSION}/bin/gcc" \
       FC="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
       CFLAGS=-fPIC FCFLAGS=-fPIC \
       ./configure --prefix="${MPISERIAL_PREFIX}" \
    && make -j"${MAKE_JOBS}" \
    && mkdir -p "${MPISERIAL_PREFIX}/lib" "${MPISERIAL_PREFIX}/include" \
    && cp libmpi-serial.a "${MPISERIAL_PREFIX}/lib/" \
    && cp mpi.h mpif.h "${MPISERIAL_PREFIX}/include/" \
    && cp libmpi-serial.a "${SERIAL_PREFIX}/lib/" \
    && (cp mpi.mod "${MPISERIAL_PREFIX}/include/" 2>/dev/null || true) \
    && rm -rf /tmp/mpi-serial-* /tmp/MPIserial_*

# ---- git (built from source; NO openssh dependency) ----
# CIME build steps clone helper repos at build time via CMake
# ExternalProject_Add: cprnc fetches genf90 (github.com/PARALLELIO/genf90),
# and PIO/CDEPS/CMEPS generate code with genf90 too -- so the image needs a
# real git at runtime. git cannot be dnf-installed here: it pulls openssh,
# whose rpm chowns a setuid helper (ssh-keysign) to a non-root id, which
# fails under rootless single-UID podman ("cpio: chown failed"). Building
# from source avoids openssh entirely -- https clones go through libcurl --
# and behaves like a stock git. NO_OPENSSL/NO_EXPAT/NO_GETTEXT/NO_TCLTK drop
# optional deps we don't ship.
#
# IT SITS HERE, BEFORE PIO, because the PIO layer below needs it AT BUILD TIME:
# with -DPIO_ENABLE_FORTRAN=ON, ParallelIO's CMakeLists.txt:137-158 shells out
# to `git clone CMake_Fortran_utils` unless USER_CMAKE_MODULE_PATH is defined,
# and that clone is what defines it. The whole block is inside
# `if (PIO_ENABLE_FORTRAN)`, which is why the old FORTRAN=OFF PIO build could
# sit ahead of git and still configure. Still placed after the toolchain, so
# editing this layer does not rebuild GCC or the netCDF stack -- it does now
# rebuild PIO, ESMF and pFUnit, which is the price of the move and is cheap
# next to those.
ARG GIT_VERSION=2.45.2
RUN cd /tmp \
    && wget -nv "https://mirrors.edge.kernel.org/pub/software/scm/git/git-${GIT_VERSION}.tar.xz" \
    && tar xf "git-${GIT_VERSION}.tar.xz" \
    && cd "git-${GIT_VERSION}" \
    && make -j"${MAKE_JOBS}" prefix=/usr/local CC=gcc \
        NO_OPENSSL=1 NO_EXPAT=1 NO_GETTEXT=1 NO_TCLTK=1 all \
    && make prefix=/usr/local \
        NO_OPENSSL=1 NO_EXPAT=1 NO_GETTEXT=1 NO_TCLTK=1 install \
    && rm -rf /tmp/git-*

# ---- ParallelIO, with the Fortran interface, against mpi-serial + netCDF ----
# ONE PIO SERVES TWO CONSUMERS. The mpiuni ESMF below links it as its external
# PIO (-lpioc), and a CTSM case build links it through the PIO_LIBDIR/PIO_INCDIR
# that cime-macros/gnu_container.cmake sets -- which needs -lpiof and pio.mod,
# hence PIO_ENABLE_FORTRAN=ON. (It was OFF before item 4b, when PIO existed only
# for ESMF and a case compiled CTSM's libraries/parallelio submodule instead.)
#
# No pnetcdf and no MPI-IO: both would reintroduce a real-MPI dependency, which
# is the whole thing this image avoids.
#
# FLAGS ARE DERECHO'S, NOT GUESSES. parallelio-serial/2.6.8.lua carries the full
# configure line in its whatis; these are the ones that matter:
#   -DPIO_ENABLE_FORTRAN=ON   -DPIO_ENABLE_TIMING=OFF
#   -DPIO_USE_MPISERIAL=ON    -DPIO_ENABLE_LOGGING=OFF
#   -DWITH_PNETCDF=OFF        -DPIO_ENABLE_TESTS=OFF
# Timing stays OFF on both sides: CIME builds its own GPTL and links it
# separately (cime/CIME/Tools/Makefile:674, ULIBS += ... -lgptl).
#
# TWO DELIBERATE DEVIATIONS from that modulefile:
#   BUILD_SHARED_LIBS -- derecho says ON, this says OFF, for the reason given at
#   "WHY /usr/local, AND WHY STATIC" above: Tools/Makefile names libpiof.a and
#   libpioc.a by path as recipe-less prerequisites, so a shared PIO stops a case
#   build with "No rule to make target".
#   PIO_ENABLE_NETCDF_INTEGRATION -- derecho says ON, this leaves it at PIO's
#   own default of OFF (CMakeLists.txt:91). The flag only adds source files to
#   libpioc/libpiof (src/clib/CMakeLists.txt:28, src/flib/CMakeLists.txt:34)
#   that nothing in CTSM or CDEPS calls, and with a static PIO the linker never
#   pulls them, so cesm.exe -- which is what the item-5 drift check compares --
#   is unaffected. Turning it on would instead need netCDF-C dispatch support
#   this image has not been checked for, and a build that quietly disabled it
#   would leave the docs claiming a parity that was not there.
#
# GENF90_PATH and USER_CMAKE_MODULE_PATH are passed so that nothing is fetched
# during cmake: src/flib/CMakeLists.txt:132-147 ExternalProject-clones genf90
# when GENF90_PATH is undefined, and CMakeLists.txt:137-158 shells out to
# `git clone CMake_Fortran_utils` when USER_CMAKE_MODULE_PATH is. genf90.pl is a
# tracked file in the release tarball; CMake_Fortran_utils is not vendored
# anywhere, so it is cloned here, once, visibly -- with the git built in the
# layer above.
#
# The clone retries. Everything else that reaches the network here goes through
# wget, which /etc/wgetrc configures to be patient (see the wgetrc step near the
# top: compute-node networking is flaky). git honours none of that, and this
# layer sits an hour of GCC into the build, so a single transient failure would
# be expensive. Five attempts with a short backoff, then give up loudly.
RUN cd /tmp \
    && for i in 1 2 3 4 5; do \
           git clone --depth 1 https://github.com/CESM-Development/CMake_Fortran_utils && break; \
           echo "CMake_Fortran_utils clone failed (attempt $i); retrying"; \
           rm -rf /tmp/CMake_Fortran_utils; \
           sleep 15; \
       done \
    && test -f /tmp/CMake_Fortran_utils/mpiexec.cmake \
    && wget -nv "https://github.com/NCAR/ParallelIO/archive/refs/tags/pio${PIO_VERSION//./_}.tar.gz" \
    && tar xzf "pio${PIO_VERSION//./_}.tar.gz" \
    && cd "ParallelIO-pio${PIO_VERSION//./_}" \
    && src="$PWD" \
    && mkdir build && cd build \
    && cmake \
        -DCMAKE_INSTALL_PREFIX="${SERIAL_PREFIX}" \
        -DCMAKE_C_COMPILER="/opt/gcc/${GCC_VERSION}/bin/gcc" \
        -DCMAKE_Fortran_COMPILER="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
        -DCMAKE_C_FLAGS=-fPIC \
        -DCMAKE_Fortran_FLAGS=-fPIC \
        -DNetCDF_PATH="${SERIAL_PREFIX}" \
        -DWITH_PNETCDF=FALSE \
        -DPIO_USE_MPIIO=FALSE \
        -DPIO_USE_MPISERIAL=TRUE \
        -DMPISERIAL_PATH="${MPISERIAL_PREFIX}" \
        -DMPI_C_INCLUDE_PATH="${MPISERIAL_PREFIX}/include" \
        -DMPI_C_LIBRARIES="${MPISERIAL_PREFIX}/lib/libmpi-serial.a" \
        -DMPI_Fortran_INCLUDE_PATH="${MPISERIAL_PREFIX}/include" \
        -DMPI_Fortran_LIBRARIES="${MPISERIAL_PREFIX}/lib/libmpi-serial.a" \
        -DPIO_ENABLE_FORTRAN=ON \
        -DPIO_ENABLE_TESTS=OFF \
        -DPIO_ENABLE_TIMING=OFF \
        -DPIO_ENABLE_LOGGING=OFF \
        -DPIO_ENABLE_EXAMPLES=OFF \
        -DBUILD_SHARED_LIBS=OFF \
        -DGENF90_PATH="${src}/scripts" \
        -DUSER_CMAKE_MODULE_PATH=/tmp/CMake_Fortran_utils \
        .. \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/ParallelIO-* /tmp/pio*.tar.gz /tmp/CMake_Fortran_utils

# Named individually so a partial install fails the image build rather than a
# later case link. libpiof.a and libpioc.a are what Tools/Makefile:667 names by
# path; pio.mod is what the CTSM sources USE.
RUN set -eu; \
    for f in "${SERIAL_PREFIX}/lib/libpiof.a" \
             "${SERIAL_PREFIX}/lib/libpioc.a" \
             "${SERIAL_PREFIX}/include/pio.mod"; do \
        [ -f "$f" ] || { echo "ERROR: $f missing after PIO install" >&2; exit 1; }; \
    done

# ---- ESMF, the only flavor (ESMF_COMM=mpiuni, BOPT=O) ----
# run_tests.py forces mpilib="mpi-serial", and so does every case this image
# builds, so executables link with plain gfortran and no MPI library. An
# ESMF_COMM=mpich ESMF linked into one of those fails -- this is what the image
# used to carry, and the failure it produced:
#
#   /usr/bin/ld: /usr/local/esmf-8.6.0-debug/lib/libesmf.so: undefined
#       reference to symbol 'MPI_Bcast'
#   /usr/bin/ld: /opt/mpich/3.4.3/lib/libmpi.so.12: error adding symbols:
#       DSO missing from command line
#
# CTSM cannot avoid pulling ESMF into these links: src/CMakeLists.txt does
# find_package(ESMF REQUIRED) and link_libraries(esmf) unconditionally.
#
# derecho keeps a separate ESMF for exactly this and this mirrors it. Its
# config_machines.xml loads the plain esmf module (esmf/8.9.1) -- NOT the
# esmf-mpi/8.9.1-debug it loads for MPI builds -- whenever
# mpilib="mpi-serial", even under DEBUG="TRUE", and that install's esmf.mk
# reports ESMF_COMM=mpiuni, ESMF_BOPT=O, and a link line with no -lmpi. So
# mpiuni + BOPT=O here; the debug mismatch is deliberate, matching derecho --
# a DEBUG=TRUE case in this container still gets this optimized ESMF, exactly
# as a DEBUG=TRUE mpi-serial case on derecho does.
#
# Both the baked-in ENV ESMFMKFILE at the bottom of this file and the
# set(ESMFMKFILE ...) in cime-macros/gnu_container.cmake now name THIS install;
# the macro's copy is no longer overriding a different default, it just keeps
# the macro self-contained.
# ESMF asks nc-config for "--libs" with no --static (build/common.mk:1583),
# which against this static prefix returns a bare -lnetcdf, for the reason
# given at the netcdf-fortran step above. ESMF skips that query entirely when
# ESMF_NETCDF_LIBS comes from the environment, so supply the complete list.
RUN cd /tmp \
    && wget -nv "https://github.com/esmf-org/esmf/archive/refs/tags/v${ESMF_VERSION}.tar.gz" \
    && tar xzf "v${ESMF_VERSION}.tar.gz" \
    && cd "esmf-${ESMF_VERSION}" \
    && sed -i 's|make VERBOSE=1 -f ./Makefile|$(MAKE) VERBOSE=1 -f ./Makefile|; s|make -f ./Makefile install|$(MAKE) -f ./Makefile install|' \
        src/Infrastructure/IO/PIO/makefile \
    && export PATH="${SERIAL_PREFIX}/bin:${PATH}" \
    && export ESMF_NETCDF_LIBS="$(${SERIAL_PREFIX}/bin/nc-config --libs --static \
            | tr ' ' '\n' | grep '^-l' | tr '\n' ' ')" \
    && export ESMF_DIR="$PWD" ESMF_COMPILER=gfortran ESMF_COMM=mpiuni ESMF_BOPT=O \
        ESMF_NETCDF=nc-config \
        ESMF_PIO=external \
        ESMF_PIO_INCLUDE="${SERIAL_PREFIX}/include" \
        ESMF_PIO_LIBPATH="${SERIAL_PREFIX}/lib" \
        ESMF_PIO_LIBS="-lpioc -lmpi-serial" \
        ESMF_INSTALL_PREFIX="/usr/local/esmf-${ESMF_VERSION}-mpiuni" \
        ESMF_INSTALL_HEADERDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/include" \
        ESMF_INSTALL_MODDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/mod" \
        ESMF_INSTALL_LIBDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/lib" \
        ESMF_INSTALL_BINDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/bin" \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/esmf-* /tmp/v${ESMF_VERSION}.tar.gz

# Encode the failures this flavor exists to prevent, at image-build time
# rather than at test link/run time: it must stay mpiuni, must not link -lmpi
# (serial unit-test links would fail), and must have PIO -- via ESMF_PIO=external,
# because ESMF force-disables an internal PIO for mpiuni.
#
# The -lmpi test splits the link line into tokens and matches EXACTLY, because
# the link line legitimately contains -lmpi-serial now (ESMF_PIO_LIBS). A
# `grep -w -- -lmpi` matches inside "-lmpi-serial" -- '-' is not a word
# character -- and failed this build for an hour on a false positive.
#
# The ESMF_COMM grep is deliberately UNANCHORED. ESMF emits no line-start
# ESMF_COMM= assignment in esmf.mk: the string appears only inside the
# ESMF_*COMPILECPPFLAGS lines (-DESMF_COMM=mpiuni) and in trailing comments.
# Verified against derecho's esmf/8.9.1 as well as 8.6.0. Anchoring this would
# read as a tightening and would instead fail every build. It still
# discriminates, since an MPI build names its own comm, not mpiuni.
RUN set -eu; \
    mk="/usr/local/esmf-${ESMF_VERSION}-mpiuni/lib/esmf.mk"; \
    [ -f "$mk" ] || { echo "ERROR: $mk missing" >&2; exit 1; }; \
    grep -q 'ESMF_COMM=mpiuni' "$mk" \
        || { echo "ERROR: $mk is not an ESMF_COMM=mpiuni build" >&2; exit 1; }; \
    if grep -E '^ESMF_F90ESMFLINKLIBS=' "$mk" | tr ' ' '\n' | grep -qx -- '-lmpi'; then \
        echo "ERROR: $mk still links -lmpi; serial unit-test links will fail" >&2; \
        exit 1; \
    fi; \
    grep -q -- '-DESMF_PIO=1' "$mk" \
        || { echo "ERROR: $mk was built WITHOUT PIO. CDEPS calls" >&2; \
             echo "       ESMF_MeshCreateFromFile, which needs it; an mpi-serial run" >&2; \
             echo "       would abort in ATM init with 'Library needed by ESMF not" >&2; \
             echo "       present', logged only to PET0.ESMF_LogFile." >&2; exit 1; }; \
    grep -qE '^#  ?ESMF_PIO: *external' "$mk" \
        || { echo "ERROR: expected ESMF_PIO=external in $mk. Internal PIO does" >&2; \
             echo "       not support mpiuni -- it fails on a missing mpi.h." >&2; exit 1; }; \
    echo "serial ESMF OK: $mk"

# ---- pFUnit (noMPI, noOpenMP) ----
# noMPI/noOpenMP matches the prebuilt derecho/casper installs
# (pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP) and is what CIME wants
# anyway: run_tests.py forces mpilib="mpi-serial" unconditionally
# (cime/scripts/fortran_unit_testing/run_tests.py, "In the switch from pFUnit3
# to pFUnit4, we have dropped support for MPI").
#
# The compiler necessarily deviates from derecho: ccs_config has a PFUNIT_PATH
# only for intel_derecho/intel_casper, so there is no gnu pFUnit on derecho to
# match. This is the same gfortran as the rest of the image (GCC_VERSION).
#
# Source is the release .tar ASSET, not the GitHub auto-generated source
# tarball and not a git clone. It is the only one that ships pFUnit's nested
# submodules already populated (extern/fArgParse -> extern/gFTL-shared ->
# extern/gFTL), so the build needs no network of its own beyond this fetch.
#
# CC/CXX/FC are named for this build only. They were an override of the MPI
# wrappers the MPICH block used to bake in (ENV CC=mpicc / FC=mpifort); that is
# gone, but naming them still pins this build to the /opt toolchain on PATH
# rather than leaving CMake to probe. CMake caches the compilers at configure
# time, so the setting does not need to extend to the make step.
#
# CMAKE_BUILD_TYPE is deliberately left unset: pFUnit's cmake/BuildType.cmake
# defaults it to Debug, which is also what run_tests.py builds by default
# (debug = not build_optimized).
#
# pFUnit installs FOUR sibling packages under this prefix -- PFUNIT-4.8,
# GFTL-1.11, GFTL_SHARED-1.7, FARGPARSE-1.6 -- and the prefix itself, not the
# PFUNIT-4.8 subdirectory, is what CIME must be pointed at; see cime-macros/.
RUN cd /tmp \
    && wget -nv "https://github.com/Goddard-Fortran-Ecosystem/pFUnit/releases/download/v${PFUNIT_VERSION}/pFUnit-v${PFUNIT_VERSION}.tar" \
    && tar xf "pFUnit-v${PFUNIT_VERSION}.tar" \
    && cd "pFUnit-v${PFUNIT_VERSION}" \
    && mkdir build && cd build \
    && CC=gcc CXX=g++ FC=gfortran F77=gfortran \
       cmake .. \
        -DCMAKE_INSTALL_PREFIX="/usr/local/pfunit-${PFUNIT_VERSION}" \
        -DSKIP_MPI=YES \
        -DSKIP_OPENMP=YES \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/pFUnit-v${PFUNIT_VERSION}*

# ---- CIME macro drop-in that points run_tests.py at the install above ----
# See the file itself for why PFUNIT_PATH cannot simply be an env var or a
# -D flag. Baked at a fixed path AND dropped in /root/.cime: the fixed path is
# what a GitHub Actions job copies from where HOME is /github/home (GitHub's
# standard container handling), while /root/.cime serves an interactive
# `podman run` with no setup -- and also CI on the CIRRUS runner, where HOME
# was measured as /root (2026-10-16).
COPY cime-macros/gnu_container.cmake /opt/ctsm-container/cime-macros/
RUN mkdir -p /root/.cime \
    && cp /opt/ctsm-container/cime-macros/gnu_container.cmake /root/.cime/

# Fail the build, not the test run, if the hardcoded default in the macro ever
# drifts from where pFUnit actually landed (e.g. a PFUNIT_VERSION bump).
# All four must resolve from the single prefix the macro names, because that
# is exactly what find_package does at test-build time under CTSM's
# CMP0074-unset CMakeLists. Checking only PFUNIT here would pass while the
# real build still failed on a missing GFTL.
RUN set -eu; \
    macro=/opt/ctsm-container/cime-macros/gnu_container.cmake; \
    d="$(sed -n 's/^ *set(PFUNIT_PATH "\([^"]*\)").*/\1/p' "$macro" | tail -1)"; \
    echo "macro default PFUNIT_PATH=$d"; \
    for pkg in PFUNIT GFTL GFTL_SHARED FARGPARSE; do \
        cfg="$(ls -d "$d"/${pkg}-*/cmake/${pkg}Config.cmake 2>/dev/null | head -1)"; \
        if [ -z "$cfg" ]; then \
            echo "ERROR: $macro points PFUNIT_PATH at $d, but there is no" >&2; \
            echo "       $d/${pkg}-*/cmake/${pkg}Config.cmake. PFUNIT_PATH must be" >&2; \
            echo "       the shared install prefix holding all four packages." >&2; \
            exit 1; \
        fi; \
        echo "  $pkg -> $cfg"; \
    done

# Same guard for the serial ESMF the macro selects when MPILIB=mpi-serial.
# set(ESMFMKFILE ...) in the macro hardcodes a version, and the macro is COPYed
# in verbatim, so an ESMF_VERSION bump that forgets to edit it leaves the macro
# pointing at a prefix that no longer exists. The esmf.mk assertions further up
# cannot catch that: they check the *templated* path, which a bump moves in
# step. Without this the image builds clean and the failure surfaces much
# later, in run-unit-tests-in-container.sh and in every mpi-serial case.
RUN set -eu; \
    macro=/opt/ctsm-container/cime-macros/gnu_container.cmake; \
    mk="$(sed -n 's/^ *set(ESMFMKFILE "\([^"]*\)").*/\1/p' "$macro" | tail -1)"; \
    echo "macro ESMFMKFILE=$mk"; \
    if [ -z "$mk" ]; then \
        echo "ERROR: no set(ESMFMKFILE ...) found in $macro." >&2; \
        exit 1; \
    fi; \
    if [ ! -f "$mk" ]; then \
        echo "ERROR: $macro points ESMFMKFILE at $mk, which does not exist." >&2; \
        echo "       ESMF_VERSION is ${ESMF_VERSION}; the serial ESMF installed" >&2; \
        echo "       at /usr/local/esmf-${ESMF_VERSION}-mpiuni. Update the" >&2; \
        echo "       hardcoded path in the macro to match the ARG." >&2; \
        exit 1; \
    fi; \
    grep -q 'ESMF_COMM=mpiuni' "$mk" \
        || { echo "ERROR: $mk is not an ESMF_COMM=mpiuni build" >&2; exit 1; }

# Same idea for MPI_SERIAL_PATH, which is the other path the macro hardcodes
# against this image. Get it wrong and nothing fails here: the case build
# quietly falls back to compiling CTSM's own libraries/mpi-serial submodule
# (cime/CIME/Tools/Makefile:919-925, the else branch), which is the vacuous
# behaviour item 4b exists to end -- or, worse, half-takes, and dies much later
# in PIO's configure with "Could NOT find MPISERIAL_C". The headers are checked
# as well as the archive, because only the dedicated prefix carries them; see
# the mpi-serial section above.
RUN set -eu; \
    macro=/opt/ctsm-container/cime-macros/gnu_container.cmake; \
    d="$(sed -n 's/^ *set(MPI_SERIAL_PATH "\([^"]*\)").*/\1/p' "$macro" | tail -1)"; \
    echo "macro MPI_SERIAL_PATH=$d"; \
    if [ -z "$d" ]; then \
        echo "ERROR: $macro sets no MPI_SERIAL_PATH, so a case build would" >&2; \
        echo "       compile CTSM's own mpi-serial submodule exactly as" >&2; \
        echo "       before this image existed." >&2; \
        exit 1; \
    fi; \
    if [ ! -f "$d/lib/libmpi-serial.a" ] || [ ! -f "$d/include/mpi.h" ]; then \
        echo "ERROR: $macro points MPI_SERIAL_PATH at $d, which has no" >&2; \
        echo "       lib/libmpi-serial.a or include/mpi.h. Expected" >&2; \
        echo "       ${MPISERIAL_PREFIX}." >&2; \
        exit 1; \
    fi

# And the same for PIO_LIBDIR/PIO_INCDIR, which the macro also hardcodes.
#
# This is NOT covered by the PIO install assertion further up. That one checks
# ${SERIAL_PREFIX}/lib/libpiof.a and friends, so it moves with the prefix and
# stays green while the macro still points somewhere stale. These two are read
# from the macro, exactly as the PFUNIT_PATH, ESMFMKFILE and MPI_SERIAL_PATH
# checks are, so a prefix change that forgets the macro fails HERE.
#
# They are literals in the macro, not $ENV{PIO}: a container has no modules, so
# gnu_container.cmake hardcodes them the way it hardcodes PFUNIT_PATH and
# ESMFMKFILE. ($ENV{PIO} is what the MACHINE files of real machines do --
# ccs_config/machines/derecho/derecho.cmake:8-10, casper.cmake:9-10,
# cirrus.cmake:6-7 -- and container.cmake does not.)
#
# The failure this prevents is cime/CIME/Tools/Makefile:667 turning a stale
# PIO_LIBDIR into a prerequisite with no recipe, and the case build dying with
# "No rule to make target '<stale>/libpiof.a'" after the sharedlib phase has
# already succeeded.
RUN set -eu; \
    macro=/opt/ctsm-container/cime-macros/gnu_container.cmake; \
    libdir="$(sed -n 's/^ *set(PIO_LIBDIR "\([^"]*\)").*/\1/p' "$macro" | tail -1)"; \
    incdir="$(sed -n 's/^ *set(PIO_INCDIR "\([^"]*\)").*/\1/p' "$macro" | tail -1)"; \
    echo "macro PIO_LIBDIR=$libdir"; \
    echo "macro PIO_INCDIR=$incdir"; \
    if [ -z "$libdir" ] || [ -z "$incdir" ]; then \
        echo "ERROR: $macro sets no PIO_LIBDIR and/or PIO_INCDIR, so a case" >&2; \
        echo "       build would compile CTSM's own libraries/parallelio" >&2; \
        echo "       submodule instead of linking this image's PIO." >&2; \
        exit 1; \
    fi; \
    for f in "$libdir/libpiof.a" "$libdir/libpioc.a" "$incdir/pio.mod"; do \
        [ -f "$f" ] || { \
            echo "ERROR: $macro points PIO_LIBDIR/PIO_INCDIR at" >&2; \
            echo "       $libdir / $incdir, where $f does not exist." >&2; \
            echo "       PIO installed under ${SERIAL_PREFIX}; update the" >&2; \
            echo "       hardcoded paths in the macro to match." >&2; \
            exit 1; }; \
        echo "  $f"; \
    done

# ---- make the installed .pc files describe a STATIC link ----
# MEASURED on the first rebuild after item 4b (2026-10-15): cprnc died at link
# with ~200 undefined H5* symbols, and the link line was
#   f95 ... -o cprnc /usr/local/lib/libnetcdff.a /usr/local/lib/libnetcdf.a
# -- the two archives and nothing else.
#
# That was cprnc v1.1.2, whose CMakeLists.txt:78 tries find_package(NetCDF),
# which cannot succeed against this stack (autotools-built, so no
# netCDFConfig.cmake, and cime ships no FindNetCDF.cmake), and falls through to
# pkg_check_modules at :91-101. The build log said so in as many words: "NetCDF
# was not found, falling back on pkg-config". pkg_check_modules fills NETCDF_LDFLAGS from the
# .pc file's Libs: field alone -- it puts the --static answer in a separate
# NETCDF_STATIC_LDFLAGS that cprnc never reads -- while netCDF's configure
# records the transitive -lhdf5_hl -lhdf5 -lz -lcurl -lxml2 in Libs.private.
#
# Before 4b this directory held the SHARED parallel stack, where none of that
# mattered: DT_NEEDED on libnetcdf.so pulled HDF5 in at link time. A static
# prefix has no such mechanism, so these .pc files are simply WRONG about how
# to link the stack they describe. Fix them rather than patching a consumer:
# every pkg-config user gets the right answer, and cime stays unmodified.
#
# THIS NO LONGER FIXES A CASE BUILD'S cprnc. cime6.5.5 pins cprnc v1.1.5,
# which tries find_library(netcdf) BEFORE pkg-config; that finds the bare
# static archives in /usr/local/lib and links them with none of their
# dependencies, so the .pc files are never consulted and the link fails exactly
# as above (measured in CI, 2026-10-16; NEXT_STEPS item 13). The fix for that
# is to build cprnc once in the image and point CCSM_CPRNC at it -- being
# proven in Dockerfile.cprnc-probe, whose cprnc build takes its link libraries
# from these .pc files. Keep this layer for that and for any other pkg-config
# user.
#
# Private libraries are appended AFTER the existing Libs:, which is the order a
# static link needs -- dependencies follow the things that depend on them.
#
# Deliberately placed this late, after everything expensive: it rewrites two
# text files and depends only on the netCDF installs, so a cache-warm rebuild
# redoes this layer and nothing before it.
RUN set -eu; \
    for pc in "${SERIAL_PREFIX}/lib/pkgconfig/netcdf.pc" \
              "${SERIAL_PREFIX}/lib/pkgconfig/netcdf-fortran.pc"; do \
        [ -f "$pc" ] || { echo "ERROR: $pc missing" >&2; exit 1; }; \
        priv="$(sed -n 's/^Libs\.private:[[:space:]]*//p' "$pc")"; \
        if [ -z "$priv" ]; then \
            echo "note: $pc carries no Libs.private; nothing to fold"; \
            continue; \
        fi; \
        awk -v priv="$priv" '/^Libs:/ { print $0 " " priv; next } { print }' \
            "$pc" > "$pc.new"; \
        mv "$pc.new" "$pc"; \
    done; \
    export PKG_CONFIG_PATH="${SERIAL_PREFIX}/lib/pkgconfig"; \
    echo "netcdf         ->$(pkg-config --libs netcdf)"; \
    echo "netcdf-fortran ->$(pkg-config --libs netcdf-fortran)"; \
    pkg-config --libs netcdf | tr ' ' '\n' | grep -qx -- '-lhdf5' \
        || { echo "ERROR: 'pkg-config --libs netcdf' still omits -lhdf5, so a" >&2; \
             echo "       cprnc link would fail on undefined H5* symbols." >&2; \
             exit 1; }; \
    pkg-config --libs netcdf-fortran | tr ' ' '\n' | grep -qx -- '-lnetcdf' \
        || { echo "ERROR: 'pkg-config --libs netcdf-fortran' omits -lnetcdf." >&2; \
             exit 1; }
# The two greps split into tokens and match EXACTLY, because -lnetcdff contains
# -lnetcdf as a substring and '-' is not a word character, so a plain grep -w
# would pass on the wrong flag. Same trap as the -lmpi/-lmpi-serial test at the
# serial ESMF assertion above.

# Runtime environment for CIME's "container" machine:
# - ESMFMKFILE: not set by the machine config; must come from the image. There
#   is one ESMF now, the ESMF_COMM=mpiuni one, so this and the macro's
#   set(ESMFMKFILE ...) name the same file.
# - PIO / PIO_VERSION_MAJOR / PIO_TYPENAME_VALID_VALUES: all three come from
#   derecho's parallelio-serial modulefile, and CIME expects all three together
#   -- cime/CIME/build_scripts/buildlib.pio:76-81 says so in as many words, and
#   82-91 is the code. PIO_VERSION_MAJOR is the one that does the work:
#   buildlib.pio reads it from the ENVIRONMENT and returns before running any
#   make, so PIO is never built from source. It cannot come from the macro -- a
#   cmake macro file produces make variables, and os.environ does not see them
#   -- and setting PIO_LIBDIR there does not substitute, because
#   buildlib.pio:115 passes its own PIO_LIBDIR on the gmake command line, which
#   beats anything Macros.make carries. PIO_TYPENAME_VALID_VALUES is not
#   optional polish either, and it fails ugly: line 90 passes
#   os.environ.get("PIO_TYPENAME_VALID_VALUES") straight into
#   _set_pio_valid_values, which does valid_values += ",nothing" at line 224 and
#   raises TypeError on None -- inside the very branch these variables exist to
#   reach. "netcdf" is the honest value: this stack has no pnetcdf and no
#   parallel netCDF4.
#     PIO itself is read by NOTHING IN THIS IMAGE, and that is deliberate
#   rather than an oversight. buildlib.pio reads only the other two
#   (os.environ at lines 82 and 90), and the $ENV{PIO} that derives
#   PIO_LIBDIR/PIO_INCDIR belongs to real machines' ccs_config files --
#   derecho.cmake:8-10, casper.cmake:9-10, cirrus.cmake:6-7 -- not to
#   container.cmake. gnu_container.cmake hardcodes both paths as literals
#   instead, because a container has no modules to export PIO in the first
#   place; its only $ENV{} read is PFUNIT_PATH. PIO is exported anyway because
#   it is the variable derecho's parallelio-serial module sets, so the image
#   satisfies the same contract a machine with modules would, and the day
#   ccs_config's container.cmake grows the derecho.cmake lines the macro can
#   drop its hardcoding with nothing else to change. The assertion above, not
#   this variable, is what actually guards those two paths.
# - CC / CXX / FC: which toolchain a CMake-driven build inside the container
#   uses. CIME's Makefile builds set their own from $(SCC)/$(SFC). cprnc is
#   built by cmake (cime/CIME/build_scripts/buildlib.cprnc), and the cmake
#   command line it generates passes no -DCMAKE_*_COMPILER at all. Since cprnc
#   v1.1.5 (cime6.5.5) cprnc's own CMakeLists.txt sets the compilers from the
#   case's Macros.cmake before project(), so what follows -- measured with
#   v1.1.2 -- no longer applies to cprnc; it still applies to any other CMake
#   build run in the container without explicit compilers. CMake then
#   falls back to its built-in name search -- `cc` for C, `f95` for Fortran --
#   and NEITHER NAME EXISTS in /opt/gcc/${GCC_VERSION}/bin, which ships `gcc`
#   and `gfortran`. PATH order cannot help: cmake searches by name first, so it
#   lands on /usr/bin/cc and /usr/bin/f95, AlmaLinux's bootstrap GCC 11.5.0.
#
#   MEASURED 2026-10-15, the rebuild after the .pc fix above: cprnc linked the
#   static HDF5 correctly and then died on undefined __truncsfhf2 /
#   __truncdfhf2 / __truncxfhf2 from libhdf5.a(H5Tconv_float.o). Those are
#   libgcc's _Float16 conversion helpers. HDF5 here is compiled by GCC
#   ${GCC_VERSION}, which has _Float16 on x86-64; GCC 11.5.0 does not (it
#   arrived in GCC 12), so its libgcc cannot resolve them. CMakeCache.txt
#   recorded CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/f95.
#
#   Before item 4b the image set ENV CC=mpicc/CXX=mpicxx/FC=mpifort, which
#   CMake reads as $ENV{CC}/$ENV{FC}; those wrapped the same GCC
#   ${GCC_VERSION}, so cprnc was built with the right toolchain by accident of
#   the MPI stack. 4b deleted that ENV with the stack, and the bare -lnetcdf
#   failure masked the consequence until it was fixed. Absolute paths rather
#   than bare names so this cannot depend on PATH.
# - CESMDATAROOT: DIN_LOC_ROOT=$CESMDATAROOT/inputdata must exist for
#   namelist generation even in build-only (--no-run) testing.
# - USER: CIME requires it; GitHub Actions container jobs don't set it.
# - PKG_CONFIG_PATH: EL9's pkg-config search path omits
#   /usr/local/lib/pkgconfig, and CIME's cprnc build locates
#   netCDF/netCDF-Fortran via pkg-config. That really is the path taken:
#   cime/CIME/non_py/cprnc/CMakeLists.txt:78 tries find_package(NetCDF) first,
#   but this stack is autotools-built (no netCDFConfig.cmake) and cime ships no
#   FindNetCDF.cmake, so the pkg_check_modules fallback at :91-101 is what runs.
#   These .pc files describe the STATIC serial stack, and as netCDF installs
#   them `pkg-config --libs netcdf` yields a bare -lnetcdf -- which broke cprnc
#   on the first rebuild after 4b. The layer above rewrites them; see it for
#   the measurement and the mechanism.
# - PKG_CONFIG_ALLOW_SYSTEM_CFLAGS: EL9's pkgconf filters
#   -I/usr/local/include out of --cflags as a "system" dir, but that is
#   only true for the system gcc; the /opt toolchain does not search
#   /usr/local/include by default.
ENV ESMFMKFILE=/usr/local/esmf-${ESMF_VERSION}-mpiuni/lib/esmf.mk \
    PIO=${SERIAL_PREFIX} \
    PIO_VERSION_MAJOR=2 \
    PIO_TYPENAME_VALID_VALUES=netcdf \
    CC=/opt/gcc/${GCC_VERSION}/bin/gcc \
    CXX=/opt/gcc/${GCC_VERSION}/bin/g++ \
    FC=/opt/gcc/${GCC_VERSION}/bin/gfortran \
    CESMDATAROOT=/opt/cesmdata \
    USER=root \
    PKG_CONFIG_PATH=/usr/local/lib/pkgconfig \
    PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1

# Fail the build, not a later cprnc link, if CC/CXX/FC name something that is
# not there. A GCC_VERSION bump moves /opt/gcc/<ver>/bin, and these are the
# only paths in the file that a bump could leave dangling without any earlier
# step noticing -- every other consumer of the toolchain either uses it during
# the build (and so fails then) or goes through PATH.
RUN set -eu; \
    for v in CC CXX FC; do \
        eval "p=\${$v}"; \
        [ -x "$p" ] || { \
            echo "ERROR: $v=$p is not an executable. GCC_VERSION is" >&2; \
            echo "       ${GCC_VERSION}; update the CC/CXX/FC ENV to match." >&2; \
            exit 1; }; \
        echo "$v=$p -> $("$p" --version | head -1)"; \
    done

RUN mkdir -p "${CESMDATAROOT}/inputdata"

# ---- build-time record of what the installed stack reports ----
# netCDF-Fortran's version and netCDF-C's filter set cannot be recovered from a
# container cesm.exe: this stack is static, so `ldd` is blind to the filter
# libraries, and nothing stamps the netCDF-Fortran version into the executable.
# DRIFT_CHECK_DESIGN.md's instrument map marks both as gaps. Item 5's drift
# check reads them back out of this file with `podman run ... cat`.
#
# A FILE, NOT A LABEL, deliberately. NEXT_STEPS item 4b asks for LABELs, but
# LABEL takes a literal string -- it cannot run nf-config -- and a LABEL fed
# from a build ARG would only repeat what this Dockerfile already says instead
# of what the build produced, which is the whole point of these two values. Do
# not "fix" this back into a LABEL.
#
# The --has-stdfilters option is confirmed to exist before it is used: some
# builds of nc-config answer an unknown option with usage text on stdout and
# exit 0, which would get recorded as if it were the filter list.
#
# AN EMPTY FILTER LIST IS NOT AN ERROR and must not fail the build. It is a
# legitimate answer -- it means no standard filters beyond what is compiled in,
# and this image has only zlib where derecho carries szip, c-blosc, zstd,
# bzip2 and lz4 (DRIFT_CHECK_DESIGN.md's instrument map). That difference is
# the kind of thing item 5 exists to REPORT, so recording it honestly is the
# whole job; failing here would turn a real finding into a broken build, at the
# very last layer, after everything expensive. The --help guard above already
# covers the case this was defending against. nf-config --version is different:
# empty there genuinely is wrong, so it stays required.
RUN set -eu; \
    nc-config --help 2>&1 | grep -q -- '--has-stdfilters' \
        || { echo "ERROR: this nc-config has no --has-stdfilters, which item" >&2; \
             echo "       5's drift check needs; see DRIFT_CHECK_DESIGN.md." >&2; \
             exit 1; }; \
    nfv="$(nf-config --version)"; \
    ncf="$(nc-config --has-stdfilters)"; \
    [ -n "$nfv" ] || { echo "ERROR: nf-config --version returned nothing" >&2; exit 1; }; \
    mkdir -p /opt/ctsm-container; \
    { echo "netcdf_fortran_version=${nfv}"; \
      echo "netcdf_c_stdfilters=${ncf}"; } > /opt/ctsm-container/stack-facts.txt; \
    cat /opt/ctsm-container/stack-facts.txt

# mpi-serial and pio are labelled because since item 4b they are what a case
# build LINKS, not just image internals, and item 5's drift check reads the
# container side of both "from the ARG" (DRIFT_CHECK_DESIGN.md's instrument
# map) -- so the ARG has to be readable from the image.
#
# Unlike the other labels, these two cannot be cross-checked against what the
# image reports, and smoke-test.sh deliberately does not try:
#   - Nothing in the image reports an mpi-serial version. The step above copies
#     libmpi-serial.a, mpi.h and mpif.h by hand and nothing else; there is no
#     config script and no version file.
#   - PIO installs lib/libpio.settings ("PIO Version:"), and that is the ONLY
#     place a built PIO states its version: pio_meta.h is generated into the
#     build tree (CMakeLists.txt:507-511), never installed, and read only by
#     tests/cunit, which this image does not build. But UPSTREAM HAS NOT BUMPED
#     THE IN-TREE VERSION SINCE 2.6.6:
#     CMakeLists.txt:11-13 still sets VERSION_PATCH 6 at the pio2_6_7 and
#     pio2_6_8 tags (checked against all three tags on 2026-10-15), and
#     derecho's own Spack-built parallelio/2.6.8 reads "PIO Version: 2.6.6"
#     in its libpio.settings -- so this is upstream, not a container
#     artifact. A build of
#     ${PIO_VERSION}=2.6.8 therefore reports 2.6.6 about itself, and a check of
#     this label against libpio.settings would fail a CORRECT image. Do not add
#     one, and do not "correct" PIO_VERSION to 2.6.6 -- the tag the Dockerfile
#     fetches is the honest version, and it is what check-derecho-versions.py
#     compares to derecho's parallelio-serial module.
LABEL gcc="${GCC_VERSION}" \
      hdf5="${HDF5_VERSION}" \
      netcdf-c="${NETCDF_C_VERSION}" \
      netcdf-fortran="${NETCDF_FORTRAN_VERSION}" \
      esmf="${ESMF_VERSION}" \
      pfunit="${PFUNIT_VERSION}" \
      mpi-serial="${MPI_SERIAL_VERSION}" \
      pio="${PIO_VERSION}" \
      org.opencontainers.image.description="CTSM build/test image for GitHub Actions CI, matching derecho's gnu stack" \
      org.opencontainers.image.source="https://github.com/ESCOMP/CTSM"

CMD ["/bin/bash"]
