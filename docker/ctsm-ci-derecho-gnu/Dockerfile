# Container image for building CTSM in GitHub Actions CI.
#
# Built from scratch on AlmaLinux 9 so every component can match what
# CIME loads for gnu builds on derecho (per
# ccs_config/machines/derecho/config_machines.xml) exactly where an open
# equivalent exists; see README.md for the version-match table and the
# accepted deviations (cray-mpich -> MPICH 3.4.3, cray-libsci ->
# reference LAPACK/BLAS).
#
# Libraries are installed under /usr/local because CIME's "container"
# machine (ccs_config/machines/container/) hard-codes
# NETCDF_PATH=/usr/local and PNETCDF_PATH=/usr/local. The toolchain
# lives under /opt. ESMF is located via ESMFMKFILE instead, so its three
# build flavors (debug, optimized, and a serial ESMF_COMM=mpiuni build for
# unit tests) live in their own prefixes.
#
# The image serves two CIME workflows: create_test builds (mpilib=mpich) and
# CTSM's Fortran unit tests (cime/scripts/fortran_unit_testing/run_tests.py,
# which forces mpilib=mpi-serial). The unit-test path additionally needs
# pFUnit and the cime-macros drop-in near the bottom of this file.
#
# All environment needed at runtime is baked in with ENV: GitHub
# Actions runs steps in non-login shells, so profile.d scripts never
# fire in CI.

FROM docker.io/almalinux:9

ARG GCC_VERSION=12.2.0
ARG MPICH_VERSION=3.4.3
ARG HDF5_VERSION=1.12.2
ARG NETCDF_C_VERSION=4.9.2
ARG NETCDF_FORTRAN_VERSION=4.6.1
ARG PNETCDF_VERSION=1.12.3
ARG ESMF_VERSION=8.6.0
ARG PFUNIT_VERSION=4.8.0

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

# ---- GCC (matches derecho's gcc/12.2.0 module) ----
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
    && wget -nv "https://ftp.gnu.org/gnu/gcc/gcc-${GCC_VERSION}/gcc-${GCC_VERSION}.tar.xz" \
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

ENV PATH=/usr/local/bin:/opt/mpich/${MPICH_VERSION}/bin:/opt/gcc/${GCC_VERSION}/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin \
    LD_LIBRARY_PATH=/usr/local/lib:/opt/gcc/${GCC_VERSION}/lib64:/opt/mpich/${MPICH_VERSION}/lib

# ---- MPICH (open equivalent of derecho's cray-mpich 8.1.27, which is
# MPICH-3.4-ABI-derived). ch4:ofi matches cray-mpich's libfabric
# lineage; the libfabric implementation is the one embedded in the
# MPICH tarball. ----
RUN cd /tmp \
    && wget -nv "https://www.mpich.org/static/downloads/${MPICH_VERSION}/mpich-${MPICH_VERSION}.tar.gz" \
    && tar xzf "mpich-${MPICH_VERSION}.tar.gz" \
    && cd "mpich-${MPICH_VERSION}" \
    && CC="/opt/gcc/${GCC_VERSION}/bin/gcc" \
       CXX="/opt/gcc/${GCC_VERSION}/bin/g++" \
       F77="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
       FC="/opt/gcc/${GCC_VERSION}/bin/gfortran" \
       FFLAGS="-fallow-argument-mismatch" \
       FCFLAGS="-fallow-argument-mismatch" \
       ./configure --prefix="/opt/mpich/${MPICH_VERSION}" \
        --with-device=ch4:ofi --enable-shared \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && printf '%s\n' /usr/local/lib "/opt/mpich/${MPICH_VERSION}/lib" \
        > /etc/ld.so.conf.d/ctsm.conf \
    && ldconfig \
    && rm -rf /tmp/mpich-*

ENV CC=mpicc \
    CXX=mpicxx \
    FC=mpifort \
    F77=mpif77

# ---- HDF5 (parallel) ----
# Only the git tag archive is still downloadable for this version; it
# includes the generated configure script.
RUN cd /tmp \
    && wget -nv "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-${HDF5_VERSION//./_}.tar.gz" \
    && tar xzf "hdf5-${HDF5_VERSION//./_}.tar.gz" \
    && cd "hdf5-hdf5-${HDF5_VERSION//./_}" \
    && ./configure --prefix=/usr/local --enable-parallel --enable-hl --with-zlib \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/hdf5-*

# ---- netCDF-C ----
# libxml2-devel is installed here rather than in the base dnf layer to
# keep the (expensive) GCC layer cache valid; fold it into the base
# layer whenever a clean rebuild happens anyway.
RUN dnf -y install libxml2-devel && dnf clean all \
    && cd /tmp \
    && wget -nv "https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF_C_VERSION}/netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
    && tar xzf "netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
    && cd "netcdf-c-${NETCDF_C_VERSION}" \
    && CPPFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib" \
       ./configure --prefix=/usr/local --enable-netcdf-4 --enable-shared \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/netcdf-c-*

# ---- netCDF-Fortran ----
RUN cd /tmp \
    && wget -nv "https://downloads.unidata.ucar.edu/netcdf-fortran/${NETCDF_FORTRAN_VERSION}/netcdf-fortran-${NETCDF_FORTRAN_VERSION}.tar.gz" \
    && tar xzf "netcdf-fortran-${NETCDF_FORTRAN_VERSION}.tar.gz" \
    && cd "netcdf-fortran-${NETCDF_FORTRAN_VERSION}" \
    && CPPFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib" \
       ./configure --prefix=/usr/local --enable-shared \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/netcdf-fortran-*

# ---- PnetCDF ----
RUN cd /tmp \
    && wget -nv "https://parallel-netcdf.github.io/Release/pnetcdf-${PNETCDF_VERSION}.tar.gz" \
    && tar xzf "pnetcdf-${PNETCDF_VERSION}.tar.gz" \
    && cd "pnetcdf-${PNETCDF_VERSION}" \
    && ./configure --prefix=/usr/local --enable-shared MPICC=mpicc MPICXX=mpicxx MPIF77=mpif77 MPIF90=mpif90 \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/pnetcdf-*

# ---- ESMF, debug flavor (BOPT=g) ----
# Derecho's gnu DEBUG builds load esmf/8.6.0-debug; this is the default
# ESMFMKFILE below since the current CI test is SMS_D (debug).
# ESMF_COMM=mpich is ESMF's setting for MPICH v3 and later.
#
# The sed patch fixes an upstream bug in ESMF's bundled-PIO recipe: it
# invokes the cmake-generated sub-make as plain `make` instead of
# `$(MAKE)`, so GNU make doesn't mark it recursive and the jobserver
# file descriptors it inherits via MAKEFLAGS are invalid ("write
# jobserver: Bad file descriptor" under parallel builds).
RUN cd /tmp \
    && wget -nv "https://github.com/esmf-org/esmf/archive/refs/tags/v${ESMF_VERSION}.tar.gz" \
    && tar xzf "v${ESMF_VERSION}.tar.gz" \
    && cd "esmf-${ESMF_VERSION}" \
    && sed -i 's|make VERBOSE=1 -f ./Makefile|$(MAKE) VERBOSE=1 -f ./Makefile|; s|make -f ./Makefile install|$(MAKE) -f ./Makefile install|' \
        src/Infrastructure/IO/PIO/makefile \
    && export ESMF_DIR="$PWD" ESMF_COMPILER=gfortran ESMF_COMM=mpich ESMF_BOPT=g \
        ESMF_NETCDF=nc-config \
        ESMF_INSTALL_PREFIX="/usr/local/esmf-${ESMF_VERSION}-debug" \
        ESMF_INSTALL_HEADERDIR="/usr/local/esmf-${ESMF_VERSION}-debug/include" \
        ESMF_INSTALL_MODDIR="/usr/local/esmf-${ESMF_VERSION}-debug/mod" \
        ESMF_INSTALL_LIBDIR="/usr/local/esmf-${ESMF_VERSION}-debug/lib" \
        ESMF_INSTALL_BINDIR="/usr/local/esmf-${ESMF_VERSION}-debug/bin" \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/esmf-* /tmp/v${ESMF_VERSION}.tar.gz

# ---- ESMF, optimized flavor (BOPT=O) ----
# For future non-debug CI tests; select it by overriding ESMFMKFILE.
# Same PIO-makefile patch as the debug flavor above.
RUN cd /tmp \
    && wget -nv "https://github.com/esmf-org/esmf/archive/refs/tags/v${ESMF_VERSION}.tar.gz" \
    && tar xzf "v${ESMF_VERSION}.tar.gz" \
    && cd "esmf-${ESMF_VERSION}" \
    && sed -i 's|make VERBOSE=1 -f ./Makefile|$(MAKE) VERBOSE=1 -f ./Makefile|; s|make -f ./Makefile install|$(MAKE) -f ./Makefile install|' \
        src/Infrastructure/IO/PIO/makefile \
    && export ESMF_DIR="$PWD" ESMF_COMPILER=gfortran ESMF_COMM=mpich ESMF_BOPT=O \
        ESMF_NETCDF=nc-config \
        ESMF_INSTALL_PREFIX="/usr/local/esmf-${ESMF_VERSION}" \
        ESMF_INSTALL_HEADERDIR="/usr/local/esmf-${ESMF_VERSION}/include" \
        ESMF_INSTALL_MODDIR="/usr/local/esmf-${ESMF_VERSION}/mod" \
        ESMF_INSTALL_LIBDIR="/usr/local/esmf-${ESMF_VERSION}/lib" \
        ESMF_INSTALL_BINDIR="/usr/local/esmf-${ESMF_VERSION}/bin" \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/esmf-* /tmp/v${ESMF_VERSION}.tar.gz

# ---- ESMF, serial flavor (ESMF_COMM=mpiuni, BOPT=O) ----
# run_tests.py forces mpilib="mpi-serial", so the unit-test executables link
# with plain gfortran and no MPI library. Both ESMF flavors in the base image
# are ESMF_COMM=mpich, and linking one of those into a serial executable fails:
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
# config_machines.xml loads esmf/8.6.0 -- NOT esmf/8.6.0-debug -- whenever
# mpilib="mpi-serial", even under DEBUG="TRUE", and that install's esmf.mk
# reports ESMF_COMM=mpiuni, ESMF_BOPT=O, and a link line with no -lmpi. So
# mpiuni + BOPT=O here; the debug mismatch is deliberate, matching derecho.
#
# Selected via ESMFMKFILE from cime-macros/gnu_container.cmake when
# MPILIB=mpi-serial. The image's baked-in ENV ESMFMKFILE still points at the
# debug mpich flavor, so ordinary create_test builds are untouched.
RUN cd /tmp \
    && wget -nv "https://github.com/esmf-org/esmf/archive/refs/tags/v${ESMF_VERSION}.tar.gz" \
    && tar xzf "v${ESMF_VERSION}.tar.gz" \
    && cd "esmf-${ESMF_VERSION}" \
    && sed -i 's|make VERBOSE=1 -f ./Makefile|$(MAKE) VERBOSE=1 -f ./Makefile|; s|make -f ./Makefile install|$(MAKE) -f ./Makefile install|' \
        src/Infrastructure/IO/PIO/makefile \
    && export ESMF_DIR="$PWD" ESMF_COMPILER=gfortran ESMF_COMM=mpiuni ESMF_BOPT=O \
        ESMF_NETCDF=nc-config \
        ESMF_INSTALL_PREFIX="/usr/local/esmf-${ESMF_VERSION}-mpiuni" \
        ESMF_INSTALL_HEADERDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/include" \
        ESMF_INSTALL_MODDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/mod" \
        ESMF_INSTALL_LIBDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/lib" \
        ESMF_INSTALL_BINDIR="/usr/local/esmf-${ESMF_VERSION}-mpiuni/bin" \
    && make -j"${MAKE_JOBS}" \
    && make install \
    && rm -rf /tmp/esmf-* /tmp/v${ESMF_VERSION}.tar.gz

# Encode the failure this flavor exists to prevent: if the serial ESMF ever
# comes out linking MPI again, fail at image-build time, not at test link time.
RUN set -eu; \
    mk="/usr/local/esmf-${ESMF_VERSION}-mpiuni/lib/esmf.mk"; \
    [ -f "$mk" ] || { echo "ERROR: $mk missing" >&2; exit 1; }; \
    grep -q 'ESMF_COMM=mpiuni' "$mk" \
        || { echo "ERROR: $mk is not an ESMF_COMM=mpiuni build" >&2; exit 1; }; \
    if grep -E '^ESMF_F90ESMFLINKLIBS=' "$mk" | grep -qw -- '-lmpi'; then \
        echo "ERROR: $mk still links -lmpi; serial unit-test links will fail" >&2; \
        exit 1; \
    fi; \
    echo "serial ESMF OK: $mk"

# ---- git (built from source; NO openssh dependency) ----
# CIME build steps clone helper repos at build time via CMake
# ExternalProject_Add: cprnc fetches genf90 (github.com/PARALLELIO/genf90),
# and PIO/CDEPS/CMEPS generate code with genf90 too -- so the image needs a
# real git at runtime. git cannot be dnf-installed here: it pulls openssh,
# whose rpm chowns a setuid helper (ssh-keysign) to a non-root id, which
# fails under rootless single-UID podman ("cpio: chown failed"). Building
# from source avoids openssh entirely -- https clones go through libcurl --
# and behaves like a stock git. NO_OPENSSL/NO_EXPAT/NO_GETTEXT/NO_TCLTK drop
# optional deps we don't ship. Placed after the toolchain so editing this
# layer doesn't rebuild those.
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

# ---- pFUnit (noMPI, noOpenMP) ----
# noMPI/noOpenMP matches the prebuilt derecho/casper installs
# (pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP) and is what CIME wants
# anyway: run_tests.py forces mpilib="mpi-serial" unconditionally
# (cime/scripts/fortran_unit_testing/run_tests.py, "In the switch from pFUnit3
# to pFUnit4, we have dropped support for MPI").
#
# The compiler necessarily deviates from derecho: ccs_config has a PFUNIT_PATH
# only for intel_derecho/intel_casper, so there is no gnu pFUnit on derecho to
# match. This is gfortran 12.2.0 -- same version as the rest of the image.
#
# Source is the release .tar ASSET, not the GitHub auto-generated source
# tarball and not a git clone. It is the only one that ships pFUnit's nested
# submodules already populated (extern/fArgParse -> extern/gFTL-shared ->
# extern/gFTL), so the build needs no network of its own beyond this fetch.
#
# CC/CXX/FC are overridden for this build only: the image bakes in the MPI
# wrappers (ENV CC=mpicc / FC=mpifort) and this is deliberately serial. CMake
# caches the compilers at configure time, so the override does not need to
# extend to the make step.
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
# what a GitHub Actions job copies from (GHA overrides HOME, so /root/.cime is
# not found there), while /root/.cime makes an interactive `podman run` work
# with no setup.
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

# Runtime environment for CIME's "container" machine:
# - ESMFMKFILE: not set by the machine config; must come from the image.
# - CESMDATAROOT: DIN_LOC_ROOT=$CESMDATAROOT/inputdata must exist for
#   namelist generation even in build-only (--no-run) testing.
# - USER: CIME requires it; GitHub Actions container jobs don't set it.
# - PKG_CONFIG_PATH: EL9's pkg-config search path omits
#   /usr/local/lib/pkgconfig, and CIME's cprnc build locates
#   netCDF/netCDF-Fortran via pkg-config.
# - PKG_CONFIG_ALLOW_SYSTEM_CFLAGS: EL9's pkgconf filters
#   -I/usr/local/include out of --cflags as a "system" dir, but that is
#   only true for the system gcc; the /opt toolchain does not search
#   /usr/local/include by default.
ENV ESMFMKFILE=/usr/local/esmf-${ESMF_VERSION}-debug/lib/esmf.mk \
    CESMDATAROOT=/opt/cesmdata \
    USER=root \
    PKG_CONFIG_PATH=/usr/local/lib/pkgconfig \
    PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1

RUN mkdir -p "${CESMDATAROOT}/inputdata"

LABEL gcc="${GCC_VERSION}" \
      mpich="${MPICH_VERSION}" \
      hdf5="${HDF5_VERSION}" \
      netcdf-c="${NETCDF_C_VERSION}" \
      netcdf-fortran="${NETCDF_FORTRAN_VERSION}" \
      pnetcdf="${PNETCDF_VERSION}" \
      esmf="${ESMF_VERSION}" \
      pfunit="${PFUNIT_VERSION}" \
      org.opencontainers.image.description="CTSM build/test image for GitHub Actions CI, matching derecho's gnu stack"

CMD ["/bin/bash"]
