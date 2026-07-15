# Container image for building CTSM in GitHub Actions CI.
#
# The software stack matches what CIME loads for gnu builds on derecho
# (per ccs_config/machines/derecho/config_machines.xml) as closely as a
# container allows; see README.md for the version-match table and the
# accepted deviations.
#
# Libraries are installed under /usr/local because CIME's "container"
# machine (ccs_config/machines/container/) hard-codes
# NETCDF_PATH=/usr/local and PNETCDF_PATH=/usr/local. ESMF is located
# via ESMFMKFILE instead, so its two build flavors (debug/optimized)
# live in their own prefixes.
#
# All environment needed at runtime is baked in with ENV: GitHub
# Actions runs steps in non-login shells, so the base image's
# /etc/profile.d/z00-build-env.sh (login shells only) never fires in CI.

FROM docker.io/ncarcisl/cisldev-x86_64-almalinux9-gcc12-mpich3:25.08

ARG GCC_ROOT=/container/gcc/12.5.0
ARG MPICH_ROOT=/container/mpich/3.4.3
ARG HDF5_VERSION=1.12.2
ARG NETCDF_C_VERSION=4.9.2
ARG NETCDF_FORTRAN_VERSION=4.6.1
ARG PNETCDF_VERSION=1.12.3
ARG ESMF_VERSION=8.6.0

# Deterministic toolchain environment for both the build layers below and
# the final image (GHA inherits image ENV). /usr/local/bin comes first so
# nc-config etc. resolve to the versions built here, not the base image's
# newer copies under /container.
ENV PATH=/usr/local/bin:${MPICH_ROOT}/bin:${GCC_ROOT}/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin \
    LD_LIBRARY_PATH=/usr/local/lib:${GCC_ROOT}/lib64:${GCC_ROOT}/lib:${MPICH_ROOT}/lib \
    CC=mpicc \
    CXX=mpicxx \
    FC=mpifort \
    F77=mpif77

SHELL ["/bin/bash", "-c"]

# LAPACK/BLAS stand in for derecho's cray-libsci (proprietary); the CIME
# "container" machine links -llapack -lblas. The Perl modules are needed
# by CTSM's build-namelist.
RUN dnf -y install dnf-plugins-core \
    && dnf config-manager --set-enabled crb \
    && dnf -y install \
        lapack-devel blas-devel \
        perl-core perl-XML-LibXML \
        zlib-devel libcurl-devel m4 patch file diffutils hostname \
    && dnf clean all \
    && printf '%s\n' /usr/local/lib "${GCC_ROOT}/lib64" "${GCC_ROOT}/lib" "${MPICH_ROOT}/lib" \
        > /etc/ld.so.conf.d/ctsm.conf \
    && ldconfig

# ---- HDF5 (parallel) ----
# Only the git tag archive is still downloadable for this version; it
# includes the generated configure script.
RUN cd /tmp \
    && wget -nv "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-${HDF5_VERSION//./_}.tar.gz" \
    && tar xzf "hdf5-${HDF5_VERSION//./_}.tar.gz" \
    && cd "hdf5-hdf5-${HDF5_VERSION//./_}" \
    && ./configure --prefix=/usr/local --enable-parallel --enable-hl --with-zlib \
    && make -j"$(nproc)" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/hdf5-*

# ---- netCDF-C ----
RUN cd /tmp \
    && wget -nv "https://downloads.unidata.ucar.edu/netcdf-c/${NETCDF_C_VERSION}/netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
    && tar xzf "netcdf-c-${NETCDF_C_VERSION}.tar.gz" \
    && cd "netcdf-c-${NETCDF_C_VERSION}" \
    && CPPFLAGS="-I/usr/local/include" LDFLAGS="-L/usr/local/lib" \
       ./configure --prefix=/usr/local --enable-netcdf-4 --enable-shared \
    && make -j"$(nproc)" \
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
    && make -j"$(nproc)" \
    && make install \
    && ldconfig \
    && rm -rf /tmp/netcdf-fortran-*

# ---- PnetCDF ----
RUN cd /tmp \
    && wget -nv "https://parallel-netcdf.github.io/Release/pnetcdf-${PNETCDF_VERSION}.tar.gz" \
    && tar xzf "pnetcdf-${PNETCDF_VERSION}.tar.gz" \
    && cd "pnetcdf-${PNETCDF_VERSION}" \
    && ./configure --prefix=/usr/local --enable-shared MPICC=mpicc MPICXX=mpicxx MPIF77=mpif77 MPIF90=mpif90 \
    && make -j"$(nproc)" \
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
    && make -j"$(nproc)" \
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
    && make -j"$(nproc)" \
    && make install \
    && rm -rf /tmp/esmf-* /tmp/v${ESMF_VERSION}.tar.gz

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
#   only true for the system gcc; the /container toolchain does not
#   search /usr/local/include by default.
ENV ESMFMKFILE=/usr/local/esmf-${ESMF_VERSION}-debug/lib/esmf.mk \
    CESMDATAROOT=/opt/cesmdata \
    USER=root \
    PKG_CONFIG_PATH=/usr/local/lib/pkgconfig \
    PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1

RUN mkdir -p "${CESMDATAROOT}/inputdata"

LABEL hdf5="${HDF5_VERSION}" \
      netcdf-c="${NETCDF_C_VERSION}" \
      netcdf-fortran="${NETCDF_FORTRAN_VERSION}" \
      pnetcdf="${PNETCDF_VERSION}" \
      esmf="${ESMF_VERSION}" \
      org.opencontainers.image.description="CTSM build/test image for GitHub Actions CI, matching derecho's gnu stack"

CMD ["/bin/bash"]
