# ctsm-ci-gh: CTSM build/test container for GitHub Actions

This image lets CIME's `container` machine (`--machine container`, compiler
`gnu`; defined in `ccs_config/machines/container/`) compile CTSM cases in
GitHub Actions, e.g. the `simple-build-create_test` job in
[.github/workflows/cirrus-testing.yml](../../.github/workflows/cirrus-testing.yml).
It is **not** the image for running on derecho itself (that will be a separate,
Cray-based image someday — hence the `-gh` suffix).

## Design goals

1. **Match derecho's gnu software stack** (what
   `ccs_config/machines/derecho/config_machines.xml` loads for gnu builds) as
   closely as a public container allows.
2. **Satisfy CIME's `container` machine unmodified**: netCDF/PnetCDF under
   `/usr/local`, LAPACK/BLAS linkable as `-llapack -lblas`, `ESMFMKFILE`,
   `USER`, and `CESMDATAROOT` preset.
3. **Work in GitHub Actions**: GHA runs steps in *non-login* shells, so the
   CISL base image's login-shell environment (`/etc/profile.d/z00-build-env.sh`)
   never fires there. Everything needed is baked in as Docker `ENV` instead.

## Version match vs. derecho gnu

| Component | derecho gnu | this image | note |
|---|---|---|---|
| GCC | 12.2.0 | 12.5.0 | from the CISL base image; exact 12.2.0 planned via a from-scratch build |
| MPI | cray-mpich 8.1.27 | MPICH 3.4.3 | cray-mpich 8.x is MPICH-3.4-ABI-derived; Cray code is proprietary |
| HDF5 | hdf5-mpi/1.12.2 | 1.12.2 (parallel) | confirmed on derecho: loaded by netcdf-mpi/4.9.2 under ncarenv/23.09 |
| netCDF-C | netcdf-mpi/4.9.2 | 4.9.2 | |
| netCDF-Fortran | bundled in netcdf/4.9.2 | 4.6.1 | confirmed on derecho: `nf-config --version` with ncarenv/23.09 + netcdf-mpi/4.9.2 |
| PnetCDF | parallel-netcdf/1.12.3 | 1.12.3 | |
| ESMF | esmf/8.6.0-debug and esmf/8.6.0 | 8.6.0, both flavors | see "ESMF flavors" below |
| BLAS/LAPACK | cray-libsci | reference `lapack`/`blas` (dnf) | libsci is proprietary |
| PIO | parallelio/2.6.2 module | none installed | CIME builds PIO from CTSM's pinned ParallelIO submodule during the case build |
| conda, nco | loaded on derecho | not included | not needed for build-only CI; revisit when tests actually run |

The CISL base image's own newer libraries (netCDF 4.9.3, PnetCDF 1.14.0,
HDF5 1.14.6, PIO 2.6.6, FFTW, heFFTe) remain under `/container/` but are
shadowed by `/usr/local` in `PATH`/`ldconfig` order. They will disappear when
the Dockerfile is converted to a from-scratch (AlmaLinux 9) build.

## ESMF flavors

Two ESMF 8.6.0 trees are installed:

- `/usr/local/esmf-8.6.0-debug` (`ESMF_BOPT=g`) — mirrors the
  `esmf/8.6.0-debug` module derecho loads for gnu `DEBUG=TRUE` builds.
  **This is the default** (`ESMFMKFILE` points here) because the current CI
  test is `SMS_D...` (debug).
- `/usr/local/esmf-8.6.0` (`ESMF_BOPT=O`) — mirrors `esmf/8.6.0`, for
  non-debug tests.

To use the optimized flavor in a workflow step:

```yaml
env:
  ESMFMKFILE: /usr/local/esmf-8.6.0/lib/esmf.mk
```

## Building the image

On an x86_64 Linux box:

```bash
podman build -t ctsm-ci-gh:dev docker/ctsm-ci-gh/
```

On Apple Silicon, add `--arch amd64` (runs under Rosetta emulation — the
ESMF builds make this take a while):

```bash
podman build --arch amd64 -t ctsm-ci-gh:dev docker/ctsm-ci-gh/
```

## Bumping versions

Library versions are `ARG`s at the top of the Dockerfile. If the base image
tag is bumped, also update `GCC_ROOT`/`MPICH_ROOT` to the paths the new base
uses under `/container/`.

## Baked-in environment (why each matters)

| Variable | Value | Why |
|---|---|---|
| `PATH` / `LD_LIBRARY_PATH` | `/usr/local` first, then `/container` gcc+mpich | non-login GHA shells; shadow base-image libs |
| `CC/CXX/FC/F77` | MPICH wrappers | same convention as the base image |
| `ESMFMKFILE` | `/usr/local/esmf-8.6.0-debug/lib/esmf.mk` | CMEPS/CDEPS builds require it; the `container` machine config does not set it |
| `CESMDATAROOT` | `/opt/cesmdata` | `DIN_LOC_ROOT=$CESMDATAROOT/inputdata` must exist even for `--no-run` |
| `USER` | `root` | CIME requires `$USER`; GHA container jobs don't set it (GHA sets `HOME=/github/home`, so `CIME_OUTPUT_ROOT=$HOME/scratch` lands somewhere writable) |
