# Next steps: ctsm-ci-derecho-gnu container

_Last updated: 2026-07-23. The image builds and validates end-to-end on
Casper; the CI repoint, Dockerfile promotion, and README are done. **Only
publishing to a registry remains.**_

## Where things stand

- **The from-scratch image builds and is VALIDATED on Casper.** `Dockerfile`
  (FROM `almalinux:9`) builds the full gnu stack — GCC 12.2.0, MPICH 3.4.3
  ch4:ofi, HDF5 1.12.2, netCDF-C 4.9.2, netCDF-Fortran 4.6.1, PnetCDF 1.12.3,
  ESMF 8.6.0 (debug + optimized), plus **git built from source**. Tagged
  `localhost/ctsm-ci-derecho-gnu:dev`. Versions match derecho's ncarenv/23.09 gnu
  stack (see README table).
  - `docker/ctsm-ci-derecho-gnu/smoke-test.sh` passes (versions, `$ESMFMKFILE`, perl
    `XML::LibXML`, and an MPI + netCDF Fortran hello world that links
    `-lnetcdff -lnetcdf -llapack -lblas` and runs under `mpiexec -n 2`).
  - `create_test --no-run --machine container
    SMS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop.container_gnu.clm-default`
    passes all phases and builds `cesm.exe`.
- The old derived-image recipe (FROM the CISL base) has been **deleted**;
  `Dockerfile` is now the from-scratch recipe (formerly `Dockerfile.scratch`).
- `.github/workflows/cirrus-testing.yml` `simple-build-create_test` is
  **repointed** at the image (placeholder `ctsm-ci-derecho-gnu:PUBLISH_TBD` until it is
  published) with the now-baked-in Perl-install step and `USER=`/
  `CESMDATAROOT=` exports removed.
- `README.md` is updated for the from-scratch build.

## Build & validate on Casper

Helper scripts (the user runs these; the image lives on node-local podman
storage):

- `docker/ctsm-ci-derecho-gnu/build-on-casper.sh` — wraps `podman build` (of
  `Dockerfile`) with the node-local `TMPDIR` rootless podman needs; has PBS
  headers for the allocation. **Use `mem` well above 64 GB** (see below).
- `docker/ctsm-ci-derecho-gnu/smoke-test.sh` — asserts versions + a link/run test.

## Casper build constraints (all handled — don't lose these)

NCAR HPC runs **rootless podman with a single uid mapping** (no
`/etc/subuid`). That shaped several fixes now baked into `Dockerfile` /
`build-on-casper.sh` (details in commit messages + in-file comments):

- **Node-local TMPDIR** (`build-on-casper.sh`): buildah's rootfs can't live
  on a parallel FS (glade scratch); set `TMPDIR=/var/tmp/$USER`.
- **No `dnf install git`**: git pulls openssh, whose rpm chowns a setuid file
  (`ssh-keysign`) to a non-root id → fails under single-uid ("cpio: chown
  failed"). git is **built from source** instead — it's needed at runtime
  because CIME's cprnc/PIO clone `genf90` via git during the build. (A
  git-less image was tried and abandoned: cprnc's CMake `ExternalProject_Add`
  requires git.)
- **`TAR_OPTIONS=--no-same-owner`**: upstream tarballs (GCC's, etc.) record
  non-root ownership; tar can't chown to unmapped ids.
- **GCC prerequisites from dnf** (`gmp-devel mpfr-devel libmpc-devel`) instead
  of `contrib/download_prerequisites` — gcc.gnu.org downloads were flaky on
  the compute node. Also set wget `retry_connrefused` for the other tarball
  fetches.
- **Memory**: the final image commit OOM-killed (SIGKILL 137) at
  `mem=64GB`; use a larger reservation (`mem=256GB`). `MAKE_JOBS` (default 16)
  bounds compile parallelism.
- Older fixes still present: ESMF bundled-PIO makefile `$(MAKE)` patch (fixes
  "write jobserver: Bad file descriptor"); `PKG_CONFIG_PATH` +
  `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS` so cprnc finds netCDF via pkg-config.

## Persisting the image (no registry yet)

Node-local storage is wiped when the allocation ends. Keep a known-good copy
on GLADE for reuse without rebuilding (writing a tarball to glade is fine —
unlike the build, it's plain file I/O):

```
podman save -o /glade/work/$USER/ctsm-ci-derecho-gnu_YYYYMMDD.tar localhost/ctsm-ci-derecho-gnu:dev
# restore later:  podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_YYYYMMDD.tar
```

## Remaining step

**Publish** the image to a registry (GHCR vs Docker Hub — undecided), then
replace the `ctsm-ci-derecho-gnu:PUBLISH_TBD` placeholder in `cirrus-testing.yml`
`simple-build-create_test` with the published ref. The repo already has a
build-and-publish-to-GHCR pattern in
`.github/workflows/docker-image-build-publish.yml` (for `ctsm-docs`) to model
it on.

Everything else (CI repoint, Dockerfile promotion, README) is done.
