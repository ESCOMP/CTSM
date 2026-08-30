# Next steps: ctsm-ci-derecho-gnu container

_Last updated: 2026-08-30. The image builds and validates end-to-end on
Casper, including pFUnit and CTSM's Fortran unit tests (55/55), and is
**published and public** at
`ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260830`. `cirrus-testing.yml` is
pointed at it. What remains is a unit-test CI job and the Phase 2 drift
cron._

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
- `.github/workflows/cirrus-testing.yml` `simple-build-create_test` runs in the
  published image, pinned to the dated tag, with the now-baked-in Perl-install
  step and `USER=`/`CESMDATAROOT=` exports removed.
- `README.md` is updated for the from-scratch build.
- **pFUnit + CTSM unit tests work in the container (2026-08-30).** Proven with
  a throwaway derived image (`Dockerfile.pfunit`, now deleted): all 55 CTSM
  unit tests pass under `run_tests.py --machine container`. Three things were
  needed, all now in `Dockerfile`:
  - pFUnit 4.8.0, noMPI/noOpenMP, at `/usr/local/pfunit-4.8.0`.
  - a third ESMF flavor, `ESMF_COMM=mpiuni` / `BOPT=O`, because the
    unit-test link is serial and an mpich ESMF fails it. derecho splits the
    same way; see README "ESMF flavors".
  - `cime-macros/gnu_container.cmake`, dropped into `$HOME/.cime`, carrying
    `PFUNIT_PATH`, `-fallow-argument-mismatch`, and the serial `ESMFMKFILE`.
    See README "Running CTSM's unit tests".

  Helper scripts: `smoke-test-pfunit.sh` (no checkout needed) and
  `run-unit-tests-in-container.sh` (the real thing).

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

The known-good save on disk today is
`/glade/work/$USER/ctsm-ci-gh_20260723.tar` (3.3 GB). It predates the
ctsm-ci-gh -> ctsm-ci-derecho-gnu rename and so restores as
`localhost/ctsm-ci-gh:dev`; re-tag it after loading, or anything referring to
`localhost/ctsm-ci-derecho-gnu:dev` will try to reach a registry named
`localhost` and fail with "pinging container registry localhost ... connection
refused":

```
podman load -i /glade/work/$USER/ctsm-ci-gh_20260723.tar
podman tag localhost/ctsm-ci-gh:dev localhost/ctsm-ci-derecho-gnu:dev
```

That save is now **out of date**: it predates the pFUnit, serial-ESMF and
cime-macros layers. It is still a useful cache for a rebuild (podman can reuse
its layers up to the first change, which is the serial ESMF), but it cannot
run the unit tests. Re-save under the new name after the next rebuild.

## Done 2026-08-30

- **Rebuilt from scratch and re-validated.** ~50 min on Casper.
  `smoke-test.sh`, `smoke-test-pfunit.sh` and
  `run-unit-tests-in-container.sh` (55/55) all pass. Saved to
  `/glade/work/$USER/ctsm-ci-derecho-gnu_20260830.tar`.
- **Published** to `ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu`, tags `20260830`
  and `latest`, package set public (verified pullable anonymously).
  `cirrus-testing.yml` pins the dated tag.

  Publishing is a **manual push from Casper**, not a workflow — see README
  "Publishing to GHCR". Decided **not** to model it on
  `.github/workflows/docker-image-build-publish.yml` (the `ctsm-docs`
  pattern): that builds on `ubuntu-latest`, which has ~14 GB free disk against
  a 3.5 GB image whose build compiles GCC and three ESMF trees from source,
  and 4 cores against a build that takes ~50 min on Casper's 16. A CI publish
  workflow is possible with a disk-reclaim step and amd64-only, but is
  follow-up work, not a blocker. The consequence: **nothing republishes
  automatically** — a Dockerfile change needs a manual rebuild, re-validate,
  push, and a tag bump in `cirrus-testing.yml`.

## Remaining steps

1. **Add a unit-test job to `cirrus-testing.yml`.** Now unblocked. It needs
   the `$HOME/.cime` copy step (GHA overrides `HOME`); see README "Running
   CTSM's unit tests".
2. **Phase 2 drift detection** (see `derecho-versions.ini`): a cron on
   Casper/Derecho reading live derecho versions, opening a GitHub issue on
   drift and emailing on success. Planned as one of the last steps.

## Worth raising upstream

`-fallow-argument-mismatch` never reaches a CTSM unit-test build on any
machine: ccs_config's `gnu.cmake` guards it on
`CMAKE_Fortran_COMPILER_VERSION >= 10`, but `src/CMakeLists.txt` includes
`CIME_initial_setup` (and so the macros) at line 4 and does not call
`project()` until line 10, so CMake has not probed the compiler yet and the
variable is empty. It has gone unnoticed because ccs_config defines
`PFUNIT_PATH` only for the intel builds on derecho, casper and izumi, so the
gnu unit-test path is effectively untravelled. The container works around it
in `gnu_container.cmake`; a real fix belongs in ccs_config or in
`src/biogeochem/ch4varcon.F90` (which calls `mpi_bcast` through an implicit
interface with both `LOGICAL` and `INTEGER`).
