# Next steps: ctsm-ci-gh container

_Last updated: 2026-07-23. The from-scratch image now **builds and validates
end-to-end on Casper**. Remaining work is publishing + repo housekeeping._

## Where things stand

- **Phase 2 (from-scratch image) BUILDS and is VALIDATED on Casper.**
  `Dockerfile.scratch` (FROM `almalinux:9`) builds the full gnu stack — GCC
  12.2.0, MPICH 3.4.3 ch4:ofi, HDF5 1.12.2, netCDF-C 4.9.2, netCDF-Fortran
  4.6.1, PnetCDF 1.12.3, ESMF 8.6.0 (debug + optimized), plus **git built
  from source**. Tagged `localhost/ctsm-ci-gh:dev`. Versions match derecho's
  ncarenv/23.09 gnu stack (see README table).
  - `docker/ctsm-ci-gh/smoke-test.sh` passes (versions, `$ESMFMKFILE`, perl
    `XML::LibXML`, and an MPI + netCDF Fortran hello world that links
    `-lnetcdff -lnetcdf -llapack -lblas` and runs under `mpiexec -n 2`).
  - `create_test --no-run --machine container
    SMS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop.container_gnu.clm-default`
    passes all phases and builds `cesm.exe`.
- **Phase 1 (derived image, `Dockerfile`, FROM the CISL base)** is still
  present and was previously verified; it will be deleted when Phase 2 is
  promoted (see Remaining steps).

## Build & validate on Casper

Helper scripts (run inside a build allocation; the user runs these — the
image lives on node-local podman storage):

- `docker/ctsm-ci-gh/build-on-casper.sh` — wraps `podman build` with the
  node-local `TMPDIR` rootless podman needs; has PBS headers for the
  allocation. **Use `mem` well above 64 GB** (see below).
- `docker/ctsm-ci-gh/smoke-test.sh` — asserts versions + a link/run test.

## Casper build constraints (all handled — don't lose these)

NCAR HPC runs **rootless podman with a single uid mapping** (no
`/etc/subuid`). That shaped several fixes now baked into `Dockerfile.scratch`
/ `build-on-casper.sh` (details in commit messages + in-file comments):

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
  `mem=64GB`; use a larger reservation (e.g. `mem=256GB`). `MAKE_JOBS`
  (default 16) bounds compile parallelism.
- Older fixes still present: ESMF bundled-PIO makefile `$(MAKE)` patch (fixes
  "write jobserver: Bad file descriptor"); `PKG_CONFIG_PATH` +
  `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS` so cprnc finds netCDF via pkg-config.

## Persisting the image (no registry yet)

Node-local storage is wiped when the allocation ends. Keep a known-good copy
on GLADE for reuse without rebuilding (writing a tarball to glade is fine —
unlike the build, it's plain file I/O):

```
podman save -o /glade/work/$USER/ctsm-ci-gh_YYYYMMDD.tar localhost/ctsm-ci-gh:dev
# restore later:  podman load -i /glade/work/$USER/ctsm-ci-gh_YYYYMMDD.tar
```

## Remaining steps

1. **Publish** to a registry (GHCR vs Docker Hub) — deferred; needed before
   CI can pull it. The repo already has a build-and-publish-to-GHCR pattern
   in `.github/workflows/docker-image-build-publish.yml` (for `ctsm-docs`).
2. **Repoint CI** — `.github/workflows/cirrus-testing.yml`
   `simple-build-create_test`: set `container.image` to the published image
   and delete the "Install Perl modules" step and the `USER=`/`CESMDATAROOT=`
   exports (all baked into the image). This is a **simple repoint**: the
   container job keeps running `actions/checkout` + `git-fleximod` in the
   container (git is in the image now), so **no** docker-run/host-checkout
   rework is needed.
3. **Promote** `Dockerfile.scratch` → `Dockerfile`; delete the derived-image
   recipe.
4. **Update `README.md`**: from-scratch build (no CISL base), GCC 12.2.0
   exact match, toolchain under `/opt`, git built from source; drop the
   "shadowed /container libraries" paragraph and the Apple-Silicon `--arch`
   note.
