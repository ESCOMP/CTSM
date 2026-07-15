# Next steps: ctsm-ci-gh container (moving work to Casper)

_Last updated: 2026-07-14. Local (Mac/emulated) work stopped deliberately;
continuing on Casper. See "Status" at the bottom for exactly where the
from-scratch build left off._

## Where things stand

- **Phase 1 (derived image) is done and fully verified.** `Dockerfile` (FROM
  `ncarcisl/cisldev-x86_64-almalinux9-gcc12-mpich3:25.08`) built as
  `ctsm-ci-gh:dev`; smoke tests passed and the CI test case
  `SMS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop.container_gnu.clm-default`
  passed CREATE_NEWCASE → MODEL_BUILD (`--no-run`) in a fresh container with
  no manual environment setup.
- **Phase 2 (from-scratch image) is written but not yet verified:**
  `Dockerfile.scratch` (FROM `almalinux:9`; builds GCC 12.2.0, MPICH 3.4.3
  ch4:ofi, HDF5 1.12.2, netCDF-C 4.9.2, netCDF-Fortran 4.6.1, PnetCDF 1.12.3,
  ESMF 8.6.0 debug+optimized). All versions confirmed against derecho's
  ncarenv/23.09 gnu stack (`nf-config` and `module list` checked on derecho
  2026-07-13).
- Two non-obvious fixes are baked into both Dockerfiles — don't lose them:
  1. ESMF bundled-PIO makefile patch (`$(MAKE)` instead of `make`), fixing
     "write jobserver: Bad file descriptor" under parallel builds.
  2. `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig` +
     `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1`, required for CIME's cprnc to find
     netCDF via pkg-config with a non-system toolchain.

## Steps on Casper

1. Get a build allocation: `-l select=1:ncpus=16:mem=64GB` (expect ~45–75 min).
   Podman storage is already node-local (`/var/tmp/<user>_podman/storage`) — fine.
2. Build (native x86_64, no `--arch` flag):
   `podman build -f docker/ctsm-ci-gh/Dockerfile.scratch -t ctsm-ci-gh:dev docker/ctsm-ci-gh/`
3. Smoke tests (non-login shell): `gcc --version` → 12.2.0; `mpichversion` →
   3.4.3; `nc-config --version`/`--prefix` → 4.9.2 at /usr/local;
   `nf-config --version` → 4.6.1; `pnetcdf-config --version` → 1.12.3;
   `test -f "$ESMFMKFILE"`; `perl -MXML::LibXML -e 1`; compile+`mpiexec -n 2`
   an MPI Fortran hello world linking `-lnetcdff -lnetcdf -llapack -lblas`.
4. End-to-end: in a container from the image, clone CTSM (branch
   `cirrus-runner-workflows`), `bin/git-fleximod update -o`, then
   `cime/scripts/create_test --no-run --machine container --output-root /root/t --test-root /root/t SMS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop.container_gnu.clm-default`
   → expect all phases PASS and `bld/cesm.exe` present.
5. Image is on node-local storage — immediately `podman save -o
   /glade/work/$USER/ctsm-ci-gh.tar ctsm-ci-gh:dev` (or push to a registry).

## Then (repo housekeeping, any machine)

6. Promote: replace `Dockerfile` with `Dockerfile.scratch` (rename), delete
   the derived-image recipe.
7. Update `README.md`: from-scratch build (no CISL base), GCC 12.2.0 exact
   match, toolchain under `/opt`, drop the "shadowed /container libraries"
   paragraph and the Apple-Silicon `--arch` note where no longer relevant.
8. Update `.github/workflows/cirrus-testing.yml` `simple-build-create_test`:
   point `container.image` at wherever the image gets published; delete the
   "Install Perl modules" step and the `USER=`/`CESMDATAROOT=` exports (all
   baked into the image now).
9. Decide where to publish the image (Docker Hub vs GHCR) — deferred earlier;
   needed before CI can actually use it.

## Status (local build, stopped 2026-07-14)

The local emulated build of `Dockerfile.scratch` was **stopped on purpose**
partway through — not because of an unresolved failure. Progress before
stopping, in order:

- GCC 12.2.0 layer: **built successfully** (~2 h emulated).
- MPICH 3.4.3 layer: first attempt failed at configure (gfortran 12 needs
  `-fallow-argument-mismatch`); **fixed in Dockerfile.scratch** (FFLAGS/FCFLAGS
  on the configure line) and **rebuilt successfully**.
- HDF5 1.12.2 layer: **built successfully**.
- netCDF-C 4.9.2 layer: failed at configure (`xml2-config` missing on minimal
  AlmaLinux 9); **fixed in Dockerfile.scratch** (`dnf install libxml2-devel`
  at the top of that layer — kept out of the base dnf layer to preserve the
  GCC layer cache) but **this fix has NOT been re-run yet**. The Casper build
  starts effectively here.
- Layers after netCDF-C (netCDF-Fortran, PnetCDF, ESMF ×2, ENV) are identical
  to the Phase-1 recipes that already passed full verification in the derived
  image, so no new surprises are expected — but they have not been built with
  the /opt toolchain yet.

So on Casper: run step 2 (build) and expect it to run start-to-finish;
the local layer cache does not transfer. Then steps 3–9.
