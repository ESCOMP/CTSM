# ctsm-ci-derecho-gnu: CTSM build/test container for GitHub Actions

This image lets CIME's `container` machine (`--machine container`, compiler
`gnu`; defined in `ccs_config/machines/container/`) compile CTSM cases in
GitHub Actions, e.g. the `simple-build-create_test` job in
[.github/workflows/cirrus-testing.yml](../../.github/workflows/cirrus-testing.yml).
Despite the name, it is **not** meant to run on the derecho HPC machine
itself — it *replicates* derecho's gnu software stack inside a container so
CTSM can be built and tested in CI. An image for running on derecho proper
would be Cray-based and separate.

It is built **from scratch** on AlmaLinux 9 (no vendor base image): every
component below is compiled or installed here so it matches derecho's gnu
stack wherever an open equivalent exists.

## Design goals

1. **Match derecho's gnu software stack** (what
   `ccs_config/machines/derecho/config_machines.xml` loads for gnu builds) as
   closely as a public container allows.
2. **Satisfy CIME's `container` machine unmodified**: netCDF/PnetCDF under
   `/usr/local`, LAPACK/BLAS linkable as `-llapack -lblas`, `ESMFMKFILE`,
   `USER`, and `CESMDATAROOT` preset.
3. **Work in GitHub Actions**: GHA runs steps in *non-login* shells, so
   profile scripts never fire. Everything needed is baked in as Docker `ENV`.

## Software stack vs. derecho gnu

The toolchain (GCC, MPICH) lives under `/opt`; libraries install under
`/usr/local` (where the `container` machine hard-codes `NETCDF_PATH` /
`PNETCDF_PATH`); ESMF is located via `ESMFMKFILE`.

| Component | derecho gnu | this image | note |
|---|---|---|---|
| GCC | 12.2.0 | 12.2.0 | exact match, built from source under `/opt/gcc` |
| MPI | cray-mpich 8.1.27 | MPICH 3.4.3 (ch4:ofi) | cray-mpich 8.x is MPICH-3.4-ABI-derived; Cray code is proprietary |
| HDF5 | hdf5-mpi/1.12.2 | 1.12.2 (parallel) | confirmed on derecho: loaded by netcdf-mpi/4.9.2 under ncarenv/23.09 |
| netCDF-C | netcdf-mpi/4.9.2 | 4.9.2 | |
| netCDF-Fortran | bundled in netcdf/4.9.2 | 4.6.1 | confirmed on derecho: `nf-config --version` with ncarenv/23.09 + netcdf-mpi/4.9.2 |
| PnetCDF | parallel-netcdf/1.12.3 | 1.12.3 | |
| ESMF | esmf/8.6.0-debug, esmf/8.6.0 (the mpi-serial build is `ESMF_COMM=mpiuni`) | 8.6.0, three flavors | see "ESMF flavors" below |
| pFUnit | 4.8.0, intel only | 4.8.0, gnu, noMPI/noOpenMP | needed by CTSM's Fortran unit tests; derecho ships no gnu pFUnit, so only the version is matched |
| BLAS/LAPACK | cray-libsci | reference `lapack`/`blas` (dnf) | libsci is proprietary |
| PIO | parallelio/2.6.2 module | none installed | CIME builds PIO from CTSM's pinned ParallelIO submodule during the case build |
| conda, nco | loaded on derecho | not included | not needed for build-only CI; revisit when tests actually run |

`git` is also built from source (it is not part of derecho's stack): CIME's
cprnc and PIO builds clone `genf90` via git during the case build. git can't
be `dnf install`ed here because that pulls in `openssh`, whose setuid
`ssh-keysign` fails to `chown` under NCAR's rootless single-UID podman — see
[NEXT_STEPS.md](NEXT_STEPS.md) for the full story.

## ESMF flavors

Three ESMF 8.6.0 trees are installed:

- `/usr/local/esmf-8.6.0-debug` (`ESMF_BOPT=g`, `ESMF_COMM=mpich`) — mirrors
  the `esmf/8.6.0-debug` module derecho loads for gnu `DEBUG=TRUE` builds.
  **This is the default** (`ESMFMKFILE` points here) because the current CI
  test is `SMS_D...` (debug).
- `/usr/local/esmf-8.6.0` (`ESMF_BOPT=O`, `ESMF_COMM=mpich`) — mirrors
  `esmf/8.6.0`, for non-debug tests.
- `/usr/local/esmf-8.6.0-mpiuni` (`ESMF_BOPT=O`, `ESMF_COMM=mpiuni`) — a
  **serial** build, used only by the unit tests. Selected automatically; you
  should not need to set `ESMFMKFILE` for it.

The serial flavor exists because `run_tests.py` forces `mpilib=mpi-serial`, so
unit-test executables link with plain `gfortran` and no MPI library — while
CTSM's `src/CMakeLists.txt` calls `find_package(ESMF REQUIRED)` and
`link_libraries(esmf)` unconditionally. Linking an `ESMF_COMM=mpich` ESMF into
such an executable fails with `libesmf.so: undefined reference to symbol
'MPI_Bcast' ... DSO missing from command line`.

derecho has exactly the same split: its `config_machines.xml` loads
`esmf/8.6.0` rather than `esmf/8.6.0-debug` whenever `mpilib="mpi-serial"`,
even under `DEBUG="TRUE"`, and that install reports `ESMF_COMM=mpiuni` with
`ESMF_BOPT=O`. The debug/optimized mismatch here is deliberate, matching
derecho.

To use the optimized MPI flavor in a workflow step:

```yaml
env:
  ESMFMKFILE: /usr/local/esmf-8.6.0/lib/esmf.mk
```

## Building the image

**On NCAR HPC (Casper), rootless podman** — use the wrapper, which sets the
node-local `TMPDIR` that rootless podman/buildah requires and carries PBS
headers for a build allocation:

```bash
docker/ctsm-ci-derecho-gnu/build-on-casper.sh
```

Request generous memory (e.g. `mem=256GB`): the build compiles the whole
toolchain and the final image commit is memory-hungry. See
[NEXT_STEPS.md](NEXT_STEPS.md) for the full set of rootless-podman build
constraints.

**On a host with unrestricted Docker/root** (a workstation or CI runner), a
plain build works:

```bash
podman build -t ctsm-ci-derecho-gnu:dev docker/ctsm-ci-derecho-gnu/   # or: docker build ...
```

Validate a built image with `docker/ctsm-ci-derecho-gnu/smoke-test.sh` (the
core stack) and `smoke-test-pfunit.sh` (pFUnit and the CIME macro plumbing).
Neither needs a CTSM checkout.

## Running CTSM's unit tests

```bash
docker/ctsm-ci-derecho-gnu/run-unit-tests-in-container.sh
```

This mounts the repo at `/ctsm` and runs
`cime/scripts/fortran_unit_testing/run_tests.py` the way
[src/README.unit_testing](../../src/README.unit_testing) describes, plus two
flags the container needs: `--machine container` (CIME cannot detect the
machine from a container's hostname) and `--enable-genf90` (CTSM's `.F90.in`
templates have no committed `.F90` counterparts, and `run_tests.py` points
CMake at the bundled `genf90.pl` only when that flag is set). Expect 55 tests.

Three container-specific settings make this work, all in
`cime-macros/gnu_container.cmake`, which is baked in at
`/opt/ctsm-container/cime-macros/` and copied to `/root/.cime/`. CIME's
`copy_local_macros_to_dir()` copies `$HOME/.cime/*.cmake` into the case, and
`Macros.cmake` includes `${COMPILER}_${MACH}.cmake` last, so that file wins
over everything ccs_config ships without shadowing anything (ccs_config has no
`gnu_container.cmake`). It sets:

- **`PFUNIT_PATH`** — `run_tests.py` scrapes this from the CMake macros and
  reads neither the environment nor `--cmake-args`, so it must come from a
  macro *file*. It names the shared install prefix, not the `PFUNIT-4.8`
  subdirectory; see the comments in that file for why.
- **`-fallow-argument-mismatch`** — ccs_config's `gnu.cmake` adds this, but
  guards it on a compiler version CMake has not probed yet at that point in a
  unit-test build.
- **`ESMFMKFILE`** — selects the serial ESMF when `MPILIB=mpi-serial`.

**In GitHub Actions**, `HOME` is overridden (to `/github/home`), so the
`/root/.cime` copy is not found. A workflow step must place it:

```yaml
- run: |
    mkdir -p "$HOME/.cime"
    cp /opt/ctsm-container/cime-macros/gnu_container.cmake "$HOME/.cime/"
```

## Publishing to GHCR

The image is published to `ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu` by pushing
from Casper, **not** by a GitHub Actions workflow. That is deliberate: the
build needs far more than a standard hosted runner's ~14 GB of free disk (the
image alone is 3.5 GB, and the build compiles GCC and three ESMF trees from
source), and at 4 cores it would take 3-4 h against a 6 h job limit. Pushing
from Casper also publishes exactly the artifact that was validated there.

The consequence is that publishing is a **manual step after any Dockerfile
change** — nothing republishes automatically.

Prerequisites, one time:

- A GitHub personal access token (classic) with `write:packages`. If the
  ESCOMP org uses SAML SSO, the token must also be **SSO-authorized for
  ESCOMP**, or the push fails with a 403 even though the token looks correct.
- `podman login ghcr.io -u <your-github-username>` and paste the token at the
  password prompt. Do not put the token in a file in this repo.

To publish, after `build-on-casper.sh` and all three validation scripts pass.
Podman's storage does not survive the session, so if you have logged out since
building, restore the image first — it saves and restores under the same
`localhost/ctsm-ci-derecho-gnu:dev` name, so nothing needs re-tagging:

```bash
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_YYYYMMDD.tar
```

Then:

```bash
tag=$(date +%Y%m%d)          # e.g. 20260830
img=ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu

podman tag localhost/ctsm-ci-derecho-gnu:dev "${img}:${tag}"
podman tag localhost/ctsm-ci-derecho-gnu:dev "${img}:latest"
podman push "${img}:${tag}"
podman push "${img}:latest"
```

**A new GHCR package is private until you change it.** After the first push,
open the package on GitHub (Packages -> ctsm-ci-derecho-gnu -> Package
settings) and set visibility to public, otherwise `cirrus-testing.yml` cannot
pull it. While it is private, a workflow needs explicit credentials:

```yaml
container:
  image: ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:<tag>
  credentials:
    username: ${{ github.actor }}
    password: ${{ secrets.GITHUB_TOKEN }}
```

That works for runs in the repo itself, but not reliably for pull requests
from forks, where `GITHUB_TOKEN` is read-limited — so public is the simpler
answer. Also link the package to this repository in the same settings page, so
it inherits the repo's visibility and shows up on the repo page.

Workflows should pin the **dated** tag, not `latest`, so a CI run is
reproducible and a republish cannot silently change what CI tested. `latest`
exists for humans.

## Bumping versions

Component versions are `ARG`s near the top of the `Dockerfile`
(`GCC_VERSION`, `MPICH_VERSION`, `HDF5_VERSION`, `NETCDF_C_VERSION`,
`NETCDF_FORTRAN_VERSION`, `PNETCDF_VERSION`, `ESMF_VERSION`,
`PFUNIT_VERSION`), plus
`GIT_VERSION` and `MAKE_JOBS` (build parallelism). Change one and rebuild.

`check-derecho-versions.py` (run in CI by `derecho-version-check.yml`, or
locally: `python docker/ctsm-ci-derecho-gnu/check-derecho-versions.py`) fails
if these ARGs drift from derecho's gnu stack in
`ccs_config/machines/derecho/config_machines.xml`. Most are read live from
that file. HDF5 and netCDF-Fortran (bundled in `netcdf-mpi`, so not standalone
modules there) and the intentional `cray-mpich`→MPICH deviation are recorded
by hand in `derecho-versions.ini`; update that file too when they change
(verify on derecho with `module show netcdf-mpi` / `nf-config --version`).

pFUnit is read live too, but from a different place: derecho has no pFUnit
module, so the check parses the version out of the `PFUNIT_PATH` that
`ccs_config/machines/derecho/intel_derecho.cmake` sets (e.g.
`pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP`). Bumping `PFUNIT_VERSION`
also means updating the hardcoded install path in
`cime-macros/gnu_container.cmake` — the Dockerfile asserts the two agree at
build time, so a mismatch fails the build rather than a later test run.

## Baked-in environment (why each matters)

| Variable | Value | Why |
|---|---|---|
| `PATH` / `LD_LIBRARY_PATH` | `/usr/local` first, then `/opt` mpich + gcc | non-login GHA shells need the toolchain on `PATH` without profile scripts |
| `CC/CXX/FC/F77` | MPICH wrappers (`mpicc`, …) | build everything against MPICH |
| `ESMFMKFILE` | `/usr/local/esmf-8.6.0-debug/lib/esmf.mk` | CMEPS/CDEPS builds require it; the `container` machine config does not set it |
| `CESMDATAROOT` | `/opt/cesmdata` | `DIN_LOC_ROOT=$CESMDATAROOT/inputdata` must exist even for `--no-run` |
| `USER` | `root` | CIME requires `$USER`; GHA container jobs don't set it |
| `PKG_CONFIG_PATH` / `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS` | `/usr/local/lib/pkgconfig` / `1` | CIME's cprnc locates netCDF via pkg-config with the non-system `/opt` toolchain |
