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
| PIO | parallelio/2.6.2 module | 2.6.2, mpi-serial only | CIME builds its own PIO from CTSM's pinned ParallelIO submodule for the case; the copy here exists solely as the mpiuni ESMF's external PIO (see "mpi-serial") |
| serial netCDF stack | netcdf/4.9.2 (loaded for mpilib=mpi-serial) | HDF5 1.12.2 + netCDF-C 4.9.2 + netCDF-Fortran 4.6.1 under `/usr/local/serial`, static | mpi-serial builds must not link the parallel, MPICH-linked netCDF |
| mpi-serial | mpi-serial/2.3.0 module | 2.5.4 under `/usr/local/mpi-serial` | only to compile the ESMF-external PIO against; the case build compiles its own from CTSM's submodule |
| conda, nco | loaded on derecho | not included | not needed to build, nor to run the single-point tests; particular testmods or post-processing may want `nco` |

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
- `/usr/local/esmf-8.6.0-mpiuni` (`ESMF_BOPT=O`, `ESMF_COMM=mpiuni`, and
  **`ESMF_PIO=external`** against the mpi-serial PIO) — the **serial** build,
  used by the unit tests *and* by every `mpi-serial` case. Selected
  automatically; you should not need to set `ESMFMKFILE` for it. The external
  PIO is what lets CDEPS read its stream meshes; see "mpi-serial".

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

## Running cases and tests

The image can build **and run** cases, not only compile them. That needs one
thing a build does not: the CESM input data. Three wrappers handle it: the
two below, plus `run-sys-tests-in-container.sh` in the next section.

```bash
# a single-point test, end to end
docker/ctsm-ci-derecho-gnu/run-test-in-container.sh

# or a case you define yourself
docker/ctsm-ci-derecho-gnu/run-case-in-container.sh \
    --case brazil_test --compset IHistClm60Bgc --res 1x1_brazil \
    --mpilib mpi-serial --run-unsupported
```

`--run-unsupported` appears in the case example but not the test one because
CIME skips the scientifically-supported check for test cases
(`cime/CIME/case/case.py`, `if not test and not run_unsupported ...`) -- so
`create_newcase` needs it for an unsupported grid/compset pairing and
`create_test` never does.

Each is a thin wrapper: **every argument passes straight through** to
`create_test` and `create_newcase` respectively. `run-case-in-container.sh`
then runs the `case.setup` -> `preview_namelists` -> `check_input_data` ->
`case.build` -> `case.submit` chain. Shared plumbing lives in
`container-common.sh`, which is sourced on *both* sides of the container
boundary -- the repo is mounted at `/ctsm`, so the in-container half sources
that same file. It therefore defines only functions, with no top-level side
effects.

### Running run_sys_tests

A third wrapper drives `./run_sys_tests` rather than `create_test` directly:

```bash
# one test by name
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc

# derecho's mpi-serial suite, filtered to this image's one compiler
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial
```

Why this rather than `run-test-in-container.sh`: `run_sys_tests` adds the
bookkeeping a real test suite needs that a bare `create_test` wrapper does
not -- a dated testroot, `cs.status` / `cs.status.fails` aggregation, a
recorded `SRCROOT_GIT_STATUS`, baseline compare/generate, and retry. It is
also what a CTSM developer actually types.

Like the other two, every argument passes straight through to
`run_sys_tests`; this wrapper only injects flags you did not pass yourself:

| Flag | Why |
|---|---|
| `--machine-name container` | picks up `MACHINE_DEFAULTS["container"]`: the no-batch launcher, `/scratch` as the testroot base, `/scratch/baselines`, and no account requirement |
| `--wait` | `run_sys_tests` otherwise backgrounds `create_test` and returns, and podman would tear the container down while it was still running |
| `--extra-create-test-args "--machine container --compiler gnu"` | `run_sys_tests` never passes `--machine` to `create_test` itself, only `--xml-machine`, so without this CIME tries to auto-detect the machine from the container hostname and fails. A value you pass with `--extra-create-test-args` is **merged** with this, not replaced |

and, only when you asked for a suite with `-s`/`--suite-name`:

| Flag | Why |
|---|---|
| `--xml-machine derecho` | whose testlist to read. The container has no entries of its own in `cime_config/testdefs/testlist_clm.xml`, and replicating derecho is the image's purpose. Override with `$XML_MACHINE` |
| `--suite-compiler gnu` | the only compiler in the image |

**The testroot holds cases, `bld/` and `run/` together**, unlike the other
two wrappers. `run_sys_tests` names it
`tests_<MMDD-HHMMSS><first 2 chars of the machine name>` (so `co` here) under
`$SCRATCH_DIR`, and passes `--output-root <testroot>` to `create_test`, so
case directories, `bld/` and `run/` all land inside it -- there is no
separate `--test-root` / `--output-root` to configure as there is for the
other wrappers.

**The `/cases` symlink dangles when read from the host.** `run_sys_tests`
also makes a symlink to the testroot in its working directory (`/cases`
inside the container), but that link points at the *container* path
(`/scratch/tests_...`), so from the host it is a dangling symlink. Use the
testroot path the wrapper prints instead.

> **`-s` skips the ctsm_pylib check.** `run_sys_tests` only checks that
> python prerequisites are importable when it has to discover the suite's
> compilers itself. Because this wrapper injects `--suite-compiler gnu`,
> that check is skipped. The practical effect is narrow: of the 20
> `derecho`/`gnu` tests in `aux_clm_mpi_serial`, only
> `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` needs python modules the
> image lacks, and it now fails on its own rather than aborting the other 19
> before any of them run. Naming that test explicitly with `-t` still fails
> up front with `ModuleNotFoundError`, because that path checks
> unconditionally.

### Getting the image onto a compute node

podman's image store here is **node-local** -- `podman info --format
'{{.Store.GraphRoot}}'` reports a path under `/var/tmp` -- so a fresh compute
node starts with an empty store. A `qsub`ed run therefore has to load the
image before the wrapper can use it:

```bash
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_YYYYMMDD.tar
```

The save restores as `localhost/ctsm-ci-derecho-gnu:dev`, which is the
wrappers' default `IMAGE_TAG`, so nothing needs re-tagging. Alternatively pull
the published image and point `IMAGE_TAG` at it:

```bash
podman pull ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260830
IMAGE_TAG=ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260830 \
    docker/ctsm-ci-derecho-gnu/run-test-in-container.sh
```

Budget disk for this. Rootless podman here uses the `vfs` storage driver,
which duplicates every layer rather than sharing them, so the on-disk cost is
several times the image's nominal size.

### Mounts

| Host | Container | Mode | Holds |
|---|---|---|---|
| the repo | `/ctsm` | rw | this checkout |
| `/glade/campaign/cesm/cesmdata/cseg/inputdata` | `/opt/cesmdata/inputdata` | **ro** | `DIN_LOC_ROOT` |
| `$HOME/cases_devcontainer` | `/cases` | rw | case directories |
| `$SCRATCH/cases_devcontainer` | `/scratch` | rw | `bld/`, `run/`, `archive/` |
| `/glade/campaign` | `/glade/campaign` | **ro** | targets of inputdata's symlinks |

**Why that last mount exists.** Much of the CESM inputdata tree is not files
but *absolute symlinks* into sibling campaign collections -- the GSWP3 datm
forcing, for instance, points at
`/glade/campaign/collections/gdex/data/d651077/cesmdata/...`. A bind mount
does not rewrite symlink targets, so with only the `inputdata` mount those
absolute paths do not exist inside the container and every such file looks
**missing**. CIME then tries to download it, the read-only mount refuses, and
the run dies in a wall of `wget failed ... No such file or directory` ending
in `ERROR: Could not find all inputdata on any server`.

Mounting the enclosing tree *at its own path* makes the symlinks resolve. This
is not a corner case: of the 744 distinct files a single-point
`IHistClm60Bgc` test first reported missing, **all 744** were symlinks and
every one pointed under `/glade/campaign`. Override the tree with
`INPUTDATA_SYMLINK_ROOT`, or set it empty to skip the mount.

The inputdata mount lands exactly on `$CESMDATAROOT/inputdata`, which the
Dockerfile creates empty, so `DIN_LOC_ROOT` resolves with **no** `xmlchange`
and no `CESMDATAROOT` override. It is read-only because nothing in a run
should write to the shared tree.

The two persistent directories carry a `_devcontainer` suffix because cases
built here use a different toolchain than a native derecho or casper build;
keeping them apart stops the two from being confused. Container paths are
fixed (`/cases`, `/scratch`) rather than `$HOME`-relative, because `HOME`
differs by context -- `/root` under podman, `/github/home` in GitHub Actions.

`DOUT_S_ROOT` is redirected to `/scratch/archive/$CASE`, since the machine
config's default (`$ENV{HOME}/archive/$CASE`) would be ephemeral container
storage.

Override with `IMAGE_TAG`, `INPUTDATA`, `CASES_DIR`, `SCRATCH_DIR`, `MAKE_J`,
`MOUNT_OPTS` or `DEFAULT_TEST`. **`DRY_RUN=1`** prints the `podman` command
and the assembled CIME argument list without running anything -- handy for
checking the mounts without an allocation.

### Injected defaults

Both wrappers add flags only where you have not passed them yourself:

| Flag | Why |
|---|---|
| `--machine container` | CIME cannot detect the machine from a container hostname |
| `--compiler gnu` | the only compiler in the image. Fills in the compiler for tests named on the command line; a test name that encodes one (`..._gnu`) still wins |
| `--xml-compiler gnu` | *test wrapper only.* The filter for `--xml-*` suite queries -- a **different** knob from `--compiler` -- so a mixed-compiler suite does not try to run its intel entries |
| `--test-root /cases` | *test wrapper only* |
| `--output-root /scratch` | `bld/` and `run/` on the persistent mount |

### mpi-serial

CTSM's single-point tests in `cime_config/testdefs/testlist_clm.xml` are
written `_Mmpi-serial`, and they run here:

```
PASS SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc.container_gnu RUN
```

**The governing constraint is that exactly ONE MPI implementation may exist in
the executable.** CTSM's mpi-serial build statically links mpi-serial, which
defines every `MPI_*` symbol the executable needs -- measured, `nm -D
--undefined-only cesm.exe` reports **zero** undefined `MPI_*`. Anything else
that drags in a second MPI breaks the run, because only mpi-serial is ever
initialized and the other MPI is called uninitialized:

```
Attempting to use an MPI routine before initializing MPICH
```

Two things in a default build would do exactly that, so both are redirected by
`cime-macros/gnu_container.cmake` on `MPILIB=mpi-serial`:

| Redirect | Why |
|---|---|
| `NETCDF_PATH=/usr/local/serial` | the default netCDF is `--enable-parallel` against MPICH (`ldd libnetcdf.so` shows `libmpi.so.12`) |
| `ESMFMKFILE` → the mpiuni ESMF | a real-MPI `libesmf.so` carries its own `DT_NEEDED` on `libmpi.so.12` |

`cime/CIME/Tools/Makefile` turns `NETCDF_PATH` into `INC_NETCDF`/`LIB_NETCDF`,
and already drops `PNETCDF_PATH` for mpi-serial. The serial stack is built
**static** so the loader cannot silently pick the parallel `.so` that
`/etc/ld.so.conf.d/ctsm.conf` puts on the default path.

**The mpiuni ESMF needs PIO, and that is the subtle part.** CDEPS calls
`ESMF_MeshCreateFromFile` to read its stream meshes, which requires ESMF to
have been built with PIO. Without it the run aborts during ATM initialization
with a message that appears **only in `PET0.ESMF_LogFile`**, never in
`cesm.log` -- so it presents as a silent crash:

```
ESMCI_mesh_create_from_file() Library needed by ESMF not present
  - This functionality requires ESMF to be built with the PIO library enabled.
```

ESMF force-disables PIO for `ESMF_COMM=mpiuni`, and `ESMF_PIO=internal` cannot
be forced back on -- ESMF's bundled PIO fails to compile without an `mpi.h`.
The supported route is `ESMF_PIO=external` against a PIO built for mpi-serial,
which is what the image does and what **derecho does too**. Read from
derecho's own install (`esmf/8.6.0` under its `mpi-serial` module hierarchy,
the one gnu mpi-serial builds load):

```
ESMF_COMM:         mpiuni
ESMF_PIO:          external
ESMF_PIO_INCLUDE:  .../parallelio/2.6.2/mpi-serial/2.3.0/gcc/12.2.0/.../include
ESMF_PIO_LIBS:     -lpioc
ESMF_NETCDF:       serial netcdf-c 4.9.2 + netcdf-fortran 4.6.1
```

So the image builds, under `/usr/local/serial`: serial HDF5, netCDF-C,
netCDF-Fortran, and PIO (`PIO_USE_MPISERIAL`, no pnetcdf, no MPI-IO); plus
mpi-serial under `/usr/local/mpi-serial` purely to compile that PIO against.
`ESMF_PIO_LIBS` is `-lpioc -lmpi-serial`, because PIO compiled against
mpi-serial's headers emits real `MPI_*` names that ESMF's renamed mpiuni stubs
do not satisfy. The Dockerfile asserts all of this at build time, so a
PIO-less or MPI-linking serial ESMF fails the build rather than a later run.

What needs no work at all, and should not be "fixed":

- CIME builds mpi-serial from CTSM's own `libraries/mpi-serial` submodule
  (`cime/CIME/config.py`).
- `<MPILIBS>mpich</MPILIBS>` does not block it: `Machines.is_valid_MPIlib()`
  special-cases `mpi-serial`.
- No `<mpirun mpilib="mpi-serial">` entry is needed; its *absence* is how CIME
  says "launch with no `mpirun`".
- It forces `NTASKS=1`, which a 1x1 grid needs, and selects
  `PIO_TYPENAME=netcdf`.

### No batch system

`BATCH_SYSTEM=none` in the machine config means CIME infers no-batch mode and
runs the model in the **foreground**; `case.submit` and `create_test` each
detect it, so neither `--no-batch` nor `--wait` is needed. `create_test`
consequently blocks until the tests finish and reports a real PASS/FAIL rather
than "submitted".

Do **not** add `--queue`: CIME asserts a queue is never combined with no-batch
mode, so it fails outright.

Because the run is in the foreground, these scripts want a compute allocation.
Both carry PBS headers, so `qsub` them directly -- or run them inside an
existing interactive session, which a single-point `Ld1` case is small enough
for.

### When input data is missing

The read-only mount means CIME cannot download anything. Left alone it fails
in a wall of `wget failed` errors (or, where the parent directory does not
exist yet, `OSError: [Errno 30] Read-only file system`), neither of which
obviously means "a file is missing", so both wrappers catch it and say so.

**Check the dangling-symlink cause first** -- see the mounts table above. It
presents as missing data but is not.

- `run-case-in-container.sh` runs `./check_input_data` (**without**
  `--download`, so it only reports and never writes) before building, and
  stops with a message naming the host path if anything is missing -- before
  wasting a build.
- `run-test-in-container.sh` cannot insert a step, because `create_test`
  drives its own phases, so it recognizes the read-only-filesystem error in
  the log afterwards and explains it.

Either way the fix is to run `./check_input_data --download` in the case from
a normal host shell, where the inputdata tree is writable, then re-run.

### If the model segfaults immediately

Look for `Setting resource.RLIMIT_STACK` in the log. The machine config asks
for an unlimited stack (Fortran automatic arrays and compiler-generated array
temporaries live there), and CIME raises only the *soft* limit, leaving the
hard limit alone. That works wherever the hard limit is already unlimited,
which is the normal case on NCAR HPC. If it is not, add
`--ulimit stack=-1:-1` to the `podman run` arguments in
`container-common.sh`.

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

- A [GitHub Personal Access Token
  (Classic)](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#personal-access-tokens-classic)
  with the `write:packages` permission. Existing ones are listed
  [here](https://github.com/settings/tokens); [this
  link](https://github.com/settings/tokens/new?scopes=write:packages) starts
  the setup for a new one. If the ESCOMP org uses SAML SSO, the token must
  also be **SSO-authorized for ESCOMP**, or the push fails with a 403 even
  though the token looks correct.
- Authenticate for the session, keeping the token out of your shell history --
  the same recipe as
  [doc/ctsm-docs_container/README.md](../../doc/ctsm-docs_container/README.md),
  "Pushing to GitHub Container Registry":

  ```bash
  # bash: in this session, commands with a leading space are not saved to history
  export HISTCONTROL=ignoreboth

  # NOTE THE LEADING SPACES, so the secret token is not written to history
     echo YOUR_PERSONAL_ACCESS_TOKEN_CLASSIC | podman login ghcr.io -u YOUR_USERNAME --password-stdin
  ```

  Do not put the token in a file in this repo.

To publish, after `build-on-casper.sh` and all three validation scripts pass.
Podman's storage does not survive the session, so if you have logged out since
building, restore the image first — it saves and restores under the same
`localhost/ctsm-ci-derecho-gnu:dev` name, so nothing needs re-tagging.

**Give the session real memory, or the load is silently OOM-killed.** Podman
here uses the `vfs` storage driver, which duplicates every layer rather than
sharing them, so unpacking the ~3.9 GB archive peaks around 54 GB. An
interactive session started without an explicit `mem=` gets a small default
and `podman load` is SIGKILLed: it prints *nothing*, exits 137, and the next
command fails with the thoroughly misleading

```
Error: localhost/ctsm-ci-derecho-gnu:dev: image not known
```

which points at the tag rather than at memory. (`podman images` may print
`Killed` too.) Get a session with headroom first:

```bash
execcasper -A <PROJECT> -l select=1:ncpus=4:mem=96GB -l walltime=02:00:00
```

Do the load, the `podman login` and the `podman push` all in that same
session -- podman's store is node-local, so an image loaded elsewhere will not
be visible:

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

### amd64 only (unlike ctsm-docs)

[doc/ctsm-docs_container/README.md](../../doc/ctsm-docs_container/README.md)
publishes a multi-architecture manifest via `podman manifest create` plus
`podman build --platform linux/amd64,linux/arm64`. **Do not try that here.**

That works for `ctsm-docs` because it is a ~241 MB image of pip/conda packages
that compiles almost nothing, so QEMU emulation costs little. This image
compiles GCC, MPICH, HDF5, netCDF-C/Fortran, PnetCDF, three ESMF trees, git
and pFUnit from source — about 50 minutes on 16 native cores. Emulated arm64
is 10-20x slower per core, i.e. days.

It is not merely slow, it is unavailable: Casper is x86_64 with no
`qemu-user-static` and no QEMU binfmt handlers registered, and registering
them needs root, which rootless podman here does not have. NCAR has no arm64
hardware to build natively on either.

Note also that an arm64 image would not be replicating derecho, which is
x86_64 — the whole premise of this image. Building one would be a separate
project aimed at local development on Apple Silicon, not a tag on this one.

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
| `CESMDATAROOT` | `/opt/cesmdata` | `DIN_LOC_ROOT=$CESMDATAROOT/inputdata` must exist even for `--no-run`; it is also the mount point for real inputdata when running (see "Running cases and tests") |
| `USER` | `root` | CIME requires `$USER`; GHA container jobs don't set it |
| `PKG_CONFIG_PATH` / `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS` | `/usr/local/lib/pkgconfig` / `1` | CIME's cprnc locates netCDF via pkg-config with the non-system `/opt` toolchain |
