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
2. **Satisfy CIME's `container` machine unmodified**: netCDF under `/usr/local`,
   LAPACK/BLAS linkable as `-llapack -lblas`, `ESMFMKFILE`, `USER`, and
   `CESMDATAROOT` preset. (That machine also sets `PNETCDF_PATH=/usr/local`,
   which nothing backs any more and nothing reads — `cime/CIME/Tools/Makefile`
   undefines it outright for `mpilib=mpi-serial`.)
3. **Work in GitHub Actions**: GHA runs steps in *non-login* shells, so
   profile scripts never fire. Everything needed is baked in as Docker `ENV`.

## Software stack vs. derecho gnu

The image carries **one software stack, and it is serial**: no MPICH, no
parallel HDF5/netCDF, no PnetCDF, one `ESMF_COMM=mpiuni` ESMF. A gnu +
mpi-serial run on derecho loads none of those, and carrying a second stack only
created ways to link the wrong half. The toolchain (GCC) lives under `/opt`;
libraries install under `/usr/local`, where the `container` machine hard-codes
`NETCDF_PATH`; ESMF is located via `ESMFMKFILE`.

A `create_test` case here must ask for mpi-serial explicitly — `_Mmpi-serial`
in the test name, or `--mpilib mpi-serial`. The `container` machine's
`<MPILIBS>` still says `mpich`, so a test that names no modifier takes that
default and finds no MPI in this image.

The **derecho gnu** column below is what
`ccs_config/machines/derecho/config_machines.xml` loads at the `ccs_config`
this repo pins (`ccs_config_cesm1.0.88`, i.e. `ncarenv/25.10`). The
**Dockerfile** column is what this directory's `Dockerfile` builds today;
`check-derecho-versions.py` passes against it.

**No published image matches this Dockerfile yet.** The serial-only rebuild
(NEXT_STEPS item 4b) has landed in this directory but has not been built on
Casper or pushed. `cirrus-testing.yml` still pins
`ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261007`, which is the two-stack
image: it was built and validated on Casper on 2026-10-07 against the previous
Dockerfile. A container case build against `:20261007` now fails, deliberately,
because `cime-macros/gnu_container.cmake` names paths only the new image has.
The pin bump and the `_Mmpi-serial` test-name change land together once the
rebuild is published. The older `:20260831` image is the 2026-08-31
`ncarenv/23.09` build (GCC 12.2.0, netCDF-C 4.9.2, ESMF 8.6.0 and so on) and is
kept only as a fallback.

| Component | derecho gnu | Dockerfile | note |
|---|---|---|---|
| GCC | 14.3.0 | 14.3.0 | built from source under `/opt/gcc` |
| HDF5 | 1.14.6 (`hdf5`, pulled in by the serial `netcdf` module via `depends_on`) | 1.14.6, serial, static | not in `config_machines.xml`, so recorded in `derecho-versions.ini` — measured 2026-10-09 against `netcdf/4.9.3` |
| netCDF-C | netcdf/4.9.3 for `mpilib="mpi-serial"` | 4.9.3, serial, static | `NETCDF_C_VERSION` is checked against the serial module, which is the one the container's scope covers |
| netCDF-Fortran | 4.6.2, bundled in the serial `netcdf` module | 4.6.2, serial, static | not in `config_machines.xml`, so recorded in `derecho-versions.ini` — `nf-config --version` on 2026-10-09 against `netcdf/4.9.3` |
| ESMF | `esmf/8.9.1` for `mpilib="mpi-serial"` (that build is `ESMF_COMM=mpiuni`) | 8.9.1, one flavor, `ESMF_COMM=mpiuni` | see "ESMF flavor" below |
| pFUnit | 4.8.0, intel only | 4.8.0, gnu, noMPI/noOpenMP | needed by CTSM's Fortran unit tests; derecho ships no gnu pFUnit, so only the version is matched |
| BLAS/LAPACK | cray-libsci 25.03.0 | reference `lapack`/`blas` (dnf) | libsci is proprietary, so this version is deliberately *not* matched — the dnf packages carry no version ARG at all. `derecho-versions.ini`'s guard records derecho's own 25.03.0 so the stand-in gets re-examined when derecho moves |
| PIO | parallelio-serial/2.6.8 module | 2.6.8, mpi-serial, Fortran on, static | both sides are what a case *links*: derecho through the module's `PIO_LIBDIR`/`PIO_INCDIR`, the container through the same two settings in `cime-macros/gnu_container.cmake` plus `PIO_VERSION_MAJOR` in the image environment. The same install is the mpiuni ESMF's external PIO |
| mpi-serial | mpi-serial/2.5.3 module | 2.5.3 under `/usr/local/mpi-serial` | likewise both sides link the installed one: derecho via the module, the container via `MPI_SERIAL_PATH`. Neither compiles CTSM's `libraries/mpi-serial` submodule |
| distro runtime | `GCC: (SUSE Linux) 13.3.1` | `GCC: (GNU) 11.5.0 (Red Hat)` | the C runtime startup objects an executable picks up from its base system. Deliberately *not* matched and not checked — closing it means rebasing the image on SLES |
| conda, nco | loaded on derecho | not included | not needed to build, nor to run the single-point tests; particular testmods or post-processing may want `nco` |

**Two accepted exceptions**, and only two: BLAS/LAPACK and the distro runtime,
both rows above. The first is *guarded* — `derecho-versions.ini` records
derecho's `cray-libsci` version and the check fails when derecho moves it, so
someone re-examines the stand-in. The second is *documented only*: nothing
short of a SLES base image would close it, so there is no useful check to
write. Everything else in the table is matched version-for-version.

A deliberate *build-flag* choice is not an accepted exception and is not
listed here: `PIO_ENABLE_NETCDF_INTEGRATION` is ON for derecho's
`parallelio-serial` and OFF here, recorded in the `Dockerfile` at the PIO
build, because the flag adds only source files nothing in CTSM or CDEPS calls
and a static link never pulls them into `cesm.exe`.

`git` is also built from source (it is not part of derecho's stack): CIME's
cprnc and PIO builds clone `genf90` via git during the case build. git can't
be `dnf install`ed here because that pulls in `openssh`, whose setuid
`ssh-keysign` fails to `chown` under NCAR's rootless single-UID podman — see
[NEXT_STEPS.md](NEXT_STEPS.md) for the full story.

## ESMF flavor

One ESMF 8.9.1 tree is installed:

- `/usr/local/esmf-8.9.1-mpiuni` (`ESMF_BOPT=O`, `ESMF_COMM=mpiuni`, and
  **`ESMF_PIO=external`** against the mpi-serial PIO) — used by the unit tests
  *and* by every case. It is both the image's default `ESMFMKFILE` and what
  `cime-macros/gnu_container.cmake` names for `MPILIB=mpi-serial`, so you
  should never need to set `ESMFMKFILE`. The external PIO is what lets CDEPS
  read its stream meshes; see "mpi-serial".

It is an `ESMF_COMM=mpiuni` build because `run_tests.py` forces
`mpilib=mpi-serial` and every case here is mpi-serial too, so executables link
with plain `gfortran` and no MPI library — while CTSM's `src/CMakeLists.txt`
calls `find_package(ESMF REQUIRED)` and `link_libraries(esmf)`
unconditionally. Linking an `ESMF_COMM=mpich` ESMF into such an executable
fails with `libesmf.so: undefined reference to symbol 'MPI_Bcast' ... DSO
missing from command line`. The image used to carry two such flavors, for
`mpich` cases it no longer builds.

derecho makes the same choice: whenever `mpilib="mpi-serial"`, its
`config_machines.xml` loads the plain `esmf` module (`esmf/8.9.1`) rather than
the `-debug` variant it loads for MPI builds — even under `DEBUG="TRUE"` — and
that install reports `ESMF_COMM=mpiuni` with `ESMF_BOPT=O` (confirmed against
`esmf/8.9.1` on 2026-10-06, and unchanged from the older `esmf/8.6.0` reading
quoted below). So a `DEBUG=TRUE` case in this container gets this optimized
ESMF; the debug/optimized mismatch is deliberate, and it is derecho's.

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

Several container-specific settings make this work, all in
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
- **`ESMFMKFILE`** and **`NETCDF_PATH`** — select the serial ESMF and the
  serial netCDF when `MPILIB=mpi-serial`.
- **`MPI_SERIAL_PATH`** and **`PIO_LIBDIR`/`PIO_INCDIR`** — point a serial case
  build at the image's mpi-serial and PIO instead of CTSM's submodules. These
  name paths the published `:20261007` does not have, so a container case build
  fails against it until `NEXT_STEPS.md` item 4b rebuilds the image. See that
  file and that item.

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
    --case brazil_test --compset IHistClm60BgcQianRsGs --res 1x1_brazil \
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

The third of those wrappers drives `./run_sys_tests` rather than
`create_test` directly:

```bash
# one test by name
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs

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
| `--machine-name ctsm-ci-container` | picks up `MACHINE_DEFAULTS["ctsm-ci-container"]`: the no-batch launcher, `/scratch` as the testroot base, `/scratch/baselines`, and no account requirement |
| `--wait` (`run_sys_tests`'s own flag, distinct from `create_test`'s -- see "No batch system" below) | `run_sys_tests` otherwise backgrounds `create_test` and returns, and podman would tear the container down, killing it mid-run, while it was still running. With `--wait`, `run_sys_tests` waits for every `create_test` it launched and exits nonzero if any of them failed |
| `--extra-create-test-args "--machine container --compiler gnu"` | `run_sys_tests` never passes `--machine` to `create_test` itself, only `--xml-machine`, so without this CIME tries to auto-detect the machine from the container hostname and fails. A value you pass with `--extra-create-test-args` is **merged** with this, not replaced |

and, only when you asked for a suite with `-s`/`--suite-name`:

| Flag | Why |
|---|---|
| `--xml-machine derecho` | whose testlist to read. The container has no entries of its own in `cime_config/testdefs/testlist_clm.xml`, and replicating derecho is the image's purpose. Override with `$XML_MACHINE` |
| `--suite-compiler gnu` | the only compiler in the image |

**The testroot holds cases, `bld/` and `run/` together**, unlike
`run-case-in-container.sh` and `run-test-in-container.sh`. `run_sys_tests`
names it `tests_<MMDD-HHMMSS><first 2 chars of the machine name>` (so `ct`
here) under `$SCRATCH_DIR`, and passes `--output-root <testroot>` to
`create_test`, so case directories, `bld/` and `run/` all land inside it --
there is no separate `--test-root` / `--output-root` to configure as there is
for the other two wrappers.

**The `/cases` symlink dangles when read from the host.** `run_sys_tests`
also makes a symlink to the testroot in its working directory (`/cases`
inside the container), but that link points at the *container* path
(`/scratch/tests_...`), so from the host it is a dangling symlink. Use the
testroot path the wrapper prints instead.

**Test output does not appear in the job log, and that silence is not a
hang.** `run_sys_tests` launches `create_test` through the no-batch job
launcher, which redirects its stdout/stderr straight to
`<testroot>/STDOUT.<testid>` and `STDERR.<testid>` -- never through this
wrapper's own tee'd log. With `--wait`, nothing streams to the log between
"Running: <create_test ...>" and the final exit code, so a silent PBS log for
a 12-hour suite is expected, not evidence of a hang. Watch those two files
under the testroot for progress instead.

**`cs.status` does not work from the host.** `create_cs_status` bakes
absolute paths into the `cs.status` / `cs.status.fails` scripts it generates
in the testroot -- `/ctsm/cime/CIME/Tools` and the container's own testroot
path -- so those generated scripts run only *inside* the container. From the
host, point `cs.status` at the `TestStatus` files positionally:

```bash
cime/CIME/Tools/cs.status <host testroot>/*/TestStatus
```

**Do not use `--test-root` on its own** -- `cs.status` consults it only inside
`for test_id in test_ids` (`cime/CIME/Tools/cs.status`), and `--test-id`
defaults to empty, so `--test-root <testroot>` with no `--test-id` prints
nothing at all and exits 0. That silence is indistinguishable from a suite
that produced no results. The `--test-id` form works if you prefer it, but
note the id carries a 3-character compiler suffix (`_NUM_COMPILER_CHARS = 3`
in `run_sys_tests.py`), so it is `<MMDD-HHMMSS>ct_gnu`, not `<MMDD-HHMMSS>ct`.

This is the same class of problem as the dangling `/cases` symlink above:
absolute container paths baked in where the host expects host paths.

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
podman pull ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261007
IMAGE_TAG=ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261007 \
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

`run-case-in-container.sh` and `run-test-in-container.sh` add flags only
where you have not passed them yourself (`run-sys-tests-in-container.sh` has
its own, separate table above, in "Running run_sys_tests"):

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

That result was recorded before the ctsm5.4.054 rebase and does not reproduce as
written. The ctsm5.4.054 rebase changed `IHistClm60Bgc` in `cime_config/config_compsets.xml` from `..._MOSART_SGLC_SWAV` to `..._MOSART_DGLC%NOEVOLVE_SWAV`, and `components/cdeps/dglc/cime_config/buildnml:63-65` refuses single-point runs (`single column mode for DGLC is not currently allowed`). `1x1_brazil` is single point, so that test now fails in SHAREDLIB_BUILD. The cdeps guard is not new -- it is in `cdeps1.0.79` too -- so the compset change is what broke it. The equivalent test today is
`SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs`.

**The governing constraint is that exactly ONE MPI implementation may exist in
the executable.** CTSM's mpi-serial build statically links mpi-serial, which
defines every `MPI_*` symbol the executable needs -- measured, `nm -D
--undefined-only cesm.exe` reports **zero** undefined `MPI_*`. Anything else
that drags in a second MPI breaks the run, because only mpi-serial is ever
initialized and the other MPI is called uninitialized:

```
Attempting to use an MPI routine before initializing MPICH
```

Two things in a default build used to do exactly that, when the image still
carried an MPI stack, and `cime-macros/gnu_container.cmake` redirects both on
`MPILIB=mpi-serial`:

| Redirect | What it used to prevent |
|---|---|
| `NETCDF_PATH=/usr/local` | the default netCDF was `--enable-parallel` against MPICH (`ldd libnetcdf.so` showed `libmpi.so.12`); the serial one lived at `/usr/local/serial` and this line was what reached it |
| `ESMFMKFILE` → the mpiuni ESMF | a real-MPI `libesmf.so` carries its own `DT_NEEDED` on `libmpi.so.12` |

Since the image was stripped to one serial stack and that stack moved to
`/usr/local`, both lines name what the image has anyway. They are kept because
they say, at the point of use, which netCDF and which ESMF an `mpi-serial`
executable may link — and because they are what makes the macro work against
any image, not only this one.

`cime/CIME/Tools/Makefile` turns `NETCDF_PATH` into `INC_NETCDF`/`LIB_NETCDF`,
and drops `PNETCDF_PATH` outright for mpi-serial (lines 214-217), which is what
makes the `container` machine's now-unbacked `PNETCDF_PATH=/usr/local`
harmless. The stack is still built **static**, now for a different reason: see
"WHY /usr/local, AND WHY STATIC" in the `Dockerfile`. There is no
`/etc/ld.so.conf.d/ctsm.conf` any more either — nothing shared installs in
`/usr/local/lib`.

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
derecho's own install (the `esmf` module under its `mpi-serial` module
hierarchy, the one gnu mpi-serial builds load; read when that was `esmf/8.6.0`
under `ncarenv/23.09`, since superseded by `esmf/8.9.1` — the arrangement is
what matters here, not the versions, and it was re-confirmed against 8.9.1 on
2026-10-06: still `ESMF_COMM=mpiuni`, `ESMF_BOPT=O`, `-DESMF_PIO=1` and
`-lpioc`, with no bare `-lmpi`):

```
ESMF_COMM:         mpiuni
ESMF_PIO:          external
ESMF_PIO_INCLUDE:  .../parallelio/2.6.2/mpi-serial/2.3.0/gcc/12.2.0/.../include
ESMF_PIO_LIBS:     -lpioc
ESMF_NETCDF:       serial netcdf-c 4.9.2 + netcdf-fortran 4.6.1
```

So the image builds, under `/usr/local`: serial HDF5, netCDF-C, netCDF-Fortran,
and PIO (`PIO_USE_MPISERIAL`, no pnetcdf, no MPI-IO, and `PIO_ENABLE_FORTRAN=ON`
so a case build has a `libpiof` and a `pio.mod` to link); plus mpi-serial under
`/usr/local/mpi-serial`. `ESMF_PIO_LIBS` is `-lpioc -lmpi-serial`, because PIO
compiled against mpi-serial's headers emits real `MPI_*` names that ESMF's
renamed mpiuni stubs do not satisfy. The Dockerfile asserts all of this at build
time, so a PIO-less or MPI-linking serial ESMF fails the build rather than a
later run.

### Linking the image's mpi-serial and PIO, not CTSM's submodules

A case build here links the installs above, the way derecho's `mpi-serial` and
`parallelio-serial` modules make derecho link its own. That takes two halves
and neither works alone:

- `cime-macros/gnu_container.cmake` sets `MPI_SERIAL_PATH`, `PIO_LIBDIR` and
  `PIO_INCDIR`, hardcoded the way `PFUNIT_PATH` and `ESMFMKFILE` are, because
  the container has no modules for `derecho.cmake`'s `$ENV{NCAR_ROOT_MPI_SERIAL}`
  and `$ENV{PIO}` to come from.
- The image exports `PIO`, `PIO_VERSION_MAJOR=2` and
  `PIO_TYPENAME_VALID_VALUES=netcdf` as `ENV`s.
  `cime/CIME/build_scripts/buildlib.pio` reads `PIO_VERSION_MAJOR` from
  `os.environ` — which no macro file reaches — and that is what makes it use an
  installed PIO instead of building one. With `MPI_SERIAL_PATH` set but
  `PIO_VERSION_MAJOR` unset, make skips the mpi-serial submodule build and PIO's
  from-source configure then dies looking for the archive that build would have
  produced.

The acceptance test is on the executable, not on the build succeeding: a freshly
built `cesm.exe` must contain neither `/ctsm/libraries/mpi-serial` nor
`/ctsm/libraries/parallelio` among its compiled-source strings (`strings -a`).
A passing version check without that is the old vacuous behaviour unchanged.

What needs no work at all, and should not be "fixed":

- `<MPILIBS>mpich</MPILIBS>` does not block mpi-serial:
  `Machines.is_valid_MPIlib()` special-cases it. But it *is* the machine
  default, so a test must say `_Mmpi-serial` (or `--mpilib mpi-serial`) to get
  it.
- No `<mpirun mpilib="mpi-serial">` entry is needed; its *absence* is how CIME
  says "launch with no `mpirun`".
- mpi-serial forces `NTASKS=1`, which a 1x1 grid needs, and selects
  `PIO_TYPENAME=netcdf`.

### No batch system

`BATCH_SYSTEM=none` in the machine config means CIME infers no-batch mode and
runs the model in the **foreground**; `case.submit` and `create_test` each
detect it, so neither `create_test`'s own `--no-batch` nor its `--wait` is
needed for `run-case-in-container.sh` or `run-test-in-container.sh`.
`create_test` consequently blocks until the tests finish and reports a real
PASS/FAIL rather than "submitted". This is a different flag from
`run_sys_tests --wait` (see "Running run_sys_tests" above), which
`run-sys-tests-in-container.sh` always injects: `run_sys_tests` backgrounds
`create_test` itself regardless of `BATCH_SYSTEM`, so it needs its own
`--wait` to block and propagate the exit status even though `create_test`
underneath it is still running in no-batch mode.

Do **not** add `--queue`: CIME asserts a queue is never combined with no-batch
mode, so it fails outright.

Because the run is in the foreground, these scripts want a compute allocation.
`run-case-in-container.sh`, `run-test-in-container.sh` and
`run-sys-tests-in-container.sh` all carry PBS headers, so `qsub` any of them
directly -- or run them inside an existing interactive session, which a
single-point `Ld1` case is small enough for.

### When input data is missing

The read-only mount means CIME cannot download anything. Left alone it fails
in a wall of `wget failed` errors (or, where the parent directory does not
exist yet, `OSError: [Errno 30] Read-only file system`), neither of which
obviously means "a file is missing", so each of these three wrappers catches
it and says so.

**Check the dangling-symlink cause first** -- see the mounts table above. It
presents as missing data but is not.

- `run-case-in-container.sh` runs `./check_input_data` (**without**
  `--download`, so it only reports and never writes) before building, and
  stops with a message naming the host path if anything is missing -- before
  wasting a build.
- `run-test-in-container.sh` cannot insert a step, because `create_test`
  drives its own phases, so it recognizes the read-only-filesystem error in
  the log afterwards and explains it.
- `run-sys-tests-in-container.sh` does the same after-the-fact recognition,
  but has to search further: `run_sys_tests` launches `create_test` through a
  launcher that redirects its stdout/stderr to `<testroot>/STDOUT.<testid>` /
  `STDERR.<testid>` rather than through this wrapper's own log (see "Running
  run_sys_tests" above), so it greps those files too whenever it can find the
  testroot.

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
two-stack image measured 3.5 GB, and the build compiles GCC, the netCDF stack
and ESMF from source), and at 4 cores it would take 3-4 h against a 6 h job
limit. The serial-only image is smaller and faster to build, but by how much is
not yet measured. Pushing
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
compiles GCC, HDF5, netCDF-C/Fortran, mpi-serial, PIO, ESMF, git and pFUnit
from source — the two-stack image took about 50 minutes on 16 native cores and
this one takes less. Emulated arm64
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
(`GCC_VERSION`, `HDF5_VERSION`, `NETCDF_C_VERSION`, `NETCDF_FORTRAN_VERSION`,
`ESMF_VERSION`, `PFUNIT_VERSION`, `MPI_SERIAL_VERSION`, `PIO_VERSION`), plus
`GIT_VERSION` and `MAKE_JOBS` (build parallelism). Change one and rebuild.

`check-derecho-versions.py` (run in CI by `derecho-version-check.yml`, or
locally: `python docker/ctsm-ci-derecho-gnu/check-derecho-versions.py`) fails
if these ARGs drift from derecho. Its subject is the software a **gnu +
mpi-serial** run on derecho links, so it reads the `mpilib="mpi-serial"` blocks
of `ccs_config/machines/derecho/config_machines.xml` — `netcdf`, `mpi-serial`,
`parallelio-serial`, `esmf`, alongside the compiler-level `gcc` — and never the
MPI ones. That is now the whole image: every version `ARG` above is covered
except `GIT_VERSION`, since git is not part of derecho's stack at all.

Two things cannot be read from `config_machines.xml` and are recorded by hand
in `derecho-versions.ini`:

- **HDF5 and netCDF-Fortran.** Neither is a standalone module: netCDF-Fortran
  is bundled in `netcdf`, HDF5 is the `hdf5` module `netcdf` pulls in via
  `depends_on`. Verify on derecho with `module show netcdf | grep -i hdf5` and
  `nf-config --version`, and update `measured_against_netcdf` in the same edit
  — it records which `netcdf` the numbers were measured against, the check
  compares it to the live one, and a re-measurement that forgets it fails
  rather than passing.
- **`cray_libsci`.** libsci is proprietary and the container substitutes dnf's
  reference `lapack`/`blas`, which carry no version ARG to compare. The guard
  instead asserts derecho's `cray-libsci` is still the recorded version, so a
  bump there makes someone re-examine the stand-in.

The check also resolves every pinned module `config_machines.xml` names for
this path in derecho's module tree under `/glade` — all but `conda/latest`,
whose "version" is not one, which is reported as skipped rather than passed —
both
`/glade/u/apps/derecho/modules` and the `/glade/u/apps/cesmdev/modules` root
that `cesmdev/1.0` adds, since the `parallelio-serial` `-debug` variant lives
only in the latter. `config_machines.xml` is CTSM's *claim* about derecho,
never verified; if it pins a version derecho has removed, a derecho run dies at
`module load`, and comparing an ARG to that version reports a match that means
nothing. That part needs `/glade`, so it runs on the Cirrus runner (on the
weekly schedule, or by `workflow_dispatch`); the PR gate runs on a hosted
runner and passes `--skip-module-tree`, which prints a `NOT CHECKED` line
rather than anything that could be mistaken for a pass.

pFUnit is read live too, but from a different place: derecho has no pFUnit
module, so the check parses the version out of the `PFUNIT_PATH` that
`ccs_config/machines/derecho/intel_derecho.cmake` sets (e.g.
`pFUnit4.8.0_derecho_Intel2023.2.1_noMPI_noOpenMP`). Bumping `PFUNIT_VERSION`
also means updating the hardcoded install path in
`cime-macros/gnu_container.cmake` — the Dockerfile asserts the two agree at
build time, so a mismatch fails the build rather than a later test run.

**`ESMF_VERSION` has the same coupling.** `gnu_container.cmake`'s
`set(ESMFMKFILE ...)` hardcodes `/usr/local/esmf-8.9.1-mpiuni/lib/esmf.mk`, and
the Dockerfile `COPY`s that file verbatim, so bumping `ESMF_VERSION` means
editing that line in the same change. The Dockerfile's other esmf.mk checks
cannot catch a miss, because they assert against the *templated*
`/usr/local/esmf-${ESMF_VERSION}-mpiuni/…` path, which a bump moves in step —
so a separate assertion reads the path back out of the macro and requires it to
exist and to be the `ESMF_COMM=mpiuni` build. A forgotten edit fails the build
instead of surfacing much later, in `run-unit-tests-in-container.sh` and in
every `mpi-serial` case. The concrete install paths quoted in "ESMF flavor"
and "Baked-in environment" move with the same bump and nothing asserts those,
so edit them in the same change.

## Baked-in environment (why each matters)

| Variable | Value | Why |
|---|---|---|
| `PATH` / `LD_LIBRARY_PATH` | `/usr/local` first, then `/opt/gcc` | non-login GHA shells need the toolchain on `PATH` without profile scripts |
| `ESMFMKFILE` | `/usr/local/esmf-8.9.1-mpiuni/lib/esmf.mk` | CMEPS/CDEPS builds require it; the `container` machine config does not set it. There is one ESMF now, so this and the CIME macro name the same file |
| `PIO` / `PIO_VERSION_MAJOR` / `PIO_TYPENAME_VALID_VALUES` | `/usr/local` / `2` / `netcdf` | what makes CIME link the image's PIO instead of building one. `buildlib.pio` reads `PIO_VERSION_MAJOR` from `os.environ`, which no macro file can reach; an unset `PIO_TYPENAME_VALID_VALUES` is a `TypeError`, not a default |
| `CESMDATAROOT` | `/opt/cesmdata` | `DIN_LOC_ROOT=$CESMDATAROOT/inputdata` must exist even for `--no-run`; it is also the mount point for real inputdata when running (see "Running cases and tests") |
| `USER` | `root` | CIME requires `$USER`; GHA container jobs don't set it |
| `PKG_CONFIG_PATH` / `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS` | `/usr/local/lib/pkgconfig` / `1` | CIME's cprnc locates netCDF via pkg-config with the non-system `/opt` toolchain. Those `.pc` files now describe a *static* stack; see the `Dockerfile` comment for what that means for cprnc |
