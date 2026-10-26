# Next steps: ctsm-ci-derecho-gnu container

_Last updated: 2026-10-16. The **published** image, `ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261015`, is item 4b's: one serial stack at `/usr/local`, no MPICH, and a case build that links the image's own mpi-serial and PIO. It builds and validates end-to-end on Casper -- `smoke-test.sh`, pFUnit and CTSM's Fortran unit tests, and single-point **runs** -- and is what `cirrus-testing.yml` pins. All four wrapper scripts have been exercised on Casper against it. (Earlier tags do NOT work with the current `gnu_container.cmake`: `:20261007` predates the serial-only strip, and `:20260830` predates the serial netCDF stack.)_

_**Only mpi-serial works in this image now.** A test name with no `_M` modifier takes the first `<MPILIBS>` entry of `ccs_config`'s container machine, which is `mpich` -- an MPI the image no longer carries. `cirrus-testing.yml` pins `_Mmpi-serial` for that reason; so does anything else that names a test._

_**Two decisions reshape the remaining work.** The container's scope is now
**gnu + mpi-serial only**, so mpich comes out of the image entirely. And the
drift check compares **built executables**, not version metadata -- because
metadata comparison is demonstrably blind to drift that exists today: a
container case build links `MPIserial_2.5.4` from CTSM's submodule while derecho
links its `mpi-serial/2.5.3` module, and `check-derecho-versions.py` reports that
as a match. Reasoning, evidence and the full design live in
`DRIFT_CHECK_DESIGN.md`; the items below assert the decisions without
re-arguing them._

_Unchecked items below are in the order they should be tackled: rebuild the
image to match the rewritten checker (4b), build the drift check (5), then
build the image in CI (8), run tests rather than only build them (9), and take
the workflows off every-push triggers (10)._

## Where things stand

- **`Dockerfile` (FROM `almalinux:9`) now builds ONE serial stack** — GCC 14.3.0, HDF5 1.14.6, netCDF-C 4.9.3, netCDF-Fortran 4.6.2 (all serial and static), mpi-serial 2.5.3, PIO 2.6.8 with the Fortran interface, ESMF 8.9.1 `ESMF_COMM=mpiuni` only, pFUnit 4.8.0, plus **git built from source**. These match derecho's `ncarenv/25.10` gnu + mpi-serial stack — see the README table. **Not yet built**: item 4b did this on the repo side only.
  - What item 4b took out, relative to the published `:20261007`: MPICH 3.4.3 ch4:ofi, PnetCDF 1.14.1, the MPI-linked netCDF-C and HDF5, and the two `ESMF_COMM=mpich` ESMF flavors.
  - `docker/ctsm-ci-derecho-gnu/smoke-test.sh` was rewritten with the image and has not been run against one. It now checks versions from labels, both `ESMFMKFILE` routes, the PIO install and its three `ENV`s, the stack-facts file, perl `XML::LibXML`, and a **serial** netCDF Fortran hello world built with `gfortran` — no `mpiexec`, no `use mpi`.
  - The most recent validated build is item 4b's, 2026-10-15, still unpublished: all four wrapper scripts and `smoke-test.sh` pass, and `cesm.exe` meets the acceptance test below. The last validated build of the *published* image was 2026-10-07: `create_test --no-run --machine container SMS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop.container_gnu.clm-default` passed all phases and built `cesm.exe`.
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
- **Node-local `XDG_RUNTIME_DIR`** — not yet done, and it bites *interactive*
  podman rather than the build script. With it unset and `/run/user/$UID`
  absent, podman falls back to `$TMPDIR`, which in an interactive session is on
  glade: shared across nodes and never cleared. A `pause.pid` left there names a
  PID that may belong to an unrelated process on the next node, and podman then
  fails every command — including `system migrate` and `system reset`, which
  die before dispatching because the rootless re-exec happens first. Deleting
  the file clears it. Set `XDG_RUNTIME_DIR=/var/tmp/${USER}_podman/run` (mode
  700) in the wrapper scripts, and in `.bashrc` guarded on unset so a session
  manager's value is never clobbered. `build-on-casper.sh` is immune only
  because it forces a node-local `TMPDIR` above.
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
  a 3.1 GB image whose build compiles GCC and ESMF from source, and 4 cores
  against a build that takes ~50 min on Casper's 16. A CI publish
  workflow is possible with a disk-reclaim step and amd64-only, but is
  follow-up work, not a blocker. The consequence: **nothing republishes
  automatically** — a Dockerfile change needs a manual rebuild, re-validate,
  push, and a tag bump in `cirrus-testing.yml`. **Revisited 2026-10-06:** this
  is now tracked as "Remaining steps" item 8, which keeps the disk and core
  constraints noted here but drops the amd64-only conclusion -- GitHub's native
  arm64 runners remove the emulation objection that stood behind it.

## Added 2026-08-31: run wrappers (VALIDATED on Casper)

Scripts that make the image do real runs, not just builds. **A single-point
CTSM case now runs to completion in the container**, reading inputdata from a
read-only campaign mount:

```
PASS SMS_D_P1_Ld1.1x1_brazil.IHistClm60Bgc.container_gnu RUN   (161.6 s)
```

- `container-common.sh` -- shared plumbing, sourced on **both** sides of the
  container boundary (the repo is mounted at `/ctsm`), hence functions only
  and no top-level side effects.
- `run-case-in-container.sh` -- thin wrapper over `create_newcase`, then
  `case.setup` -> `preview_namelists` -> `check_input_data` -> `case.build` ->
  `case.submit`. Verified on Casper: exit 0, having written
  `*.clm2.r.1850-01-06-00000.nc` to the archive mount.
- `run-test-in-container.sh` -- thin wrapper over `create_test` (no
  `--no-run`). Verified PASS.

See README "Running cases and tests". Findings worth not re-deriving:

- **inputdata is full of absolute symlinks, and that breaks naively.** Much of
  the tree points into sibling campaign collections (e.g. GSWP3 datm forcing
  -> `/glade/campaign/collections/gdex/data/...`). A bind mount does not
  rewrite symlink targets, so with only the `inputdata` mount they dangle and
  every such file looks *missing*; CIME then tries to download it and dies in
  a wall of `wget failed`. Of the 744 distinct files the first run reported
  missing, **all 744** existed on the host, **all 744** were symlinks, and
  every one pointed under `/glade/campaign`. Hence the second read-only mount
  of `/glade/campaign` at its own path (`INPUTDATA_SYMLINK_ROOT`).
- **podman's image store is node-local** (`/var/tmp/...`), so a `qsub`ed run
  must `podman load` the image first; the tarball restores as
  `localhost/ctsm-ci-derecho-gnu:dev`, already the wrappers' default. Load
  takes ~75 s.
- **No batch flags are needed** -- `BATCH_SYSTEM=none` makes CIME infer
  no-batch mode, which also makes `create_test` block for a real PASS/FAIL, so
  `--wait` is redundant and `--queue` is a hard error.
- **mpi-serial needs no `ccs_config` change** -- CIME builds it from CTSM's
  `libraries/mpi-serial` submodule, `is_valid_MPIlib()` special-cases it, and
  a missing `<mpirun mpilib="mpi-serial">` entry *is* how CIME says "no
  launcher". Do not "fix" the machine config.
- `--compiler` and `--xml-compiler` are different knobs; suite queries filter
  on the latter.
- **`RLIMIT_STACK` was a non-issue** -- no `--ulimit` needed.
- Peak memory for load+build+run was ~54-56 GB, so request well above 64 GB.

## Resolved 2026-08-31: mpi-serial runs work

`_Mmpi-serial` tests -- how CTSM writes all its single-point tests -- now run
in the container:

```
PASS SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc.container_gnu RUN
```

(Recorded before the ctsm5.4.054 rebase. That alias now resolves to
`DGLC%NOEVOLVE`, which cdeps refuses to run single-point, so the test above no
longer builds; `SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs` is the
equivalent today. See README "mpi-serial".)

with no regression: 55/55 unit tests and the mpich `_P1` run still pass.

**The one rule that explains every failure along the way: exactly ONE MPI
implementation may exist in the executable.** CTSM's mpi-serial build
statically defines every `MPI_*` symbol (`nm -D --undefined-only cesm.exe`
reports zero undefined `MPI_*`), so anything dragging in a second MPI leaves
that MPI uninitialized -> `Attempting to use an MPI routine before
initializing MPICH`.

Three attempts were eliminated by evidence, and should not be revisited:

| Attempt | Killed by |
|---|---|
| real-MPI ESMF for cases | `libesmf.so` has its own `DT_NEEDED` on `libmpi.so.12`; ESMF's calls go to MPICH while CTSM's go to mpi-serial. Also unfixable by exporting symbols -- ESMF is compiled against MPICH headers |
| mpiuni ESMF as shipped | no PIO, and CDEPS needs `ESMF_MeshCreateFromFile` |
| `ESMF_PIO=internal` + mpiuni | ESMF's bundled PIO: `pio.h:16: fatal error: mpi.h`. Matches ESMF's own note that internal PIO "does not support mpiuni mode" |

**The working recipe, which is what derecho does.** derecho's install is
readable from Casper, and its `esmf/8.6.0` under the `mpi-serial` module
hierarchy (spack hash `him6`) reports:

```
ESMF_COMM: mpiuni     ESMF_PIO: external
ESMF_PIO_INCLUDE: .../parallelio/2.6.2/mpi-serial/2.3.0/gcc/12.2.0/.../include
ESMF_PIO_LIBS: -lpioc          ESMF_NETCDF: serial netcdf-c + netcdf-fortran
```

Now in `Dockerfile`: serial HDF5 + netCDF-C + netCDF-Fortran + PIO 2.6.2
(`PIO_USE_MPISERIAL`, no pnetcdf/MPI-IO) under `/usr/local/serial`, static;
mpi-serial 2.5.4 under `/usr/local/mpi-serial` purely to compile that PIO
against; and the mpiuni ESMF rebuilt with `ESMF_PIO=external` and
`ESMF_PIO_LIBS="-lpioc -lmpi-serial"`. In `cime-macros/gnu_container.cmake`,
`MPILIB=mpi-serial` selects that ESMF and `NETCDF_PATH=/usr/local/serial`.

Versions and paths in that paragraph are as of 2026-08-31 and have since been superseded twice -- by item 6 (the version bumps) and by item 4b (the serial stack moved to `/usr/local` and PIO gained its Fortran interface). The *arrangement* it describes is what is still worth reading.

Gotchas worth not rediscovering:

- **ESMF errors go to `PET0.ESMF_LogFile`, not `cesm.log`.** The PIO failure
  presented as a totally silent crash for hours because nothing looked there.
  Check it first for any run that dies during initialization.
- mpi-serial's `make install` is broken upstream (`$(INSTALL)` and
  `$(MKINSTALLDIRS)` are never defined), so the artifacts are copied by hand,
  as CIME's own `buildlib.mpi-serial` does.
- PIO compiled against mpi-serial emits real `MPI_*` symbol names, which
  ESMF's renamed mpiuni stubs do not satisfy -- hence `-lmpi-serial` in
  `ESMF_PIO_LIBS`. Only the archive is copied into the serial stack's `lib`
  (`/usr/local/serial/lib` then, `/usr/local/lib` since item 4b); the headers
  stay out, in `/usr/local/mpi-serial/include`, so they cannot shadow anything
  on the case build's include path.
- The image bakes a snapshot of `gnu_container.cmake`, which silently shadowed
  edits during debugging. The wrapper scripts now prefer the mounted repo's
  copy.

`Dockerfile.serial-netcdf` was the throwaway probe used to prove all this. It
has been folded into `Dockerfile` and deleted, exactly as `Dockerfile.pfunit`
was before it. Rebuilt and re-validated from scratch on 2026-08-31; saved to
`/glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar`, **published** as
`ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260831`, and pinned in
`cirrus-testing.yml`.

One publishing gotcha, since it presents as something else entirely: if
`podman load` is OOM-killed, it exits **137** having printed nothing, and the
next `podman tag` fails with the misleading `image not known`. `vfs` storage
duplicates every layer, so a 3.9 GB archive peaks around 54 GB (measured on the
two-stack image; 4b's is smaller, so treat this as an upper bound) -- load it from
an allocation with real memory, e.g.
`execcasper -A <PROJECT> -l select=1:ncpus=4:mem=96GB -l walltime=02:00:00`.

## Added 2026-08-31: run_sys_tests wrapper (written, not yet validated on Casper)

The fourth wrapper script in this directory (after `run-case-in-container.sh`,
`run-test-in-container.sh` and `run-unit-tests-in-container.sh`),
`run-sys-tests-in-container.sh`, drives `./run_sys_tests` rather than
`create_test` directly -- the bookkeeping a real test suite needs that the
bare `create_test` wrapper does not provide: a dated testroot,
`cs.status`/`cs.status.fails` aggregation, a recorded `SRCROOT_GIT_STATUS`,
baseline compare/generate, and retry. See README "Running run_sys_tests".

Four small CTSM-side changes made this possible, none touching `ccs_config`
or `testlist_clm.xml`:

- **`MACHINE_DEFAULTS["ctsm-ci-container"]`**
  (`python/ctsm/machine_defaults.py`, commit `75727441f`). Without it,
  `run_sys_tests` falls through to its unknown-machine branch: `scratch_dir`
  and `baseline_dir` come back `None`, making `--testroot-base` mandatory and
  `--compare`/`--generate` unusable. Named `ctsm-ci-container` rather than
  the shorter `container` because `ccs_config/machines/container/` already
  defines a generic CIME machine named `container` (`MACH="container"` in its
  `config_machines.xml`); a `MACHINE_DEFAULTS["container"]` key would collide
  with that unrelated machine. `ctsm-ci-container` now exists in
  `MACHINE_DEFAULTS`, giving the no-batch launcher, `/scratch` as the
  testroot base, `/scratch/baselines`, and no account requirement. One
  consequence: `run_sys_tests` derives the testid prefix as the first two
  characters of the machine name, so the testroot is
  `tests_<MMDD-HHMMSS>ct`, not `...co`.
- **`run_sys_tests --xml-machine`** (`python/ctsm/run_sys_tests.py`, commit
  `97ff51c85`). Separates "what machine am I" from "whose testlist do I
  read" -- previously `--xml-machine` was hardcoded to the machine name, so a
  machine with no `testlist_clm.xml` entries of its own (the container)
  could not use suite mode (`-s`) at all. Defaults to the machine name, so no
  existing invocation changes. Now also rejected up front when passed
  without `--suite-name`, since it is meaningless on the `-t`/`-f` paths.
- **`run_sys_tests --wait`** (`python/ctsm/run_sys_tests.py` +
  `joblauncher/`, commit `0be3c0f8e`, corrected post-review). The no-batch
  launcher `Popen`s `create_test` and returns without waiting for it -- right
  on a login node, where the job should outlive the shell, but fatal in a
  container: the wrapper's shell exits, podman tears down the PID namespace,
  and kills any `create_test` still running underneath it. `--wait` now waits
  for every `create_test` process `run_sys_tests` launched -- including when
  a suite spans multiple compilers and launches more than one -- and exits
  nonzero if any of them failed. Also now rejected up front on a launcher
  that cannot wait (e.g. qsub), rather than failing only after `create_test`
  has already been dispatched. Opt-in, so no existing invocation changes.
- **The no-batch job launcher's wait method** (`joblauncher/`, same commit as
  `--wait`, corrected post-review) now waits for every process it launched,
  not just the most recently launched one, and returns nonzero if any of
  them failed. The base launcher class's own version of this method already
  existed (since commit `0be3c0f8e`) purely to raise for launcher types that
  cannot wait at all (e.g. qsub); that behavior is unchanged here, but the
  base class also gained a separate predicate that `--wait`'s new pre-flight
  check (above) queries directly, so an incompatible launcher is now rejected
  before `create_test` is dispatched rather than only discovered afterward
  via the raise.

`--xml-machine` and `--wait` are both upstreamable on their own merits, not
container-specific plumbing masquerading as general features: each defaults
to prior behavior for every existing machine and invocation.

**No image rebuild was needed for any of this.** The wrapper and all four
CTSM-side changes run entirely from the repo mounted read-write at `/ctsm`
(see `container-common.sh`); nothing here touches the `Dockerfile`.

**Not yet run on Casper.** The wrapper is written, syntax-checked, and its
argument-injection logic -- `--machine-name ctsm-ci-container`, `--wait`, the
`--extra-create-test-args` merge, and the suite-only `--xml-machine`/
`--suite-compiler` injections -- is verified on the login node under
`DRY_RUN=1`, which makes `ctsm_podman_run` print the assembled `podman run`
command and `run_sys_tests` argument line instead of running either. It has
not yet been run against the actual container on a Casper compute node; see
"Remaining steps" below.

## Remaining steps

_Each item below includes updating `README.md` wherever it changes what the
README describes -- the component table, the deviation list, the install
paths. That is part of the item, not a follow-up. Item 12 is a different
kind of documentation and is tracked on its own._

1. ✅ **Add a unit-test job to `cirrus-testing.yml`.** Now unblocked. It needs
   the `$HOME/.cime` copy step (GHA overrides `HOME`); see README "Running
   CTSM's unit tests".
2. ✅ **Wire runs into CI.** `simple-build-create_test` currently runs on
   `ubuntu-latest`, which has no `/glade` at all, so a run job has to move to
   `runs-on: gha-runner-ctsm` and bind-mount inputdata into the container job.
   Whether container jobs on that runner see `/glade` is **unknown** -- the
   existing `list-glade-cesm-input` job proves only that the *host* does. Worth
   a probe workflow in the style of `probe-derecho-modules.yml`. Now that
   `run_sys_tests --wait` blocks on every launched test and exits nonzero if
   any of them failed, CI has the exit code it needs to fail the job on a
   test failure -- that piece is no longer a gap.
3. ✅ **Validate `run-sys-tests-in-container.sh` on Casper.** See VALIDATION_2026-09.md.
4. **Make the checker and the image agree on the right subject.** Two steps, in
   order: fix what the check asks, then rebuild the image to satisfy it.

   **4a ✅ -- `check-derecho-versions.py` rewritten.** It now compares the
   Dockerfile ARGs against the modules a gnu + mpi-serial run on derecho loads
   -- the `mpilib="mpi-serial"` blocks of `config_machines.xml`, never the MPI
   ones. The MPI-flavor checks, the `cray-mpich` guard, the `<MPILIBS>`
   assertion and the two-stack twin machinery are gone. `cray-libsci` is now
   guarded; since its stand-in (dnf `lapack`/`blas`) carries no version ARG,
   `deviation` mode no longer requires one. `MPICH_VERSION` and
   `PNETCDF_VERSION` stayed in the `Dockerfile`, simply unchecked, until 4b
   deleted them -- which it has. Reasoning and evidence:
   `DRIFT_CHECK_DESIGN.md`.

   Three things a later reader should not have to re-derive:

   - **The snapshot numbers did not move; only the provenance did.** HDF5
     1.14.6 and netCDF-Fortran 4.6.2 are what the serial `netcdf/4.9.3` gives,
     the same values `netcdf-mpi` gave. The ini key is now
     `measured_against_netcdf`, re-measured on Casper 2026-10-09.
   - **derecho has two module roots and the existence check needs both.**
     `/glade/u/apps/derecho/modules` and `/glade/u/apps/cesmdev/modules`, the
     latter contributed by `cesmdev/1.0`. `parallelio-serial/2.6.8-debug`,
     which `config_machines.xml` loads for `DEBUG="TRUE"`, exists *only* under
     the CSEG root, so a check that walked the system root alone would have
     failed on its first run.
   - **Three search tiers, and no others.** `<system>/environment/<name>/` (for
     `ncarenv` and `cesmdev` themselves), then `<root>/<env>/Core/` and
     `<root>/<env>/gcc/<gccver>/` under each root -- exactly what ncarenv's and
     gcc's modulefiles put on `MODULEPATH`. A module found only under
     `<root>/<env>/cray-mpich/<ver>/gcc/<gccver>/` must **not** count: that
     directory is never on a serial run's `MODULEPATH`, and accepting it would
     hide the drift the check exists for.

   The existence check needs `/glade`, so `derecho-version-check.yml` now picks
   its runner per event: `gha-runner-ctsm` for `schedule` (weekly) and
   `workflow_dispatch`, `ubuntu-latest` plus `--skip-module-tree` for
   `pull_request` and `push`. `test_check_derecho_versions.py` (stdlib
   `unittest`, temp-directory fixtures) mutation-tests each check and runs as
   its own CI step before the check itself.

   **4b ✅ -- the image is serial-only and links its own mpi-serial and PIO.** Built, validated on Casper, published as `:20261015` and pinned in `cirrus-testing.yml`, all 2026-10-15. Everything in this item needed one Casper rebuild, so the repo side landed as one change. The only bullet left below is the upstream `ccs_config` PR, which is deliberately separate work shared with CESM. Reasoning and evidence: `DRIFT_CHECK_DESIGN.md`.

   What landed:

   - `Dockerfile`: MPICH, PnetCDF, the MPI-linked netCDF-C and HDF5 builds and both `ESMF_COMM=mpich` ESMF flavors are gone, with `MPICH_VERSION` and `PNETCDF_VERSION`. The serial stack is the only stack and moved from `/usr/local/serial` to `/usr/local`, where `ccs_config/machines/container/container.cmake:4` hardcodes `NETCDF_PATH`. `PNETCDF_PATH=/usr/local` survives in that machine config and is harmless: `cime/CIME/Tools/Makefile:214-217` undefines it outright for `mpilib=mpi-serial`. PIO is rebuilt with `-DPIO_ENABLE_FORTRAN=ON`, and the image exports `PIO`, `PIO_VERSION_MAJOR=2` and `PIO_TYPENAME_VALID_VALUES=netcdf`.
   - `cime-macros/gnu_container.cmake`: `NETCDF_PATH` moved to `/usr/local`. `ESMFMKFILE`, `MPI_SERIAL_PATH`, `PIO_LIBDIR` and `PIO_INCDIR` were already written against these paths and did not change.
   - The `Dockerfile` grew a `MPI_SERIAL_PATH` macro-path assertion beside the `PFUNIT_PATH` and `ESMFMKFILE` ones, and a PIO install assertion on `libpiof.a`/`libpioc.a`/`pio.mod`.
   - `smoke-test.sh` lost the MPI hello-world and the mpich/pnetcdf checks and gained checks on the PIO install, the three PIO `ENV`s and the stack-facts file.
   - `README.md`'s accepted-exception list is two entries: `cray-libsci` (guarded) and the distro runtime (documented only).
   - `Dockerfile.pio-fortran`, the throwaway probe, is deleted; its PIO layer, flag list, assertions and `ENV`s are folded into the `Dockerfile`.

   Four things a later reader should not have to re-derive:

   - **The `git` block moved above the PIO block, deliberately.** With `-DPIO_ENABLE_FORTRAN=ON`, ParallelIO's `CMakeLists.txt:137-158` shells out to `git clone CMake_Fortran_utils` unless `USER_CMAKE_MODULE_PATH` is defined, and the whole block sits inside `if (PIO_ENABLE_FORTRAN)` -- which is why the old `FORTRAN=OFF` build could sit ahead of a git that is built from source. The clone is done explicitly in the PIO layer and the path passed in, so nothing is fetched behind cmake's back. Cost of the move: editing the git layer now rebuilds PIO, ESMF and pFUnit too.
   - **`/etc/ld.so.conf.d/ctsm.conf` is gone, not rehomed.** It was written by the MPICH block and listed `/usr/local/lib` plus the mpich lib dir. After the move nothing shared installs in `/usr/local/lib`: HDF5, netCDF-C and netCDF-Fortran are `--disable-shared --enable-static`, PIO is `BUILD_SHARED_LIBS=OFF`, mpi-serial is a hand-copied `.a`, git installs to `bin/`, `libexec/` and `share/`, and pFUnit and ESMF have their own prefixes, which that file never listed. `ENV LD_LIBRARY_PATH` still carries `/usr/local/lib` if anything shared ever does appear.
   - **The build-time stack facts are a FILE, not `LABEL`s.** This item originally asked for `LABEL`s recording `nf-config --version` and `nc-config --has-stdfilters`. `LABEL` takes a literal string and cannot run a command, and a `LABEL` fed from a build `ARG` would only repeat what the `Dockerfile` says rather than what the build produced -- which is the entire point of those two values. They are written to `/opt/ctsm-container/stack-facts.txt` as `key=value`, both required non-empty, and item 5 reads them with `podman run ... cat`.
   - **`/usr/local/mpi-serial` has no `mpi.mod`, and never will.** mpi-serial's `fort.F90` is bare subroutines with `include "mpif.h"`, no `module mpi`, so the `Dockerfile`'s best-effort `cp mpi.mod` always no-ops. PIO's `FindMPISERIAL` asks for `INCLUDE_NAMES mpi.mod mpif.h` and `cmake/LibFind.cmake:239-249` breaks on the first hit, so `mpif.h` alone satisfies the Fortran component. Not a blocker; do not go looking for the missing `.mod`.

   What is left, in order:

   - ✅ **Built and validated on Casper, 2026-10-15.** All four wrapper scripts and `smoke-test.sh` pass against the rebuilt image. Three defects surfaced during that rebuild and are fixed in the `Dockerfile` and the wrappers: the `.pc` files describing a shared link for a static stack, `CC`/`CXX`/`FC` going missing with the MPI stack so cmake fell back to the distro compiler, and the `run_sys_tests` wrapper injecting neither half of the required compare/generate argument pairs. Each is written up in its own bullet or at its code site. ✅ **Published as `:20261015` and pinned in `cirrus-testing.yml`, with `_Mmpi-serial` added to that job's test name in the same commit** -- the two are one switch, since the pin alone would ask the new image for MPICH and the modifier alone would run against an image whose paths `gnu_container.cmake` no longer matches. **Item 4b is complete.**
   - **Deleting `ENV CC/CXX/FC` regressed cprnc onto the distro compiler, and the runtime `ENV` now sets them back.** Measured 2026-10-15, on the rebuild after the `.pc` fix below: cprnc linked the static HDF5 correctly and then died on undefined `__truncsfhf2`/`__truncdfhf2`/`__truncxfhf2` from `libhdf5.a(H5Tconv_float.o)` -- libgcc's `_Float16` helpers, which GCC 11.5.0 does not carry because x86-64 `_Float16` arrived in GCC 12. `CMakeCache.txt` recorded `CMAKE_Fortran_COMPILER:FILEPATH=/usr/bin/f95`. cprnc is built by cmake (`buildlib.cprnc`), the generated cmake line passes no `-DCMAKE_*_COMPILER`, and `env_mach_specific.sh` sets no compilers -- so cmake used its built-in name search for `cc` and `f95`, neither of which exists in `/opt/gcc/<ver>/bin` (it ships `gcc` and `gfortran`), and PATH order cannot help because cmake searches by name before path. Before 4b the image set `CC=mpicc`/`FC=mpifort`, which cmake reads as `$ENV{CC}`/`$ENV{FC}`; those wrapped the right GCC, so cprnc was correct by accident of the MPI stack. The `.pc` failure masked this one until it was fixed. Absolute paths now, with a build-time assertion that each is executable, since a `GCC_VERSION` bump is the one edit that could leave them dangling silently.
   - **`cprnc` did fail on that first build, and the `.pc` fix is now in the `Dockerfile`.** Measured 2026-10-15: `cprnc` died at link with ~200 undefined `H5*` symbols, its link line reading `f95 ... -o cprnc /usr/local/lib/libnetcdff.a /usr/local/lib/libnetcdf.a` -- the two archives and nothing else -- after the build log announced `NetCDF was not found, falling back on pkg-config`. `cime/CIME/non_py/cprnc/CMakeLists.txt:78` tries `find_package(NetCDF)`, which cannot succeed against an autotools-built stack with no `FindNetCDF.cmake` in cime, and falls through to `pkg_check_modules` at `:91-101`, which reads the `.pc` file's `Libs:` alone; netCDF's configure records the transitive `-lhdf5_hl -lhdf5 -lz -lcurl -lxml2` in `Libs.private`. Before 4b this directory held the *shared* parallel stack and `DT_NEEDED` covered it. The fix folds `Libs.private` into `Libs` in the installed `netcdf.pc` and `netcdf-fortran.pc`, so every pkg-config consumer gets a correct static link line and cime stays unmodified; the layer asserts `pkg-config --libs netcdf` contains `-lhdf5` afterwards. Leaving it to be discovered rather than pre-empting it cost one build and bought a measurement instead of a guess -- the same trade would be right again.
   - ✅ **Revalidated 2026-10-15** with all four wrapper scripts (`run-unit-tests-`, `run-test-`, `run-case-`, `run-sys-tests-in-container.sh`), `smoke-test.sh` and `smoke-test-pfunit.sh`. All pass.
   - ✅ **Acceptance test, on the EXECUTABLE.** A freshly built `cesm.exe` must contain neither `/ctsm/libraries/mpi-serial` nor `/ctsm/libraries/parallelio` among its compiled-source strings (`strings -a`). A passing version check without that is the old vacuous behaviour unchanged. The probe image met this on 2026-10-09, and **the real image met it on 2026-10-15**: both a `SMS_D_Ld1_Mmpi-serial.1x1_brazil` test build and a standalone `IHistClm60BgcQianRsGs` case build came back with zero hits for either submodule path and non-zero for `/usr/local/mpi-serial`.
   - ✅ **Republished to GHCR** as `:20261015`, 2026-10-15.
   - ✅ **Repinned `cirrus-testing.yml` to `:20261015` and added `_Mmpi-serial`**, together in one commit, as required: the tests carried no `_M` modifier and so took the machine default of `mpich`, which the image no longer has, while the old pin named an image whose paths `gnu_container.cmake` no longer matches. Only `simple-build-create_test` changed; `simple-build-run_sys_tests` runs in an unrelated mpich image. **Not yet exercised by a passing CI run** -- that job is currently red for a reason unrelated to 4b (item 13), and it builds a different test (`SMS_D_Ld3_Mmpi-serial.f10_f10_mg37.I1850Clm50BgcCrop`, global grid, `--no-run`) than anything validated on Casper, on `ubuntu-latest` rather than `gha-runner-ctsm`.
   - Consider an upstream `ccs_config` PR so `machines/container/config_machines.xml` declares `mpi-serial` in `<MPILIBS>`. mpi-serial works there today despite only `mpich` being listed, but the file is lying -- and now that mpich is gone from the image, `<MPILIBS>mpich</MPILIBS>` names the one thing the container cannot do. Shared with CESM, so it is separate work.
5. **The drift check.** One scheduled thing that compares *built executables* --
   container against derecho, gnu + mpi-serial -- given the latest CTSM master
   tag and derecho's current stack. Version metadata is not enough; it is what
   lets today's checker pass while the two sides link different mpi-serials.
   Reasoning, the instrument map and the measurements behind it:
   `DRIFT_CHECK_DESIGN.md`.

   **2026-10-26: CCPP-502 answered.** The CIRRUS admin added read-only mounts of `/glade/u/apps/derecho` and `/glade/u/apps/cesmdev` to `gha-runner-ctsm`, at their own paths -- not all of `/glade/u/apps`. Those two are the checker's two module roots, and derecho's `25.10` Spack tree and `spack.lock` live under the first, so they cover what the plan reads; no symlink under them escapes to an unmounted path except `cesmdev/modules/{23.06,23.09,24.12}`, which point into `/glade/u/apps/cseg` and matter only if ccs_config goes back to an env older than `25.10` (checked from a login node, 2026-10-26). **Not yet observed from the runner.** Appendix A of the plan (the cron fallback) is not needed unless the mounts turn out not to work.

   **2026-10-16: implementation plan written, execution deferred until CCPP-502 is answered.** The plan is `docs/superpowers/plans/2026-10-16-drift-check.md`. It records the decisions taken that day -- the executable comparison runs per CI build as a `strings` assertion in `cirrus-testing.yml`; the scheduled monitor compares the published image's as-built facts and the recipe against derecho's live module tree and `spack.lock`, and never builds a derecho binary; a hosted-runner heartbeat catches a silent self-hosted runner; the cron fallback is built only if the mount is refused -- and its Task 1 rewrites this item and `DRIFT_CHECK_DESIGN.md` to match. Until that task runs, the text below predates those decisions where the two differ.

   - **One checker script, two workflows.** The PR gate reads the PR branch and
     only in-repo files, runs on any hosted runner, and blocks the merge. The
     scheduled monitor reads the latest master tag and derecho's live stack,
     runs on `gha-runner-ctsm` (else a cron on an NCAR machine), and opens an
     issue with its own `GITHUB_TOKEN`. Sharing one script is what keeps their
     definition of "matching" from drifting apart.
   - **The cron fallback does not need a PAT, and should not have one.** Split
     it so the cron never talks to GitHub: the cron reads derecho's module tree
     and `spack.lock` and writes a small timestamped JSON report somewhere under
     `/glade/campaign/cesm/cesmdata/cseg`, which the Cirrus runner **already**
     mounts read-only; the scheduled workflow reads that file through the
     existing mount and opens the issue with `GITHUB_TOKEN`. No new mount, no
     stored credential -- the cron's only privilege is writing a file it owns.
     Two things this needs: the report carries a timestamp and the workflow
     fails when it is stale, or a dead cron is silent; and the path is chosen
     deliberately, not inside `inputdata`.

     This is the fallback's shape whether or not the `/glade/u/apps` mount
     lands. CIRRUS's own GitHub Actions runners already require a PAT, but it
     lives in NCAR's OpenBao secret store (`<ucar email>/github_pat`,
     fine-grained, per-repo, one-year expiry), not in a dotfile -- so "no
     credential on glade" matches the platform's own posture rather than
     fighting it. A GitHub App is not an improvement here: its private key is
     still a secret at rest.
   - **The two sides need different instruments**, because the container links
     the serial stack statically and derecho links it dynamically. No single
     tool covers both; the per-component mapping is in
     `DRIFT_CHECK_DESIGN.md`.
   - **Live reads need `/glade`, not Lmod.** Resolve `default` symlinks rather
     than taking the highest version, take the HDF5 pairing from `netcdf`'s
     `depends_on`, and run `nf-config` by absolute path. Loading derecho's stack
     off derecho is impossible -- `cray-mpich` needs Cray PE -- and unnecessary.
     `probe-derecho-modules.yml` encodes this.
   - **`gha-runner-ctsm` cannot reach `/glade/u/apps`, so it cannot host these
     reads as things stand.** Measured 2026-10-09 by running the probe:
     `/glade/u/apps/derecho/modules` is not visible from the runner. The
     hostname (`gha-runner-ctsm-*-runner-*`) says it is a Kubernetes pod, and
     glade reaches it as exactly **one** NFS mount, nothing else:

     ```
     gladedm1.ucar.edu:/gpfs/csfs1/cesm/cesmdata/cseg
         -> /glade/campaign/cesm/cesmdata/cseg   nfs4 ro
     ```

     `campaign` is the only entry under `/glade`. So the question was never
     "does the runner see glade" but "which parts", and the fix is another
     mount of the same shape rather than this item's cron fallback.

     **Requested 2026-10-09: <https://jira.ucar.edu/browse/CCPP-502>; answered 2026-10-26 with mounts of `/glade/u/apps/derecho` and `/glade/u/apps/cesmdev` (see the top of this item).** Until
     it was answered, neither the scheduled half of `derecho-version-check.yml`
     nor this item's live reads can run on the Cirrus runner. The cron
     fallback above is unblocked either way, so a refusal costs a different
     shape, not the feature.

     **Ask for read-only `/glade/u/apps`, one mount.** If it must be
     narrower, `/glade/u/apps/derecho` and `/glade/u/apps/cesmdev`. Do **not**
     ask for paths below those: everything the check wants lives under an
     `ncarenv` version directory (`/glade/u/apps/derecho/25.10/envs/public/`)
     or in a sibling `modules/` tree, so a narrower mount would have to be
     re-requested every time ccs_config bumps `ncarenv` -- which is one of the
     events this check exists to notice. Size is not the reason to narrow: it
     is an NFS mount, not a copy, and the check reads a handful of small
     files. The Spack **install** tree is not needed at all.

     Shape to hand the admin: the existing mount re-exports a subtree of the
     `csfs1` GPFS filesystem over NFSv4 from `gladedm1.ucar.edu`. This is the
     same thing for the `glade_user` filesystem, which is what `/glade/u`
     is.

     Until a mount or the cron lands, the scheduled half of
     `derecho-version-check.yml` would fail every week -- correctly, since an
     unreachable module tree is a loud failure by design, but uselessly.
   - **`[snapshot]`'s HDF5 and netCDF-Fortran need a live source here.**
     `check-derecho-versions.py` guards them only against a `netcdf` *version*
     bump. derecho can repoint an unchanged `netcdf/4.9.3` at a different HDF5
     and nothing in that check notices:
     `/glade/u/apps/derecho/modules/24.12/gcc/12.4.0/netcdf/4.9.3.lua` has
     `depends_on("hdf5/1.12.3")`, while
     `/glade/u/apps/derecho/modules/25.10/gcc/14.3.0/netcdf/4.9.3.lua` has
     `depends_on("hdf5/1.14.6")` -- and `hdf5` is never named in
     `config_machines.xml`, so the module-existence check does not look at it
     either. The scheduled check must read both live. The PR gate cannot: it
     runs `--skip-module-tree` and has no `/glade`.

     **Read them from `spack.lock`, not from `nf-config`.** The environment
     file `/glade/u/apps/derecho/<ncarenv>/envs/public/spack.lock` is the
     concretized DAG as JSON, and the 4-character token in a module's install
     prefix is a unique prefix of that spec's dag hash. So the serial gnu
     `netcdf/4.9.3`, whose prefix ends `/gcc/14.3.0/s6hn`, resolves to exactly
     one spec there, whose dependencies read (verified 2026-10-09):

     ```
     netcdf-fortran 4.6.2   netcdf-c 4.9.3   hdf5 1.14.6
     ```

     That is better than the `nf-config` plan in every way that matters here:
     one text file instead of the Spack install tree, no executing anything
     out of derecho's prefix, all three versions from one read, and the hash
     link pins them to the module the serial path actually loads rather than
     assuming the environment holds only one version of each.
   - **derecho's defaults moving is not drift and is not reported.**
     `config_machines.xml` pins exact versions, so a derecho run still uses the
     pinned one and the container still matches. Whether CTSM is keeping up with
     derecho's defaults has a different owner.
   - **PR gate hygiene:** drop `push: branches: ['**']` (see item 10) and add a
     staleness guard that fails only when the base branch moved in a way this
     check can see -- `ccs_config`, `Dockerfile`, `derecho-versions.ini`, the
     checker. Use the same path list the workflow declares in `paths:`, so guard
     and trigger stay defined in one place; `fetch-depth: 0`; compare against
     `github.event.pull_request.head.sha`, **not** `HEAD`, which always passes
     on a merge ref; read `base.ref` rather than hardcoding `master`.
   - **Scheduled workflows only fire from the default branch**, so this does
     nothing until merged to `master` and cannot be schedule-tested on a branch.
     Give it `workflow_dispatch` too.
   - A self-hosted runner that is offline makes a scheduled job queue rather
     than fail, so a missed run is silent. Worth a staleness check on the last
     successful run, or accept it.
   - Report only on drift; the Actions run history is the "it ran" record.
6. ✅ **The image is up to derecho's current stack.** The ARGs are bumped,
   `check-derecho-versions.py` passes all 10 checks, and the image was rebuilt
   and validated on Casper on 2026-10-07 -- `smoke-test.sh`,
   `smoke-test-pfunit.sh`, `run-unit-tests-in-container.sh` and
   `run-test-in-container.sh` all pass. Saved at
   `/glade/work/$USER/ctsm-ci-derecho-gnu_20261007.tar` (4.1 GB), published as
   `ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261007`, and pinned in
   `cirrus-testing.yml`. The published image's config digest is
   `sha256:be3badb712...`, matching the Casper build log, so what is on GHCR is
   demonstrably the artifact that was validated rather than a re-tag.

   | ARG | was | now (= derecho) |
   |---|---|---|
   | `GCC_VERSION` | 12.2.0 | 14.3.0 |
   | `NETCDF_C_VERSION` | 4.9.2 | 4.9.3 |
   | `NETCDF_FORTRAN_VERSION` | 4.6.1 | 4.6.2 |
   | `HDF5_VERSION` | 1.12.2 | 1.14.6 |
   | `PNETCDF_VERSION` | 1.12.3 | 1.14.1 |
   | `ESMF_VERSION` | 8.6.0 | 8.9.1 |
   | `MPI_SERIAL_VERSION` | 2.5.4 | 2.5.3 |
   | `PIO_VERSION` | 2.6.2 | 2.6.8 |

   `MPICH_VERSION` stays 3.4.3: cray-mpich 8.1.32 is still 8.1.x and still
   MPICH-3.4-ABI-derived, so only `[deviation_guard]` moved (8.1.27 ->
   8.1.32). `PFUNIT_VERSION` was already current.

   Two coupled edits went with the ARGs. `gnu_container.cmake`'s hardcoded
   `set(ESMFMKFILE ...)` moved to `esmf-8.9.1-mpiuni` -- the coupling the
   build-time assertion now guards. And **HDF5's upstream tag scheme changed**:
   `hdf5-1_12_2` became `hdf5_1.14.6`, so the URL, the tarball name and the
   extracted directory all had to change. Every download URL and extracted
   directory name for the new versions was checked against upstream; HDF5 was
   the only break. GCC's `ftp.gnu.org` URL could not be reached from a login
   node, but `mirrors.kernel.org` confirms 14.3.0 exists, and the same URL
   fails for 12.2.0 too, so that is egress and not a bad link.

   **Do not merge this ahead of the rebuild.** The check now passes against the
   Dockerfile, which is a statement about the recipe, not about what is on
   GHCR. Merging before the image is published and `cirrus-testing.yml` is
   repointed puts master in a state where CI is green while pulling an image
   that matches nothing in the repo.

   The rebuild needs the full revalidation chain (`smoke-test.sh`,
   `smoke-test-pfunit.sh`, `run-unit-tests-in-container.sh`, and at least one
   run wrapper), then a manual republish and a tag bump in
   `cirrus-testing.yml` -- nothing republishes automatically until item 8.

   **On the GCC 12.2.0 -> 14.3.0 jump, which was expected to be the risky
   part:** it was not. Nothing in the stack hit a GCC 14 diagnostic. The two
   things that did break were unrelated to the compiler, and both are recorded
   where they bite: HDF5 changed its upstream tag scheme between these releases
   (`hdf5-1_12_2` -> `hdf5_1.14.6`), and netCDF-C 4.9.3 dropped the transitive
   libraries from `nc-config --libs`, which only matters where the stack is
   linked statically. The original reasoning is kept below because it is still
   the right thing to check first on the next compiler jump.

   It is two major releases, and
   every other library in the image gets recompiled against it, so a new
   diagnostic anywhere in that chain stops the build. Expect the trouble on the
   **C** side, not the Fortran one. Per GCC 14's porting notes
   (https://gcc.gnu.org/gcc-14/porting_to.html), several long-standing warnings
   are errors by default in GCC 14: `-Wimplicit-function-declaration`,
   `-Wincompatible-pointer-types`, `-Wint-conversion` and `-Wreturn-mismatch`.
   That is what breaks old autotools `configure` scripts and old C sources. The
   `-fallow-argument-mismatch` workaround in `gnu_container.cmake` (see "Worth
   raising upstream" below) is *not* a concern: that file sets the flag
   unconditionally, precisely because ccs_config's version guard cannot fire
   there, and neither gfortran's argument-mismatch behavior (unchanged since
   10) nor the flag itself moved in GCC 14. Budget for a debug cycle, not a
   single clean build.

7. ✅ **`[snapshot]` can no longer go stale silently.** The two snapshot
   versions were re-measured on 2026-10-06 against `netcdf-mpi/4.9.3` under
   `ncarenv/25.10`: HDF5 **1.14.6** (now its own `hdf5-mpi` module, pulled in by
   `netcdf-mpi` through `depends_on`, where under `ncarenv/23.09` it was
   bundled) and netCDF-Fortran **4.6.2**. Both had drifted from the recorded
   1.12.2 / 4.6.1, so the two checks that had been printing ✅ were reporting
   stale as green. They now fail, and the ARG bump belongs in the item 6
   rebuild.

   `snapshot` mode now has the guard `deviation` already had.
   `derecho-versions.ini` records `measured_against_netcdf_mpi` alongside the
   values; the check reads derecho's live `netcdf-mpi` from
   `config_machines.xml` and, when the two differ, reports the snapshot as
   unusable and refuses to compare either ARG against it. A ✅ now means
   "measured against the module derecho has today" rather than "equal to a
   number someone typed once".

   Verified by replaying the original failure: with the pre-bump values
   (1.12.2 / 4.6.1), their matching ARGs, and provenance recorded as 4.9.2, the
   old code printed two passes and the new code fails all three lines. A
   missing provenance key raises rather than passing.

   What this does **not** cover, and item 5 still has to: drift in a module
   `config_machines.xml` does name, where derecho moves and nothing in the repo
   changes. This guard only fires once some commit makes the check run.

8. **Build the image in CI: build on PRs, build and publish on merges to
   master, for x86_64 and arm64.** Today the image is built by hand on Casper
   and pushed by hand (README "Publishing"), which is why item 6's rebuild is a
   manual errand and why `cirrus-testing.yml`'s pinned tag has to be bumped by
   hand too. Trigger on changes under `docker/ctsm-ci-derecho-gnu/**` (excluding
   the `.md` files) and on the workflow file itself; `workflow_dispatch` as
   well, since a rebuild is sometimes wanted with no file change (a base-image
   or upstream-tarball refresh).

   **Both architectures build natively, one job per architecture** --
   `ubuntu-latest` for x86_64 and `ubuntu-24.04-arm` for arm64. No QEMU: this
   image compiles GCC, MPICH, HDF5, netCDF-C/Fortran, PnetCDF, three ESMF trees,
   git and pFUnit from source, and emulating any of that is hours-to-days of
   wall clock. The shape that fits is the standard native multi-arch build:
   each job builds and, on master only, pushes by digest
   (`outputs: type=image,push-by-digest=true,name-canonical=true`), and a final
   job joins the digests with `docker buildx imagetools create`. On a PR the
   push and the join are skipped and the jobs only have to build.

   Note that the existing `docker-image-*.yml` workflows are not a template:
   they build the `ctsm-docs` image, which is being retired in favor of an
   external one, and `docker-image-build-publish.yml` asks a single
   `ubuntu-latest` job for `platforms: linux/amd64,linux/arm64`, which is the
   emulated form this must not use. The `platforms:` option is not itself the
   problem -- it is how a native per-architecture job names what it is
   building, and both jobs here will set it. What makes a build emulated is
   naming a platform the builder has no native node for, which is exactly what
   one x86_64 runner asked for `linux/arm64` does.

   **CIRRUS has only x86_64 runners today** (CIRRUS admin, 2026-10-16); Mac
   hardware is on order to provide organizational osx arm64 and linux arm64,
   with no date given. Two consequences. First, if the capacity measurement
   below sends the build to `gha-runner-ctsm`, that fallback covers x86_64
   only; arm64 stays on `ubuntu-24.04-arm` and so keeps whatever capacity
   problem pushed x86_64 off the hosted runners, until the CIRRUS arm64 runners
   exist. Revisit this when they land -- that would give the arm64 job a native
   self-hosted home. Second, the workaround the admin offered for the meantime
   -- add `platforms: linux/arm64` to the build-push action and build arm64 on
   an x86_64 runner -- is the QEMU route under another name. It is the right
   advice for an ordinary image and the wrong one for this one, which would be
   emulating the compilation of GCC, netCDF, ESMF and pFUnit. (The suggested
   snippet also omits `docker/setup-qemu-action`, without which that build does
   not run at all on an x86_64 runner -- which is itself the tell that it is
   emulation.)

   **The feasibility spike is written:**
   `.github/workflows/probe-runner-capacity.yml`, `workflow_dispatch` only,
   informational. Phase A (default, ~2 min) reports what a runner actually has
   on both hosted architectures and on `gha-runner-ctsm`, before and after
   reclaiming the preinstalled toolchains (hosted only -- the Cirrus runner is
   shared, persistent hardware); Phase B (`full_build: true`) attempts the real
   build and reports how far it gets and how long it takes. Run Phase A first:
   if the reclaimed disk on the hosted runners is not comfortably above what the
   build needs, the choice is between the Cirrus runner and the base-image split
   below.

   A job aimed at an offline or busy self-hosted runner queues rather than
   failing, and `timeout-minutes` does not bound queue time, so that row sitting
   pending while the hosted ones finish is itself the answer about availability,
   not a hung workflow.

   **The open feasibility question is capacity, and it should be measured
   before the workflow is designed in detail.** The build is about 50 minutes on
   16 native cores (README "Publishing"); GitHub-hosted runners are much
   smaller, and a job is killed at 6 hours. **Disk is the tighter limit**: a
   hosted runner has roughly 14 GB free against a 3.1 GB image whose build
   unpacks and compiles GCC, ESMF and the rest of the stack, so a
   reclaim step (dropping the runner's preinstalled Android/.NET/Haskell trees)
   and aggressive cleanup between layers are both likely required. Two more
   things follow: `MAKE_JOBS` (default 16) must be set from the runner's actual
   core count or the build thrashes, and layer caching is not optional. Prefer registry-backed cache
   (`cache-from`/`cache-to` with `type=registry`) over the GitHub Actions cache,
   which is capped at 10 GB per repo and small against this image. If a cold
   build will not fit in 6 hours even with cache, the fallback is to split the
   Dockerfile: a rarely-changing base image holding the compilers and libraries,
   rebuilt on demand, and a thin top layer rebuilt per change.

   **Do this after 4b, not before.** Stripping the image to serial removes
   MPICH, PnetCDF, the MPI-linked netCDF/HDF5 and two of the three ESMF trees,
   which cuts both the build time and the peak disk the measurement above is
   about. Measuring capacity against today's image would answer a question that
   no longer applies.

   **The arm64 image exists only for local development on Apple Silicon. The
   x86_64 image remains the sole canonical one**, because replicating derecho
   is this image's premise and derecho is x86_64. That settles several things:

   - `check-derecho-versions.py` and the whole version-matching argument in this
     directory describe the x86_64 image. The arm64 image is the same recipe on
     different hardware, not "the derecho stack", and its README section should
     say so rather than leaving a reader to assume a test passing there means
     anything about derecho.
   - `cirrus-testing.yml` keeps pinning x86_64 explicitly -- a per-arch tag or a
     digest, not the multi-arch manifest -- so no CI run can land on arm64 even
     if it is someday run on an arm64 runner.
   - A multi-arch manifest under one tag is still the right shape despite that,
     and is the reason to prefer it: an Apple Silicon developer runs the same
     `podman pull` as everyone else and gets the arm64 image without having to
     know the tag scheme, while CI's explicit pin keeps the canonical path
     unambiguous. A pull of the manifest tag fetches only the matching
     architecture -- the index is metadata, and the other architecture's layers
     are never transferred.
   - **The README must carry this caveat when the step lands:** on Apple
     Silicon, which image you get is decided by the architecture of the `podman
     machine` VM, not by macOS. The default ASi machine is arm64, which is the
     intent; but a machine created with x86_64 emulation silently resolves the
     same tag to the amd64 image and runs it emulated. That presents as "the
     container is mysteriously slow", not as an architecture mistake, so it
     needs saying outright, along with `podman pull --platform linux/amd64` as
     the way to ask for x86_64 deliberately -- which is also why the per-arch
     tags are a convenience and a pinning mechanism rather than a capability.
   - Validation scales to the role. The full chain (`smoke-test.sh`,
     `smoke-test-pfunit.sh`, `run-unit-tests-in-container.sh`, a run wrapper)
     stays an x86_64 gate. For arm64 the bar is lower but not zero: run the
     smoke tests on the arm64 runner before its digest joins the manifest, so a
     plainly broken developer image cannot be published. The arm64 runner is the
     only place that can be done at all.

   The Dockerfile's own build-time assertions (the pFUnit prefix, the two
   esmf.mk checks) are architecture-independent and run unchanged on both --
   the cheapest evidence that an arm64 build is coherent.

   **`README.md` has to be reconciled with this.** Its "Publishing" section
   currently says, of multi-architecture manifests, **"Do not try that here"**,
   for two reasons. The first -- QEMU is 10-20x slower per core, and is
   unavailable on Casper anyway, which has no binfmt handlers and no root to
   register them -- remains exactly right *for building on Casper* and should be
   kept as such, not deleted; it simply does not apply to a native arm64 runner.
   The second, that arm64 is not derecho, is answered by labeling rather than by
   not building, as above.

   GHCR needs no new secret: the package is already public, so `packages: write`
   on `GITHUB_TOKEN` is enough.
9. **Make the `create_test` CI job actually run the test, not just build it.**
   `simple-build-create_test` in `cirrus-testing.yml` passes `--no-run`, so CI
   compiles `cesm.exe` and stops; nothing in CI has ever run the model. Not yet
   thought through -- it overlaps item 2 (inputdata and the `/glade` question
   on `gha-runner-ctsm`).
10. **Require manual runs on PRs instead of running on every push.**
    `cirrus-testing.yml` triggers on `[push, pull_request]`, so every push to
    every branch takes the self-hosted runner. Wanted: these run when asked
    for, not automatically. Not yet thought through -- which trigger, and
    which workflows it should cover.

    `derecho-version-check.yml`'s `push: branches: ['**']` has a second,
    stronger reason to go, and item 5 removes it: a `push` event checks out the
    branch itself with no merge commit, so it compares a stale `ccs_config`
    against a stale `Dockerfile` and passes. That is the one place the check
    gives a wrong answer rather than merely running too often.
11. **If not already there: Add a way for users to customize bind mount targets**
12. **Add documentation.** User-facing docs beyond this directory's README;
    separate from the README updates folded into the items above.
13. **Both `cirrus-testing.yml` jobs are currently red, for two unrelated reasons.** Neither is caused by item 4b; both predate the `:20261015` repin, and the first is confirmed to predate it.

    **`simple-build-run_sys_tests` dies before it reaches CTSM at all**, at the `git checkout -b tmp` that opens its "Call run_sys_tests" step:

    ```
    fatal: detected dubious ownership in repository at '/__w/CTSM/CTSM'
    Error: failed to run script step: command terminated with non-zero exit code: exit status 128
    Error: Executing the custom container implementation failed. Please contact your self hosted runner administrator.
    ```

    That is git's `safe.directory` protection: the checked-out workspace is owned by a different uid than the user git runs as inside the container. No workflow in this repo sets `safe.directory` anywhere. The fix is a `git config --global --add safe.directory "$GITHUB_WORKSPACE"` -- but **put it inside the same step**, not in a step of its own: this job runs under the k8s container hook (`/home/runner/k8s/index.js`), where `HOME` can differ between steps, so a `--global` config written in an earlier step need not be visible later. That is also the likely reason `actions/checkout`, which normally registers the workspace itself, did not prevent this. Worth checking whether the submodules need `'*'` rather than just the workspace, though `bin/git-fleximod update -o` in the preceding step appears to survive as-is.

    **`simple-build-create_test` fails building cprnc: `buildlib.cprnc FAILED`.** Measured 2026-10-16 from the uploaded test directory: the link line is `gfortran ... -o cprnc /usr/local/lib/libnetcdff.a /usr/local/lib/libnetcdf.a`, then undefined `H5*` symbols. Compiler (`/opt/gcc/14.3.0`), the genf90 clone and the `PKG_CONFIG_*` environment were all fine. The cause is the cprnc version: cime6.5.5 pins cprnc **v1.1.5**, whose `CMakeLists.txt` added a `find_library(netcdf)` fallback that runs *before* pkg-config and links the bare static archives, so the `.pc` fix never comes into play. **The 4b Casper validation never saw this** because this checkout's `cime/CIME/non_py/cprnc` was stale at v1.1.2 (`git submodule status` showed `+9cbdce4`); run `git -C cime submodule update CIME/non_py/cprnc` before any revalidation.

    Fix chosen (2026-10-26): build cprnc once in the image and point the case at it via `CCSM_CPRNC`, as derecho does. CIME can only override a machine value by replacing the whole machine file, so `cime-machines/container/config_machines.xml` is a stopgap copy of ccs_config's with that one line added, baked into `/root/.cime/container/`; `cirrus-testing.yml` fails if it drifts from ccs_config's. The proper fix is the same line upstream in ccs_config. Status: proven 2026-10-26 with a throwaway probe image on top of `:20261015` -- the case recorded `CCSM_CPRNC=/usr/local/bin/cprnc`, built no cprnc, and `SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs` passed every phase through `SHORT_TERM_ARCHIVER`, with the `strings` acceptance test clean (0/0, 6 `/usr/local`). Then folded into `Dockerfile`. **Rebuilt and validated on Casper 2026-10-26** (`build-and-validate.sh`, log `build-and-validate.log.20261026132645`): smoke test, pFUnit smoke test, unit tests, `create_test`, `run_sys_tests` and a case all pass; the `create_test` case recorded `CCSM_CPRNC=/usr/local/bin/cprnc` and built no cprnc; `strings` acceptance test 0/0, 6 `/usr/local`. Saved as `/glade/work/$USER/ctsm-ci-derecho-gnu_20261026.tar`. **Still to do:** republish and repin `cirrus-testing.yml`.

    Also measured on that run: `HOME` was `/root` inside the CIRRUS job container, not `/github/home`, so the image's baked `/root/.cime` is what CI uses; and the job container sees `/glade/campaign/cesm/cesmdata/cseg` (now used as `CESMDATAROOT`).

    Related but distinct: item 9 is about that job building without running, and item 10 about when these jobs fire at all. Both assume the jobs pass first.

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
