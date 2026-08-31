# Next steps: ctsm-ci-derecho-gnu container

_Last updated: 2026-08-31. The image builds and validates end-to-end on
Casper -- pFUnit and CTSM's Fortran unit tests (55/55), plus single-point
**runs** with both mpi-serial and mpich -- and is **published and public** at
`ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260831`, which `cirrus-testing.yml`
pins. (The earlier `:20260830` tag predates the serial netCDF stack and does
NOT work with the current `gnu_container.cmake`.) All three wrapper scripts
have now been exercised on Casper. What is left: a unit-test CI job, wiring
runs into CI, the version-check question below, and the Phase 2 drift cron._

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

Gotchas worth not rediscovering:

- **ESMF errors go to `PET0.ESMF_LogFile`, not `cesm.log`.** The PIO failure
  presented as a totally silent crash for hours because nothing looked there.
  Check it first for any run that dies during initialization.
- mpi-serial's `make install` is broken upstream (`$(INSTALL)` and
  `$(MKINSTALLDIRS)` are never defined), so the artifacts are copied by hand,
  as CIME's own `buildlib.mpi-serial` does.
- PIO compiled against mpi-serial emits real `MPI_*` symbol names, which
  ESMF's renamed mpiuni stubs do not satisfy -- hence `-lmpi-serial` in
  `ESMF_PIO_LIBS`. Only the archive is copied into `/usr/local/serial/lib`;
  the header stays out so it cannot shadow the case build's own `mpi.h`.
- The image bakes a snapshot of `gnu_container.cmake`, which silently shadowed
  edits during debugging. All three wrapper scripts now prefer the mounted
  repo's copy.

`Dockerfile.serial-netcdf` was the throwaway probe used to prove all this. It
has been folded into `Dockerfile` and deleted, exactly as `Dockerfile.pfunit`
was before it. Rebuilt and re-validated from scratch on 2026-08-31; saved to
`/glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar`, **published** as
`ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260831`, and pinned in
`cirrus-testing.yml`.

One publishing gotcha, since it presents as something else entirely: if
`podman load` is OOM-killed, it exits **137** having printed nothing, and the
next `podman tag` fails with the misleading `image not known`. `vfs` storage
duplicates every layer, so a 3.9 GB archive peaks around 54 GB -- load it from
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

1. **Add a unit-test job to `cirrus-testing.yml`.** Now unblocked. It needs
   the `$HOME/.cime` copy step (GHA overrides `HOME`); see README "Running
   CTSM's unit tests".
2. **Wire runs into CI.** `simple-build-create_test` currently runs on
   `ubuntu-latest`, which has no `/glade` at all, so a run job has to move to
   `runs-on: gha-runner-ctsm` and bind-mount inputdata into the container job.
   Whether container jobs on that runner see `/glade` is **unknown** -- the
   existing `list-glade-cesm-input` job proves only that the *host* does. Worth
   a probe workflow in the style of `probe-derecho-modules.yml`. Now that
   `run_sys_tests --wait` blocks on every launched test and exits nonzero if
   any of them failed, CI has the exit code it needs to fail the job on a
   test failure -- that piece is no longer a gap.
3. **Validate `run-sys-tests-in-container.sh` on Casper.** Written, syntax
   checked, and dry-run-verified on the login node only (see "Added
   2026-08-31: run_sys_tests wrapper" above) -- it has not yet been run
   against the actual container on a compute node. Each step below is
   ordered by what it actually proves; do them in order and do not skip one
   because a later step looks like it would cover it too.

   1. `cime/CIME/scripts/query_testlists.py --xml-machine derecho
      --xml-category aux_clm_mpi_serial --xml-compiler gnu --count` -- real
      suite resolution against derecho's testlist, on the host, in seconds,
      no container or allocation needed. This is the step that actually
      calls `get_tests_from_xml`; nothing below does.
   2. The wrapper with `-s aux_clm_mpi_serial --dry-run`. This does **not**
      prove suite resolution -- the wrapper always injects
      `--suite-compiler gnu`, which makes `run_sys_tests` skip
      `_get_compilers_for_suite`, the only caller of `get_tests_from_xml`,
      and `--dry-run` stops `create_test` from running at all. What it does
      prove: `run_sys_tests` imports and runs under the container's
      `python3`; `create_machine("ctsm-ci-container")` resolves;
      `git`/`bin/git-fleximod status` succeed against the bind-mounted
      `/ctsm`; and the testroot is named as predicted
      (`tests_<MMDD-HHMMSS>ct`).
   3. One known-good test through the wrapper --
      `-t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc`, then `echo $?`.
      This is the **only** step that proves `--wait` actually blocks and
      that the exit status propagates through podman to PBS -- the entire
      reason the `--wait` work exists, and until this runs it is untested
      outside unit tests against a fake launcher. Confirm the testroot
      appears under `$SCRATCH/cases_devcontainer/` and the case PASSes.
   4. A deliberately failing test, to confirm the exit status is nonzero --
      `-t FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` is a convenient
      choice, since it is already known to fail here for a documented reason
      (missing python modules; see the `-s` / ctsm_pylib note in README
      "Running run_sys_tests"). Confirm `echo $?` is nonzero.
   5. The full `-s aux_clm_mpi_serial`. Judge this run by `cs.status.fails`
      inside the testroot, **not** by the wrapper's exit code: a nonzero
      exit is *expected* on a first full run of this suite, for reasons
      unrelated to this change --
      `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` needs python modules
      the image lacks (same as step 4), and the suite's NEON/`CLM_USRDAT`
      and FATES entries need user datasets and FATES build support that this
      change does not touch. Do **not** use `-s clm_short` as a substitute
      "quick suite" -- it has exactly two derecho/gnu entries,
      `ERP_D_P64x2_Ld3.f10_f10_mg37.I1850Clm50BgcCrop` and
      `ERS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop`, neither mpi-serial, and
      `P64x2` wants 64 MPI tasks against a machine config with
      `MAX_MPITASKS_PER_NODE=4` inside an 8-cpu PBS reservation; it will
      fail for reasons that have nothing to do with this change.

   **What to watch for**, highest-risk failure modes first:
   - `_record_git_status` can abort `run_sys_tests` up front if git's
     dubious-ownership / `safe.directory` check trips on the bind-mounted
     repo. Mitigation: pass `--skip-git-status`.
   - The silent job log described in README "Running run_sys_tests" (Test
     output does not appear in the job log): with `--wait`, nothing streams
     to the PBS log between "Running: <create_test ...>" and the final exit
     code, so a long quiet job is expected, not a hang -- watch
     `<testroot>/STDOUT.<testid>` / `STDERR.<testid>` instead.
   - `MAX_MPITASKS_PER_NODE=4` together with `GMAKE_J=4` means CIME may build
     and run up to 4 tests at once inside the 8-cpu PBS reservation these
     wrappers request -- expect concurrency, not one test running at a time.
   - The `podman load` OOM (exit 137, prints nothing, `podman tag`/`podman
     run` then fail with the misleading "image not known") already
     documented above under "Resolved 2026-08-31: mpi-serial runs work"
     applies here too: load the image from a session with real memory before
     running any of the above.
4. **Decide whether `MPI_SERIAL_VERSION` and `PIO_VERSION` belong in
   `check-derecho-versions.py`.** Both are new `Dockerfile` ARGs that nothing
   checks, so they can drift from derecho silently -- they are the only ARGs in
   that position. Neither is a plain `direct` check:
   - derecho's mpi-serial is **2.3.0**; the container deliberately uses
     **2.5.4**, the version CTSM's own `libraries/mpi-serial` submodule pins.
     That is the `deviation` shape (assert derecho still says 2.3.0, recorded
     under `[deviation_guard]`, without comparing the ARG) -- the same shape as
     `cray-mpich` -> MPICH.
   - `PIO_VERSION` 2.6.2 *does* match derecho, but its parallelio module lives
     in the mpi-serial hierarchy and is absent from `config_machines.xml`, so
     it would have to be a `snapshot` entry.

   The prior question is a design call, not a mechanical addition, which is why
   this was left undone: for the serial stack, do we **track derecho** or
   **track CTSM's submodules**? They disagree today, and the answer decides
   which mode each ARG gets.
5. **Phase 2 drift detection** (see `derecho-versions.ini`): a cron on
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
