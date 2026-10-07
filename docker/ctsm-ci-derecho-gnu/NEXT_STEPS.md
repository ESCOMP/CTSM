# Next steps: ctsm-ci-derecho-gnu container

_Last updated: 2026-10-06. The image builds and validates end-to-end on
Casper -- pFUnit and CTSM's Fortran unit tests (55/55), plus single-point
**runs** with both mpi-serial and mpich -- and is **published and public** at
`ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20260831`, which `cirrus-testing.yml`
pins. (The earlier `:20260830` tag predates the serial netCDF stack and does
NOT work with the current `gnu_container.cmake`.) All three wrapper scripts
have now been exercised on Casper. The version-check question is now settled
(item 4), and settling it showed the image has fallen behind derecho: as of
`ccs_config_cesm1.0.88` the checker had been failing 9 of its 10 checks. The
ARGs are now bumped, all 10 pass, and the image has been **rebuilt and
validated on Casper (2026-10-07)** — but it is **not yet published**, and
`cirrus-testing.yml` still pins the August image, so this must not merge until
both are done. What is left: publishing and repointing CI (item 6), building
the image in CI instead of by hand (item 8), and the Phase 2 drift cron
(item 5)._

## Where things stand

- **The from-scratch image builds and is VALIDATED on Casper.** `Dockerfile`
  (FROM `almalinux:9`) builds the full gnu stack — GCC 12.2.0, MPICH 3.4.3
  ch4:ofi, HDF5 1.12.2, netCDF-C 4.9.2, netCDF-Fortran 4.6.1, PnetCDF 1.12.3,
  ESMF 8.6.0 (debug + optimized), plus **git built from source**. Tagged
  `localhost/ctsm-ci-derecho-gnu:dev`. Versions matched derecho's
  `ncarenv/23.09` gnu stack as of the 2026-08-31 build; derecho is now on
  `ncarenv/25.10` (`config_machines.xml:41`) and the two no longer agree — see
  the README table and item 6.
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
4. ✅ **`MPI_SERIAL_VERSION` and `PIO_VERSION` are now checked against
   derecho.** They were the only `Dockerfile` ARGs that nothing checked, so
   they could drift silently. `check-derecho-versions.py` now checks both
   `direct` against derecho's serial stack -- the `mpilib="mpi-serial"`
   `<modules>` blocks of `config_machines.xml` -- and that immediately showed
   **both are stale**: `MPI_SERIAL_VERSION=2.5.4` vs derecho's
   `mpi-serial/2.5.3`, and `PIO_VERSION=2.6.2` vs derecho's
   `parallelio-serial/2.6.8`. Clearing them needs an image rebuild and
   revalidation, which is a separate decision; the ARGs are untouched for now.

   Deciding that the yardstick is derecho, not CTSM's `.gitmodules`, is the
   governing policy for the whole image: it should be up to date with what you
   get from CTSM's current `ccs_config` plus derecho's default software stack.
   The tracing below is why that needed a decision at all. Its mpi-serial half
   was recorded wrong the first time and is corrected here (2026-10-01); the
   PIO half was rewritten too, gaining the citations it had lacked, but its
   conclusion is the same as before.

   The question it answers is *what versions get used when doing a serial CTSM
   test on derecho*. Traced 2026-10-01:

   - **mpi-serial: 2.5.3, from derecho's own `mpi-serial/2.5.3` module** --
     *not* from CTSM's `libraries/mpi-serial` submodule, as an earlier version
     of this item claimed. `MPI_SERIAL_PATH` is a **cmake macro, not an XML
     setting**, so looking for it in `config_machines.xml` finds nothing and
     proves nothing. `ccs_config/machines/derecho/derecho.cmake:4` sets it
     unconditionally:

     ```cmake
     set(MPI_SERIAL_PATH "$ENV{NCAR_ROOT_MPI_SERIAL}")
     ```

     `module load mpi-serial/2.5.3` is what populates
     `NCAR_ROOT_MPI_SERIAL`. `cime/CIME/BuildTools/configure.py` turns the
     cmake macros into `Macros.make`, and `cime/CIME/Tools/Makefile:911-918`
     then does, under `ifeq ($(MPILIB),mpi-serial)` / `ifdef
     MPI_SERIAL_PATH`:

     ```make
     MPISERIAL = $(MPI_SERIAL_PATH)/lib/libmpi-serial.a
     MLIBS += -L$(MPI_SERIAL_PATH)/lib -lmpi-serial
     ```

     The `else` branch -- CTSM's submodule, built into
     `$(INSTALL_SHAREDPATH)/lib/libmpi-serial.a` -- is reached only when the
     variable is *empty*, i.e. where no mpi-serial module is loaded. That is
     the **container's** case, not derecho's. So derecho links 2.5.3 into
     `cesm.exe`.
   - **PIO: 2.6.8, from derecho's `parallelio-serial/2.6.8` module.** Same
     shape: `derecho.cmake:8-11` sets `PIO_LIBDIR`/`PIO_INCDIR` from
     `$ENV{PIO}` when that module has been loaded, and
     `cime/CIME/Tools/Makefile:463-472` honors an external `PIO_LIBDIR`,
     falling back to `$(INSTALL_SHAREDPATH)/lib` only when it is unset. CTSM's
     `libraries/parallelio` submodule is pinned at `pio2_6_8`, so the two agree
     today.

   **That makes the conclusion stronger than "policy says so."** 2.5.3 is what
   a serial CTSM build on derecho actually links, so `MPI_SERIAL_VERSION=2.5.4`
   is genuine drift from derecho on the merits, not merely by the yardstick
   convention; the same holds for `PIO_VERSION=2.6.2` against 2.6.8.

   Two facts recorded here previously were wrong, both because the ctsm5.4.054
   rebase moved `ccs_config` from `ccs_config_cesm1.0.48` to
   `ccs_config_cesm1.0.88`: derecho's mpi-serial is **2.5.3**, not 2.3.0, and
   its `parallelio-serial/2.6.8` **is** in `config_machines.xml`, not absent
   from it.

   **Neither ARG reaches a CTSM case build _in the container_.** They exist
   only to give the mpiuni ESMF an external PIO: the container's mpi-serial
   installs under its own prefix specifically so it cannot shadow the case
   build -- nothing sets `MPI_SERIAL_PATH` there, so the Makefile's `else`
   branch above compiles `libraries/mpi-serial` -- and its PIO installs under
   `${SERIAL_PREFIX}` without setting `PIO_LIBDIR`, so a case build in the
   container compiles `libraries/parallelio` for itself. **On derecho both of
   the corresponding modules _are_ linked**, per the tracing above, which is
   why derecho rather than the container's own internals decides what counts
   as drift here.

   That made CTSM's `.gitmodules` fxtags a plausible yardstick, and an earlier
   version of this item picked them. It is superseded twice over: what the
   image promises is to *mirror derecho*, and a serial library the image ships
   under a version derecho no longer has is drift whether or not the
   *container's* case build links it -- and, per the correction above,
   derecho's case build does link its own. So both ARGs are compared to
   derecho, and the fxtags are recorded here only as context:

   | ARG | value | derecho (serial) | CTSM fxtag |
   |---|---|---|---|
   | `MPI_SERIAL_VERSION` | 2.5.4 | 2.5.3 ❌ | `MPIserial_2.5.4` |
   | `PIO_VERSION` | 2.6.2 | 2.6.8 ❌ | `pio2_6_8` |

   Both therefore want an image rebuild, now tracked as item 6 below along
   with the four other stale ARGs and the `cray-mpich` deviation guard that the
   same check reports; no fourth check mode was needed.
5. **Phase 2 drift detection** (see `derecho-versions.ini`): a cron on
   Casper/Derecho reading live derecho versions, opening a GitHub issue on
   drift and emailing on success. Planned as one of the last steps.
6. **Bring the image up to derecho's current stack: rebuild.** The ARGs are
   bumped, `check-derecho-versions.py` passes all 10 checks, and the image
   **has been rebuilt and validated on Casper (2026-10-07)**: `smoke-test.sh`,
   `smoke-test-pfunit.sh`, `run-unit-tests-in-container.sh` and
   `run-test-in-container.sh` all pass. The saved image is
   `/glade/work/$USER/ctsm-ci-derecho-gnu_20261007.tar` (4.1 GB).

   **What remains: publish it and repoint CI.** The published
   `:20260831` image is still the 2026-08-31 `ncarenv/23.09` build, and
   `cirrus-testing.yml:55` still pins it, so this item stays open until the new
   image is pushed to GHCR and that pin is bumped. Until then the repo
   describes a stack that nothing in CI actually runs.

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
   external one, and their single-step `platforms: linux/amd64,linux/arm64`
   is the emulated form this must not use.

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
   hosted runner has roughly 14 GB free against a 3.5 GB image whose build
   unpacks and compiles GCC, three ESMF trees and the rest of the stack, so a
   reclaim step (dropping the runner's preinstalled Android/.NET/Haskell trees)
   and aggressive cleanup between layers are both likely required. Two more
   things follow: `MAKE_JOBS` (default 16) must be set from the runner's actual
   core count or the build thrashes, and layer caching is not optional. Prefer registry-backed cache
   (`cache-from`/`cache-to` with `type=registry`) over the GitHub Actions cache,
   which is capped at 10 GB per repo and small against this image. If a cold
   build will not fit in 6 hours even with cache, the fallback is to split the
   Dockerfile: a rarely-changing base image holding the compilers and libraries,
   rebuilt on demand, and a thin top layer rebuilt per change.

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
