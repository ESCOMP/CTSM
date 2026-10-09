# Serial-only container, and a drift check that compares binaries

## Context

The goal: **one** thing, on a schedule (GitHub Workflow on CIRRUS, else a cron
on derecho), that checks — given the latest CTSM master tag and derecho's
current stack — that a CTSM **standalone run, gnu + mpi-serial only**, uses the
same software in the container as on derecho.

Two things changed during planning, and both shrink the work:

1. **The image drops mpich entirely.** Scope becomes serial-only by
   construction rather than by policy. The accepted-exception list loses MPI and
   keeps cray-libsci plus the distro runtime.
2. **The check compares built executables**, not version strings in metadata —
   because metadata comparison is demonstrably blind to drift that exists today.

### The evidence that forced (2)

NEXT_STEPS item 4 already states that `MPI_SERIAL_VERSION` and `PIO_VERSION`
never reach a case build in the container, then compares those ARGs to derecho
anyway. That hides live drift:

| | what a gnu+mpi-serial run actually links | from |
|---|---|---|
| container | `MPIserial_2.5.4` | CTSM submodule, compiled during the case build |
| derecho | `mpi-serial/2.5.3` | derecho's module, via `MPI_SERIAL_PATH` |

**`check-derecho-versions.py` is green on this today**, because it compares
`MPI_SERIAL_VERSION=2.5.3` against derecho's `2.5.3` — a perfect match
describing a copy of mpi-serial no case build ever links.

PIO has identical structure and currently agrees (`pio2_6_8` vs
`parallelio-serial/2.6.8`).

Do **not** try to read mpi-serial's version from the binary.
`libraries/mpi-serial/mpi.c:16` hardcodes `"mpi-serial 2.5.0"` regardless of
tag. PIO embeds no version string at all.

---

## Closing the mpi-serial gap: change the mechanism

**Decided.** Set `MPI_SERIAL_PATH` and `PIO_LIBDIR`/`PIO_INCDIR` in
`cime-macros/gnu_container.cmake`, mirroring
`ccs_config/machines/derecho/derecho.cmake:4,8-11`, so a container case build
**links the image's mpi-serial and PIO** instead of compiling CTSM's submodules
— the same mechanism derecho uses. CIME's `Tools/Makefile:911-918` and
`463-472` honor both when set and fall back to the submodules only when they are
empty, so this is the switch the container has been leaving unset.

Consequences:

- the difference stops existing rather than being reported, and the failure mode
  moves from "red check nobody can clear" (CTSM pins the submodule repo-wide,
  the module is NCAR's) to one file in this repo
- `MPI_SERIAL_VERSION` and `PIO_VERSION` become load-bearing instead of
  vestigial, so the checker's existing ARG-vs-module comparison becomes
  meaningful rather than vacuous — no fxtag mode is needed
- it reverses a deliberate earlier decision. NEXT_STEPS records that the image's
  mpi-serial installs under its own prefix *specifically* so it cannot shadow
  the case build. That reasoning was written when CTSM's submodules were the
  yardstick; under "mirror derecho" it points the other way.

**A suspected safety argument for this was checked and does not hold**, so it is
not among the reasons. The worry was two mpi-serial implementations in one
process — the image's 2.5.3 inside the mpiuni ESMF, CTSM's 2.5.4 in the
executable — which README's mpi-serial section warns against. Measured
2026-10-09:

- `cesm.exe` defines 176 `MPI_*` symbols and imports none; the only mpi-serial
  source compiled in is `/ctsm/libraries/mpi-serial`, and `/usr/local/mpi-serial`
  is referenced zero times
- `libesmf.so` (28 MB, 20835 dynamic symbols) **neither exports nor imports any
  `MPI_*` symbol**, so its mpiuni stubs have internal linkage and never
  participate in symbol resolution with the executable

Two copies of MPI code exist in the process but cannot interfere, which is what
matters. (Scope of that claim: `nm -D` sees only the dynamic symbol table, so it
establishes non-interaction, not that only one copy of the code exists.)

The image's PIO was built only as the mpiuni ESMF's external PIO, so confirming
it is a complete enough install for a case build is part of Phase 4's
revalidation, not an assumption.

---

## Settled

**Scope: serial only.** Drop mpich from the image. A gnu+mpi-serial run on
derecho loads `gcc/14.3.0`, `cray-libsci/25.03.0`, `ncarcompilers/1.2.0`,
`netcdf/4.9.3`, `mpi-serial/2.5.3`, `parallelio-serial/2.6.8` (`-debug` when
`DEBUG=TRUE`), `esmf/8.9.1`, plus the common `ncarenv/25.10`, `cesmdev/1.0`,
`conda/latest`, `nco`, `craype`, `cmake`.

**Accepted exceptions, now two:**

| | derecho | container | treatment |
|---|---|---|---|
| BLAS/LAPACK | `cray-libsci/25.03.0` | reference `lapack`/`blas` | **guarded** — record the version, fail when derecho moves |
| distro runtime | `GCC: (SUSE Linux) 13.3.1` | `GCC: (GNU) 11.5.0 (Red Hat)` | **documented only** — C runtime startup objects; unfixable short of rebasing on SLES |

cray-libsci is nearly free to guard: it's pinned in `config_machines.xml`'s
`compiler="gnu"` block, so it's one `[deviation_guard]` entry plus one `CHECKS`
row using the existing `deviation` mode.

**`ccs_config` naming a module derecho no longer has must FAIL.** The checker
never contacts derecho — "read live" means from `config_machines.xml`, which is
CTSM's *claim*, never verified. If ccs_config pins a version derecho has
removed, a derecho run dies at `module load`, there is no run to compare
against, and the check currently reports a match.

**derecho's defaults moving is NOT drift, and is not reported.**
`config_machines.xml` pins exact versions, so if derecho's default moves while
the pinned version still exists, a derecho run still uses the pinned one and the
container still matches. "Is CTSM keeping up with derecho?" has a different
owner; mixing it in would make a red check ambiguous about who must act.

**One script, two workflows.** Different subjects, intrinsic refs:

| | ref | derecho's side from | runner | on failure |
|---|---|---|---|---|
| PR gate | the PR branch (merge ref) | `ccs_config`, in-repo | any hosted | block the merge |
| scheduled monitor | latest master tag | live, from `/glade` | `gha-runner-ctsm`, else cron on derecho | open an issue |

CTSM tags every merge to master (`git describe origin/master` → bare
`ctsm5.4.050`), so "latest master tag" and "master HEAD" are the same commit in
practice. Sharing one script keeps their definition of "matching" from drifting
apart, and means the PR trigger exercises most of the scheduled path's logic.

Retire the "Phase 1 / Phase 2" vocabulary — "Phase 1" is defined nowhere.

---

## The instrument map

Measured 2026-10-08 against a `_Mmpi-serial` `container_gnu` binary under
`cases_devcontainer/` and a `_Mmpi-serial` `derecho_gnu` binary under
`tests_1007-104555de/`.

**The two sides are not symmetric** — the container statically links the serial
stack on purpose, derecho links it dynamically. No single tool covers both.

| component | derecho | container |
|---|---|---|
| gcc | `readelf -p .comment` | `readelf -p .comment` |
| ESMF | RPATH → `esmf/8.9.1` | RPATH → `/usr/local/esmf-8.9.1-mpiuni/lib` (version **and** flavor) |
| netCDF-C | RPATH; `ldd` path | `strings` → `4.9.3 of Oct  7 2026 …` |
| HDF5 | RPATH; `ldd` path | `strings` → `HDF5 library version: 1.14.6` |
| netCDF-Fortran | `nf-config --version` by absolute path | ❌ gap — needs a build-time LABEL |
| mpi-serial | RPATH → `mpi-serial/2.5.3` | version from the ARG (load-bearing after the switch); the binary's job is to **confirm the mechanism** — see below |
| PIO | RPATH → `parallelio/2.6.8` | same |
| BLAS/LAPACK | RPATH → `libsci/25.03.0` | `ldd` → `liblapack.so.3` (ABI only) |
| compression filters | `ldd` closure + `libnetcdff.so.7`'s RPATH (szip 2.1.1, c-blosc 1.21.6, zstd 1.5.7, bzip2 1.0.8, lz4 1.10.0) | ❌ static → `ldd` blind; needs a build-time LABEL |

Tool properties that constrain the design:

- `readelf -d` lists **direct** `NEEDED` only (13 libraries on derecho's
  binary); `ldd` walks the **transitive closure** (35). Every compression filter
  — and `libhdf5` itself — appears only in the closure. `readelf` cannot make
  the filter comparison.
- RPATH is **per file**. `cesm.exe`'s RPATH gives the netCDF *bundle* (4.9.3);
  `netcdf-fortran/4.6.2` appears only in `libnetcdff.so.7`'s own RPATH. That
  works but depends on NCAR's Spack bundle packaging — prefer `nf-config
  --version`, a supported interface.
- RPATH also records directories for **statically** linked libraries, which is
  why derecho's mpi-serial version is visible there and nowhere in `ldd`.
- `strings`, `readelf` and `readlink` only read the file and work from any
  machine with glade. **Only `ldd` resolves**, so it needs the right environment
  — inside the image for the container binary, modules loaded for derecho's.
- Version banners live wherever the library's *code* lives: inside `cesm.exe`
  when static, inside the `.so` when dynamic.

**Same version ≠ same build.** derecho's netCDF 4.9.3 was built 2025-10-27, the
container's 2026-10-07.

**Confirming the mpi-serial/PIO switch took effect.** Don't expect RPATH to
record them in the container: both are static archives, so unlike derecho's
Spack build nothing adds an `-rpath` entry for them. The reliable signal is the
compiled-source path in debug info. Today `strings -a cesm.exe` reports
`/ctsm/libraries/mpi-serial` — CTSM's submodule, compiled during the case build.
After the switch that string must be **absent**, replaced by the image's own
build path. That transition is the acceptance test for the mechanism; the
version numbers then follow from the ARGs, which the Dockerfile already
controls.

---

## Work

### Phase 1 — strip the image to serial

`docker/ctsm-ci-derecho-gnu/Dockerfile`:

- remove the MPICH build and `MPICH_VERSION`; drop `/opt/mpich/...` from `PATH`
- remove PnetCDF entirely and `PNETCDF_VERSION` — parallel I/O only
- remove the MPI-linked netCDF-C and HDF5 builds; the serial stack becomes the
  only one. It can move from `/usr/local/serial` to `/usr/local`, which is where
  `ccs_config/machines/container/config_machines.xml` hardcodes `NETCDF_PATH`.
  Check what `PNETCDF_PATH=/usr/local` does once PnetCDF is gone.
- remove two of the three ESMF flavors; keep only `ESMF_COMM=mpiuni`
- drop the netCDF-C twin assertion — there is no second stack to agree with

`cime-macros/gnu_container.cmake` and the `ESMFMKFILE` path move with the
prefix change. The Dockerfile's existing macro-path assertion must be updated in
the same commit, not after.

Same file, same commit: add `MPI_SERIAL_PATH` and `PIO_LIBDIR`/`PIO_INCDIR`
pointing at the image's installs. Hardcode them the way `PFUNIT_PATH` and
`ESMFMKFILE` already are — derecho reads them from `$ENV{NCAR_ROOT_MPI_SERIAL}`
and `$ENV{PIO}`, which the container has no modules to set. **Extend the
Dockerfile's macro-path assertion to cover them**, so a wrong path fails the
image build rather than surfacing as a link error in someone's case, or — worse
— as a silent fall-back to compiling the submodule if the variable ends up empty.

`smoke-test.sh`: drop the MPI hello-world and the mpich version check.

### Phase 2 — record what the binary can't carry

netCDF-Fortran's version and the compression-filter set cannot be read from the
container binary. `smoke-test.sh` already reads versions from image LABELs set
at build time; extend that mechanism so the Dockerfile also records
`nf-config --version` and `nc-config --has-stdfilters`. The scheduled check then
reads a manifest rather than needing to run things inside the image.

### Phase 3 — rewrite the checker

`docker/ctsm-ci-derecho-gnu/check-derecho-versions.py`:

- **remove** the netcdf-mpi / parallel-netcdf / esmf-mpi checks, the cray-mpich
  deviation guard, and the `<MPILIBS>` assertion
- **re-aim** every remaining module check at the `mpilib="mpi-serial"` blocks
- **re-measure** `[snapshot]` against the serial `netcdf` module and the plain
  `hdf5` it pulls in via `depends_on`, not `netcdf-mpi`/`hdf5-mpi`. Move
  `SNAPSHOT_SOURCE_MODULE` and the provenance key with it.
- **add** the cray-libsci deviation guard
- **add** an existence check: every module `config_machines.xml` names for the
  gnu + mpi-serial path must resolve in derecho's module tree
- **keep** the existing `direct` comparison for `MPI_SERIAL_VERSION` and
  `PIO_VERSION` against the `mpilib="mpi-serial"` modules. No fxtag mode is
  needed: once `gnu_container.cmake` sets the paths, those ARGs describe what a
  case build actually links. Replace the comment explaining that they don't
  reach a case build — it will no longer be true.
- rewrite the module docstring around the new subject

### Phase 4 — rebuild, revalidate, republish

Same shape as NEXT_STEPS item 6. The image changes substantially, so this is a
full rebuild on Casper, re-run of all four wrapper scripts, republish to GHCR
and a new pin in `cirrus-testing.yml`.

### Phase 5 — CI and workflows

- `cirrus-testing.yml`: add `_Mmpi-serial` to the test names. They currently
  carry no `_M` modifier, so they take the machine default — which is mpich.
- Consider an upstream `ccs_config` change so
  `machines/container/config_machines.xml` declares `mpi-serial` in `<MPILIBS>`.
  mpi-serial works there today despite only `mpich` being listed (CIME needs no
  launcher for it), but the file would be lying. Shared with CESM, so it's a
  separate PR.
- `derecho-version-check.yml`: keep as the PR gate, drop
  `push: branches: ['**']`, add the staleness guard below.
- new `.github/workflows/derecho-drift.yml`: `schedule` + `workflow_dispatch`,
  `runs-on: gha-runner-ctsm`, resolve the latest master tag, read derecho's
  stack from `/glade`, open an issue with `GITHUB_TOKEN`.

**Staleness guard on the PR gate.** A `pull_request` event checks out the merge
commit, so the post-merge state is already evaluated; the gap is that the check
doesn't re-run when the base moves. That's a **stale green, not a wrong answer**
— the message must say so. Fail only when the base moved in a way this check can
see (`ccs_config`, `Dockerfile`, `derecho-versions.ini`, the checker):

- `git diff --name-only <merge-base> origin/<base.ref>` filtered to the same
  path list the workflow declares in `paths:`, so guard and trigger stay defined
  in one place
- `fetch-depth: 0` on `actions/checkout`
- compare against `github.event.pull_request.head.sha`, **not** `HEAD` — on the
  merge ref, base is always an ancestor of `HEAD`, so the naive
  `merge-base --is-ancestor` form silently always passes
- read `github.event.pull_request.base.ref`; don't hardcode `master`

**Live reads need no modules.** Use the method proven in
`probe-derecho-modules.yml`: resolve `default` symlinks, take the HDF5 pairing
from `netcdf`'s `depends_on`, run `nf-config` by absolute path. Loading
derecho's stack off-derecho is impossible (`cray-mpich` needs Cray PE) and
unnecessary.

### Phase 6 — docs

`README.md` deviation table and `NEXT_STEPS.md`: two exceptions, not three;
mpich gone; the phase vocabulary retired.

---

## Verification

- **The mpi-serial case is the acceptance test, and it has two halves.** The
  checker must still report `MPI_SERIAL_VERSION` matching derecho's
  `mpi-serial/2.5.3` — but that is only meaningful once a freshly built
  `cesm.exe` no longer contains `/ctsm/libraries/mpi-serial` in its
  compiled-source strings. Check the binary first; a green check without that is
  the old vacuous behaviour unchanged. Same for PIO.
- Mutation-test each check: perturb the ARG, the module version, the snapshot
  value and a pinned-module name in turn; each failure must be reported
  distinctly rather than collapsing or passing.
- Exercise the staleness guard both ways — a base that moved an irrelevant file
  must pass, one that moved `ccs_config` must fail — and confirm it isn't the
  always-passing `HEAD` form.
- After the rebuild, re-run `smoke-test.sh`, `smoke-test-pfunit.sh`,
  `run-unit-tests-in-container.sh` and `run-test-in-container.sh` on Casper.
- Build a current `_Mmpi-serial` `derecho_gnu` case so the Tier 3 comparison has
  a derecho-side binary from the same CTSM version as the container's.
- Run the scheduled workflow by `workflow_dispatch` from the branch; the
  schedule itself can only be confirmed after merge to master.

## Unrelated, found along the way

podman's runtime state defaults onto `$TMPDIR`, which in an interactive Casper
session is glade — shared across nodes and never cleared. A `pause.pid` written
Aug 30 on another node broke every podman command, including `system migrate`
and `system reset`, which fail before dispatching because the rootless re-exec
happens first. A stale-but-*dead* PID is handled fine, so the trigger was PID
reuse.

`build-on-casper.sh:83` already forces a node-local `TMPDIR` and so was never
affected — which is why this only appeared in hand-run podman. The fix is to set
`XDG_RUNTIME_DIR` explicitly, node-local, beside the graphroot: in `.bashrc`
guarded on unset, and in the wrapper scripts, which run where `.bashrc` is
silent. Container plumbing, not drift-checker work.
