# run_sys_tests in the ctsm-ci-derecho-gnu container

**Date:** 2026-08-31

> **Corrections applied after the final review.** Two things this document
> names were changed during implementation and are NOT accurate below:
>
> - The `MACHINE_DEFAULTS` key is **`ctsm-ci-container`**, not `container`.
>   `ccs_config/machines/container/config_machines.xml` already defines an
>   unrelated generic CIME machine with `MACH="container"`. The testid prefix
>   is therefore `ct`, not `co`.
> - `--wait` waits for **every** launched `create_test`, not just the last
>   one, and exits nonzero if any of them failed. Waiting only for the last
>   would let a container tear down its PID namespace and kill the tests still
>   running -- the exact failure `--wait` exists to prevent.

## Problem

`docker/ctsm-ci-derecho-gnu/` has three wrappers that run CTSM inside the
container: `run-case-in-container.sh`, `run-test-in-container.sh` and
`run-unit-tests-in-container.sh`. They cover `create_newcase`, `create_test`
and the pFUnit unit tests. Missing is CTSM's own test driver,
`run_sys_tests`, which is what a CTSM developer actually types.

`run_sys_tests` adds bookkeeping that the bare `create_test` wrapper does not
have: a dated testroot, `cs.status` / `cs.status.fails` aggregation, a
recorded `SRCROOT_GIT_STATUS`, baseline compare/generate, and retry. Those are
the reason to want it in the container, not the launching -- the container is
`BATCH_SYSTEM=none`, so there is nothing to submit to.

Four things block it today.

1. **`container` is not in `MACHINE_DEFAULTS`** (`python/ctsm/machine_defaults.py`).
   `get_machine_name()` is just the hostname, so the lookup misses and
   `create_machine` logs "machine %s not recognized; using generic no-batch
   settings", leaving `scratch_dir` and `baseline_dir` as `None`. That makes
   `--testroot-base` mandatory and `--compare`/`--generate` unusable.

2. **`run_sys_tests` never passes `--machine` to `create_test`.** It passes
   only `--xml-machine`, for the test-list query. CIME would therefore try to
   auto-detect the machine inside the container -- the exact failure the other
   three wrappers avoid by injecting `--machine container`.

3. **`testlist_clm.xml` has no `container` entries** (verified: zero). Suite
   mode resolves the suite through `--xml-machine`, which `run_sys_tests`
   hardcodes to the machine name, so `-s <suite>` would raise `No tests found
   for suite X on machine container`.

4. **The no-batch launcher backgrounds `create_test` and does not wait.**
   `JobLauncherNoBatch.run_command_impl` calls `subprocess.Popen` with a
   `preexec_fn` that ignores `SIGHUP`, and `run_sys_tests` returns without
   waiting. In a container the wrapper then exits, podman tears down the PID
   namespace, and `create_test` dies with it.

## Decisions

Each was chosen against alternatives; the alternatives are recorded so they
are not revisited.

### Topology: `run_sys_tests` runs inside the container

A fourth wrapper, `run-sys-tests-in-container.sh`, in the same shape as the
existing three.

*Rejected:* running `run_sys_tests` on the Casper/Derecho host and having each
`create_test` become a `qsub`ed job that invokes podman. That keeps
`run_sys_tests` in its designed role as a launcher, but needs a new job-launcher
class under `python/ctsm/joblauncher/` -- the qsub launcher builds its command
with no hook to wrap it in `podman run` -- and it can never work in GitHub
Actions, which has no PBS. Wiring runs into CI is `NEXT_STEPS.md` step 2.

### Suite lookup: split "what machine am I" from "whose testlist do I read"

Add an `--xml-machine` option to `run_sys_tests`, defaulting to the machine
name (today's behavior). The wrapper then passes `--machine-name container`
with `--xml-machine derecho`, so the container runs derecho's tests -- which is
the image's whole purpose, since it replicates derecho's gnu stack.

The motivating case is `aux_clm_mpi_serial`: 47 entries, of which 20 are
`derecho`/`gnu`. The rest -- 16 `derecho`/`intel` and 11 `izumi`/`nag` -- the
container cannot run, and the two injected flags exclude them independently:
`--xml-machine derecho` drops izumi, `--suite-compiler gnu` drops intel. That
category is why mpi-serial support was on the critical path for the image.

*Rejected:* passing `--machine-name derecho` and overriding the rest on the
command line. Zero code change, but `run_sys_tests` would log "Running on
machine: derecho" when it is not, testids would start with `de`, and we would
inherit derecho's `MACHINE_DEFAULTS` (account required, qsub launcher,
`/glade/derecho/scratch`) only to override every field.

*Rejected:* adding `<machine name="container" .../>` entries to
`testlist_clm.xml`. The most CTSM-idiomatic answer and it would give CI a named
suite, but it edits a heavily shared file needing upstream buy-in, and adds a
machine name meaningful only inside this container.

### Teardown: a `--wait` option on `run_sys_tests`

`--wait` calls `machine.job_launcher.wait_for_last_process_to_complete()` after
launching and exits with `create_test`'s status. A failing test therefore makes
`run_sys_tests` exit nonzero, since `create_test`'s code is propagated verbatim.

`wait_for_last_process_to_complete()` already exists on `JobLauncherNoBatch`,
is public, and is documented. It is called only by its own unit tests in
`python/ctsm/test/joblauncher/test_unit_job_launcher_no_batch.py`; no
production code path calls it. Those tests also already exercise the
no-process case, which is why returning `0` when nothing was launched is safe.

*Rejected:* solving teardown entirely in the wrapper. Bash cannot `wait` on the
orphan -- `wait` consults bash's own jobs table and a reparented process was
never a bash job -- so the wrapper would have to poll `/proc`:

```bash
pid=$(pgrep -f "scripts/create_test .*--test-id ${testid}" | head -1)
while [ -n "${pid}" ] && kill -0 "${pid}" 2>/dev/null; do sleep 15; done
```

That works (podman gives the container its own PID namespace, so `pgrep` sees
only our processes, and there is no startup race because the process exists
before `run_sys_tests` returns), but it discards the exit code. Recovering
pass/fail would mean parsing `cs.status.fails` output -- inventing a contract
against CIME's status formatting -- and it conflates "three tests FAILed" with
"create_test was OOM-killed halfway". For CI, a nonzero exit is the entire
interface.

### Baselines: container-local on scratch

`baseline_dir = /scratch/baselines`, i.e. host `$SCRATCH/cases_devcontainer/baselines`.
`--generate` then `--compare` works out of the box and stays writable.

Comparing against derecho's real gnu baselines remains one flag away, because
`/glade/campaign` is already mounted read-only for the inputdata-symlink fix:
`--baseline-root /glade/campaign/cgd/tss/ctsm_baselines`.

*Rejected:* defaulting to derecho's baselines. That mount is read-only so
`--generate` would fail against it, and container-vs-derecho answers will
differ anyway, making failure the norm rather than the signal.

## Design

### Flow

```
host: run-sys-tests-in-container.sh [run_sys_tests args]
  |
  |  ctsm_host_setup          (container-common.sh: mounts, image tag)
  |  ctsm_print_host_dirs     (case/build/run/archive, host-side paths)
  v
podman run  -v repo:/ctsm  -v inputdata:/opt/cesmdata/inputdata:ro
            -v $HOME/cases_devcontainer:/cases
            -v $SCRATCH/cases_devcontainer:/scratch
            -v /glade/campaign:/glade/campaign:ro
  |
  |  ctsm_container_setup     (.cime macro from the mounted repo, git identity)
  v
/ctsm/run_sys_tests --machine-name container --wait \
                    --extra-create-test-args "--machine container --compiler gnu" \
                    [--xml-machine derecho --suite-compiler gnu]  [user args]
  |
  |  MACHINE_DEFAULTS["container"] -> no-batch launcher, /scratch, no account
  |  testroot = /scratch/tests_<testid_base>
  |  cs.status + cs.status.fails written there
  v
create_test (foreground under --wait) -> exit code -> podman -> wrapper -> qsub
```

### Components

**`python/ctsm/machine_defaults.py`** -- new entry:

```python
"container": MachineDefaults(
    job_launcher_type=JOB_LAUNCHER_NOBATCH,
    scratch_dir=os.path.join(os.path.sep, "scratch"),
    baseline_dir=os.path.join(os.path.sep, "scratch", "baselines"),
    account_required=False,
    create_test_retry=0,
    create_test_queue=CREATE_TEST_QUEUE_UNSPECIFIED,
    job_launcher_defaults={},
),
```

Unlike every other entry these paths are not `get_user()`-derived: inside the
container the mounts are always at `/scratch` and `/cases` regardless of who
ran it. That asymmetry gets a comment in the file.

**`python/ctsm/run_sys_tests.py`** -- two options:

- `--xml-machine` (default `None`, falling back to `machine.name`), plumbed
  through `run_sys_tests()` to `_run_test_suite()`, into the `--xml-machine`
  element of `test_args`, and into `_get_compilers_for_suite()` so its
  "No tests found" message names the machine actually queried.
- `--wait`, plumbed to a call after launching, whose result becomes
  `sys.exit(rc)` in `main()`. `main()` currently never calls `sys.exit`.

**`python/ctsm/joblauncher/job_launcher_no_batch.py`** --
`wait_for_last_process_to_complete()` returns `self._process.wait()` instead of
discarding it, and `0` when nothing was launched.

**`python/ctsm/joblauncher/job_launcher_base.py`** -- the same method raising a
clear "`--wait` is only supported with the no-batch job launcher" error, so
`--wait` with qsub fails loudly instead of silently doing nothing.

**`docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh`** -- new wrapper.
PBS headers, sources `container-common.sh`, full argument passthrough.
Injects only what the user did not supply, via `ctsm_has_arg`:

| Injected | When | Why |
|---|---|---|
| `--machine-name container` | always | `MACHINE_DEFAULTS` lookup, testid prefix, cs.status |
| `--wait` | always | keep the container alive; propagate exit status |
| `--extra-create-test-args "--machine container --compiler gnu"` | always, merged not clobbered | `run_sys_tests` never passes `--machine` itself |
| `--xml-machine derecho` | only with `-s`/`--suite-name` | whose testlist to read; `XML_MACHINE` env overrides |
| `--suite-compiler gnu` | only with `-s`/`--suite-name` | only compiler in the image |

`--extra-create-test-args` is a single string, so a user-supplied value must be
concatenated with the injected one, not replaced.

### Error handling

Reuse the missing-inputdata scrape already in `run-test-in-container.sh`
(`Could not find all inputdata on any server|wget failed|Read-only file
system|Errno 30`) and the same host-side remediation message.

Beyond that, `--wait` means `create_test`'s exit code reaches `qsub`, so
failures are visible without scraping.

### Two accepted behaviors, documented not fixed

**`--suite-compiler gnu` skips the ctsm_pylib check.** `_check_py_env` is
called from inside `_get_compilers_for_suite`, which only runs when no compiler
was given. Injecting `--suite-compiler gnu` therefore bypasses it. The
consequence is specific and acceptable: of the 20 gnu tests in
`aux_clm_mpi_serial`, exactly one (`FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon`)
needs python modules the image lacks, and it will now fail on its own rather
than aborting all 20 before any run. On the `-t` path `_check_py_env` *is*
unconditional, so naming that test explicitly still raises `ModuleNotFoundError`
up front.

**The convenience symlink dangles from the host.** `_make_testroot` links the
testroot into the working directory (`/cases`), so
`$HOME/cases_devcontainer/tests_<id>` points at `/scratch/tests_<id>` -- valid
inside the container, broken outside it. The same class of thing as the
`DOUT_S_ROOT` redirect in `run-case-in-container.sh`. Fixing it means either
`cd`ing somewhere that collides with the link name or patching CTSM.

## Testing

**Unit tests**, run in `ctsm_pylib`. `python/ctsm/test/test_unit_run_sys_tests.py`
already drives command construction through a fake job launcher;
`test_createTestCommands_testsuiteSpecifiedCompilers` is the model for:

- `--xml-machine` absent -> `--xml-machine` equals the machine name (no regression)
- `--xml-machine derecho` with `--machine-name container` -> the derecho testlist
  is queried while the machine stays `container`
- `test_unit_machine.py`: the `container` entry yields a no-batch launcher,
  `/scratch`, `/scratch/baselines`, and no account requirement

**On Casper, by `qsub`**, in order:

1. `-t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc` -- already known-good
   through `run-test-in-container.sh`, so it isolates the `run_sys_tests` layer.
2. `-s aux_clm_mpi_serial --dry-run` -- confirms suite resolution against
   derecho's testlist without running anything.
3. A real short suite, confirming `cs.status.fails`, the testroot layout, and a
   propagated nonzero exit when a test fails.

**No image rebuild is required.** `run_sys_tests`, `python/ctsm` and `cime` all
come from the repo mounted at `/ctsm`; the image supplies only the compiled
stack and the python interpreter. This work therefore does not trigger the
rebuild/revalidate/republish/tag-bump cycle.

## Out of scope

- Adding `container` entries to `testlist_clm.xml`.
- A host-side podman job launcher.
- Wiring any of this into `cirrus-testing.yml` (`NEXT_STEPS.md` step 2, which
  still needs the `/glade`-visibility probe on `gha-runner-ctsm`).
- ~~Waiting for more than one launched process; the wrapper always forces a
  single compiler.~~ **REVERSED after the final whole-change review
  (2026-08-31).** The stated reason does not hold: the wrapper injects
  `--suite-compiler gnu` only when the user has not passed their own, so a
  single compiler is a convention, not an enforcement -- and
  `--job-launcher-nobatch` reaches the same code path on any machine.

  The consequence is worse than "the other exit statuses are discarded",
  which is how this was first characterized. `JobLauncherNoBatch.run_command_impl`
  reassigns `self._process` on every call without waiting, so
  `_run_test_suite`'s per-compiler loop leaves N `create_test` processes
  running concurrently and only the last-launched one is ever waited on. On a
  login node the others survive, because `preexec_fn` ignores SIGHUP and they
  are reparented. **In a container they do not**: `run_sys_tests` exits, podman
  tears down the PID namespace, and the still-running `create_test` processes
  are killed mid-test -- exactly the failure `--wait` exists to prevent.

  So `--wait` now waits for **every** launched process and returns nonzero if
  any of them failed. The launcher tracks a list of processes rather than one,
  and the method is named for that. See the post-review fix commits.
