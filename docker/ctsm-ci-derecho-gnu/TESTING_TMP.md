# Testing the cirrus-runner-workflows branch

**Validate `run-sys-tests-in-container.sh` on Casper.**

Written, syntax checked, and dry-run-verified on the login node only (see "Added 2026-08-31: run_sys_tests wrapper" in `NEXT_STEPS.md`) -- it has not yet been run against the actual container on a compute node. Each step below is ordered by what it actually proves; do them in order and do not skip one because a later step looks like it would cover it too.

## 1. Suite resolution via `--xml-machine`

**Done 2026-09-30; recorded because nothing else covers this path.**
   `run_sys_tests` itself, on the host, from `ctsm_pylib` -- no container and no allocation -- with `--xml-machine` and deliberately *no* `--suite-compiler`:
   ```
   ./run_sys_tests --machine-name ctsm-ci-container --xml-machine derecho \
       -s aux_clm_mpi_serial --skip-compare --skip-generate \
       --dry-run --skip-git-status -v
   ```
   This is the only way to reach `_get_compilers_for_suite`, and so `get_tests_from_xml`, through `run_sys_tests`' own code: the wrapper always injects `--suite-compiler gnu` (there is no way to suppress it), which skips that call entirely. Result: resolves `['gnu', 'intel']` from derecho's testlist and builds two `create_test` commands carrying `--xml-machine derecho`, with testids `<MMDD-HHMMSS>ct_gnu` and `..._int` under testroot `/scratch/tests_<MMDD-HHMMSS>ct`.
   An earlier version of this step ran `query_testlists.py` instead. That was near-useless: it calls `get_tests_from_xml` directly as CIME's own script, so it exercised CIME and the testlist data rather than any of this branch's code.
   One host-vs-container difference this surfaces: on the host, `create_machine` finds an account and adds `--project`. Inside the container there is none, so `--project` is absent there.

## Setup for sections 2-5

Sections 2-5 run the wrapper against the real container, so they need a compute node with the image loaded. The wrapper reads `"$@"`, so it takes no arguments through `qsub -v`; run it from an **interactive** session, which is also how the other wrappers were validated. Give the session enough walltime for whichever section you are on -- section 5 is the long one.

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=04:00:00
```

Then once per session, because podman's storage is node-local and does not survive it:

```bash
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows
```

The tarball restores as `localhost/ctsm-ci-derecho-gnu:dev`, the wrappers' default `IMAGE_TAG`, so nothing needs re-tagging. Every command below runs from the repo root, and `<testroot>` below means the path the wrapper prints as `testroot` in its `host-side directories:` block.

## 2. `run_sys_tests` starting up inside the container

The wrapper with `-s aux_clm_mpi_serial --dry-run`. This does **not** prove suite resolution -- the wrapper always injects `--suite-compiler gnu`, which makes `run_sys_tests` skip `_get_compilers_for_suite`, the only caller of `get_tests_from_xml`, and `--dry-run` stops `create_test` from running at all. What it does prove: `run_sys_tests` imports and runs under the container's `python3`; `create_machine("ctsm-ci-container")` resolves; `git`/`bin/git-fleximod status` succeed against the bind-mounted `/ctsm`; and the testroot is named as predicted (`tests_<MMDD-HHMMSS>ct`).

**Run it:**

```bash
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial --dry-run -v
```

`-v` is load-bearing: `run_sys_tests` logs the assembled `create_test` command at INFO, so without it you see only the `Testroot:` line and cannot check what would have been run.

**Check:** exit status 0. The wrapper's `host-side directories:` block names a testroot `tests_<MMDD-HHMMSS>ct`. The logged `Running: <.../create_test ...>` line carries `--xml-category aux_clm_mpi_serial --xml-machine derecho --xml-compiler gnu`, `--output-root /scratch/tests_<MMDD-HHMMSS>ct` and `--baseline-root /scratch/baselines` (that pair is `MACHINE_DEFAULTS["ctsm-ci-container"]` resolving), and `--machine container --compiler gnu` from the injected `--extra-create-test-args`. Unlike the same dry run on the host, there should be **no** `--project`, since the container has no account. Nothing should have been created under `$SCRATCH/cases_devcontainer`.

## 3. `--wait` blocking, and the exit status reaching PBS

One known-good test through the wrapper -- `-t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc`, then `echo $?`. This is the **only** step that proves `--wait` actually blocks and that the exit status propagates through podman to PBS -- the entire reason the `--wait` work exists, and until this runs it is untested outside unit tests against a fake launcher. Confirm the testroot appears under `$SCRATCH/cases_devcontainer/` and the case PASSes.

**Run it:**

```bash
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc
echo "wrapper exit status: $?"
```

**Check:** the status is 0, *and* it is printed only after the test has finished. If it comes back in seconds, `--wait` is not blocking and the rest of this section proves nothing. Then read the result from the host -- the generated `cs.status` in the testroot bakes in container paths and will not run here:

```bash
cime/CIME/Tools/cs.status <testroot>/*/TestStatus
```

Confirm the case PASSes. While it runs, progress is in `<testroot>/STDOUT.<MMDD-HHMMSS>ct` and the matching `STDERR.*`, not in the PBS log.

## 4. A nonzero exit status when a test fails

A deliberately failing test, to confirm the exit status is nonzero. The test must fail *inside* `create_test`, after `run_sys_tests` has launched it; a failure that `run_sys_tests` catches first never reaches `--wait` and so proves nothing about the propagation path.

**Run it** with a name `create_test` itself will reject. The obvious candidate, `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon`, is the wrong choice here: in the `-t` branch `_check_py_env` runs before `_run_create_test` (`python/ctsm/run_sys_tests.py:273`) and aborts on any name containing `FSURDATMODIFYCTSM`, so nothing is launched.

```bash
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcNOSUCHCOMPSET
echo "wrapper exit status: $?"
```

**Check:** the status is nonzero, and `<testroot>/STDERR.<MMDD-HHMMSS>ct` shows `create_test` rejecting the compset. This is the half of the `--wait` contract section 3 cannot show: that a nonzero `create_test` status survives the trip back through podman to PBS.

The `FSURDATMODIFYCTSM` run is still worth doing once, as a check of that early-abort path rather than of `--wait` (missing python modules; see the `-s` / ctsm_pylib note in README "Running run_sys_tests"):

```bash
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon
echo "wrapper exit status: $?"
```

Expect a nonzero status and `ModuleNotFoundError: modify_fsurdat can't be loaded` before any testroot contents appear.

## 5. A full suite end to end

The full `-s aux_clm_mpi_serial`. Judge this run by the suite's failures (command below), **not** by the wrapper's exit code: a nonzero exit is *expected* on a first full run of this suite, for reasons unrelated to this change -- `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` needs python modules the image lacks (same as section 4's second run), and the suite's NEON/`CLM_USRDAT` and FATES entries need user datasets and FATES build support that this change does not touch. Do **not** use `-s clm_short` as a substitute "quick suite" -- it has exactly two derecho/gnu entries, `ERP_D_P64x2_Ld3.f10_f10_mg37.I1850Clm50BgcCrop` and `ERS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop`, neither mpi-serial, and `P64x2` wants 64 MPI tasks against a machine config with `MAX_MPITASKS_PER_NODE=4` inside an 8-cpu PBS reservation; it will fail for reasons that have nothing to do with this change.

**Run it:**

```bash
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial
echo "wrapper exit status: $?"
```

Budget for hours, not minutes, and size the `execcasper` walltime accordingly -- the wrapper's own PBS header asks for `walltime=12:00:00`. `qsub`-ing the script directly does not work for this: it takes no arguments that way, so it would run with neither `-s` nor `-t`.

**Check:** ignore the exit status, for the reasons above, and judge by the failures:

```bash
cime/CIME/Tools/cs.status --fails-only <testroot>/*/TestStatus
```

Expect `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` plus the NEON / `CLM_USRDAT` and FATES entries to appear there; anything else is a finding.

## What to watch for

Highest-risk failure modes first:
- `_record_git_status` can abort `run_sys_tests` up front if git's dubious-ownership / `safe.directory` check trips on the bind-mounted repo. Mitigation: pass `--skip-git-status`.
- The silent job log described in README "Running run_sys_tests" (Test output does not appear in the job log): with `--wait`, nothing streams to the PBS log between "Running: <create_test ...>" and the final exit code, so a long quiet job is expected, not a hang -- watch `<testroot>/STDOUT.<testid>` / `STDERR.<testid>` instead.
- `MAX_MPITASKS_PER_NODE=4` together with `GMAKE_J=4` means CIME may build and run up to 4 tests at once inside the 8-cpu PBS reservation these wrappers request -- expect concurrency, not one test running at a time.
- The `podman load` OOM (exit 137, prints nothing, `podman tag`/`podman run` then fail with the misleading "image not known") documented in `NEXT_STEPS.md` under "Resolved 2026-08-31: mpi-serial runs work" applies here too: load the image from a session with real memory before running any of the above.