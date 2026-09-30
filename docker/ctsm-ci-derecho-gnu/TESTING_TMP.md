# Testing the cirrus-runner-workflows branch

**Validate `run-sys-tests-in-container.sh` on Casper.**

Written, syntax checked, and dry-run-verified on the login node only (see "Added 2026-08-31: run_sys_tests wrapper" in `NEXT_STEPS.md`) -- it has not yet been run in full against the actual container on a compute node. Each section below is ordered by what it actually proves; do them in order and do not skip one because a later one looks like it would cover it too.

Every section has the same shape: **What this tests**, then **Run it**, then **Check**, then **Result**. Each `Run it` block is complete on its own, including its `execcasper` and `podman load`, so sections 2-6 can be run in one session or in separate ones. A section whose header carries a ✅ is finished; one without it is either unrun or only partly done, and its **Result** says which.

## 1. Suite resolution via `--xml-machine` ✅

**What this tests:** that `--xml-machine` reaches the testlist query and resolves a suite's compilers. This is the only way to exercise `_get_compilers_for_suite`, and so `get_tests_from_xml`, through `run_sys_tests`' own code: the wrapper always injects `--suite-compiler gnu` (there is no way to suppress it), which skips that call entirely. It runs on the host, from `ctsm_pylib`, with no container and no allocation.

An earlier version of this section ran `query_testlists.py` instead. That was near-useless: it calls `get_tests_from_xml` directly as CIME's own script, so it exercised CIME and the testlist data rather than any of this branch's code.

**Run it:**

```bash
module load conda
conda activate ctsm_pylib
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows

./run_sys_tests --machine-name ctsm-ci-container --xml-machine derecho \
    -s aux_clm_mpi_serial --skip-compare --skip-generate \
    --dry-run --skip-git-status -v
```

**Check:** the compilers resolve from derecho's testlist, and one `create_test` command is built per compiler, each carrying `--xml-machine derecho`, under a testroot named `tests_<MMDD-HHMMSS>ct`.

**Result: passed, 2026-09-30.** Resolved `['gnu', 'intel']` and built two `create_test` commands with testids `<MMDD-HHMMSS>ct_gnu` and `..._int` under testroot `/scratch/tests_<MMDD-HHMMSS>ct`. One host-vs-container difference this surfaced: on the host, `create_machine` finds an account and adds `--project`; inside the container there is none, so `--project` is absent there (confirmed in section 2).

## 2. `run_sys_tests` starting up inside the container ✅

**What this tests:** that `run_sys_tests` imports and runs under the container's `python3`; that `create_machine("ctsm-ci-container")` resolves; that `git` / `bin/git-fleximod status` succeed against the bind-mounted `/ctsm`; and that the testroot is named as predicted (`tests_<MMDD-HHMMSS>ct`).

It does **not** prove suite resolution -- the wrapper always injects `--suite-compiler gnu`, which makes `run_sys_tests` skip `_get_compilers_for_suite`, the only caller of `get_tests_from_xml`, and `--dry-run` stops `create_test` from running at all. That is section 1's job.

**Run it:**

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=01:00:00

# inside the session; podman's storage is node-local and does not survive it.
# The tarball restores as localhost/ctsm-ci-derecho-gnu:dev, the wrappers'
# default IMAGE_TAG, so nothing needs re-tagging.
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows

docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial --dry-run -v --skip-compare --skip-generate
echo "wrapper exit status: $?"
```

`-v` is load-bearing: `run_sys_tests` logs the assembled `create_test` command at INFO, so without it you see only the `Testroot:` line and cannot check what would have been run.

**Check:** exit status 0. The wrapper's `host-side directories:` block names a testroot `tests_<MMDD-HHMMSS>ct`. The logged `Running: <.../create_test ...>` line carries `--xml-category aux_clm_mpi_serial --xml-machine derecho --xml-compiler gnu`, `--output-root /scratch/tests_<MMDD-HHMMSS>ct` and `--baseline-root /scratch/baselines` (that pair is `MACHINE_DEFAULTS["ctsm-ci-container"]` resolving), and `--machine container --compiler gnu` from the injected `--extra-create-test-args`. Unlike the same dry run on the host, there should be **no** `--project`, since the container has no account. Nothing should have been created under `$SCRATCH/cases_devcontainer`.

**Result: passed, 2026-09-30, exit 0.** The assembled command was:

```
/ctsm/cime/scripts/create_test --test-id 0930-220320ct_gnu --output-root /scratch/tests_0930-220320ct --xml-category aux_clm_mpi_serial --xml-machine derecho --xml-compiler gnu --baseline-root /scratch/baselines --retry 0 --machine container --compiler gnu
```

Every item on the check list held: `--output-root` and `--baseline-root` are `MACHINE_DEFAULTS["ctsm-ci-container"]` resolving, `--machine container --compiler gnu` is the injected `--extra-create-test-args`, `_gnu` on the test id is `_NUM_COMPILER_CHARS = 3`, there is no `--project`, and nothing was created under `$SCRATCH/cases_devcontainer`.

## 3. `--wait` blocking, and the exit status reaching PBS ✅

**What this tests:** three things -- that `--wait` blocks rather than returning as soon as `create_test` is launched; that a nonzero `create_test` status survives the trip back out through podman to PBS; and that a *passing* test returns 0, so that a nonzero status means something. All three were untested outside unit tests against a fake launcher.

The test to use is `SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs`. **Not `SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc`**, which no longer builds: the ctsm5.4.054 rebase changed that alias in `cime_config/config_compsets.xml` from `..._MOSART_SGLC_SWAV` to `..._MOSART_DGLC%NOEVOLVE_SWAV`, and `components/cdeps/dglc/cime_config/buildnml:63-65` refuses single-point runs (`single column mode for DGLC is not currently allowed`). `1x1_brazil` is single point, so it fails in SHAREDLIB_BUILD. The cdeps guard is not new -- it is in `cdeps1.0.79` too -- so the compset change is what broke it. `IHistClm60BgcQianRsGs` (`HIST_DATM%QIA_CLM60%BGC_SICE_SOCN_SROF_SGLC_SWAV`) stays on SGLC and is what derecho's `aux_clm_mpi_serial` entry for this grid uses now; `SMS_D_Ld1` keeps the one-day debug shape of the original.

**Run it:**

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=04:00:00

# inside the session; podman's storage is node-local and does not survive it.
# The tarball restores as localhost/ctsm-ci-derecho-gnu:dev, the wrappers'
# default IMAGE_TAG, so nothing needs re-tagging.
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows

docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs --skip-compare --skip-generate
echo "wrapper exit status: $?"
```

**Check:** the status is 0. Then read the result from the host -- the generated `cs.status` in the testroot bakes in container paths and will not run here:

```bash
cime/CIME/Tools/cs.status $SCRATCH/cases_devcontainer/tests_<MMDD-HHMMSS>ct/*/TestStatus
```

Confirm the case PASSes. While it runs, progress is in `$SCRATCH/cases_devcontainer/tests_<MMDD-HHMMSS>ct/STDOUT.<MMDD-HHMMSS>ct` and the matching `STDERR.*`, not in the PBS log. The wrapper prints the exact path as `testroot` in its `host-side directories:` block.

**Result: passed, 2026-09-30.** All three settled, across two runs.

A run with the old `IHistClm60Bgc` test name exited **100** -- `create_test`'s own status for a failed test, returned by `--wait`, carried back out through podman and reported by the wrapper. That settled blocking (a status of 100 could only arrive after `create_test` finished) and nonzero propagation. The test failed in SHAREDLIB_BUILD for the DGLC reason above, unrelated to the wrapper.

A second run, with `IHistClm60BgcQianRsGs`, exited **0** and every phase PASSed under testroot `tests_0930-220608ct`:

```
PASS ... SHAREDLIB_BUILD time=81
PASS ... MODEL_BUILD time=15
PASS ... RUN time=48
PASS ... MEMLEAK insufficient data for memleak test
PASS ... SHORT_TERM_ARCHIVER
```

Together those give what CI needs: the wrapper exits nonzero *if and only if* a test failed. The pair also rules out a wrapper hardwired to one status, which either run alone would have been consistent with. Incidentally, SHAREDLIB_BUILD passing is the direct check that `IHistClm60BgcQianRsGs` really is on SGLC, so the section 3 substitution is sound; the whole test took about two and a half minutes, so the four-hour walltime above is far more than this section needs on its own.

## 4. A nonzero exit status when a test fails ✅

**What this tests:** that a test failing *inside* `create_test`, after `run_sys_tests` has launched it, produces a nonzero exit status. A failure that `run_sys_tests` catches first never reaches `--wait` and so proves nothing about the propagation path -- which rules out the obvious candidate, `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon`: in the `-t` branch `_check_py_env` runs before `_run_create_test` (`python/ctsm/run_sys_tests.py:273`) and aborts on any name containing `FSURDATMODIFYCTSM`, so nothing is launched. Use a name `create_test` itself will reject instead.

The `FSURDATMODIFYCTSM` run is still worth doing once, as a check of that early-abort path rather than of `--wait` (missing python modules; see the `-s` / ctsm_pylib note in README "Running run_sys_tests").

**Run it:**

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=01:00:00

# inside the session; podman's storage is node-local and does not survive it.
# The tarball restores as localhost/ctsm-ci-derecho-gnu:dev, the wrappers'
# default IMAGE_TAG, so nothing needs re-tagging.
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows

docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcNOSUCHCOMPSET --skip-compare --skip-generate
echo "create_test-rejects exit status: $?"

docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon --skip-compare --skip-generate
echo "early-abort exit status: $?"
```

**Check:** both statuses are nonzero, and for different reasons. For the first, `$SCRATCH/cases_devcontainer/tests_<MMDD-HHMMSS>ct/STDOUT.<MMDD-HHMMSS>ct` (the wrapper prints the exact path) shows `create_test` rejecting the compset -- its own stdout, not STDERR, which carries only CIME's python-version warning. For the second, expect `ModuleNotFoundError: modify_fsurdat can't be loaded`; since nothing is launched, the status comes from the uncaught exception rather than from `create_test`, so it should be 1 rather than 100. The testroot is still created -- `_make_testroot` and `_record_git_status` run before `_check_py_env` -- so it holds `SRCROOT_GIT_STATUS` and `cs.status.fails`; the tell that nothing launched is the absence of a case directory and of any `STDOUT.`/`STDERR.` files.

If the second returns **2**, check the command: `run_sys_tests` requires one of `-c`/`--skip-compare` and one of `-g`/`--skip-generate`, and argparse exits 2 when they are missing, before `_check_py_env` is ever reached.

**Result: passed, 2026-09-30.** Both halves, and they fail by different mechanisms as intended.

The bogus-compset run exited **100** with `create_test` launched and failing in CREATE_NEWCASE under testroot `tests_0930-221656ct`:

```
ERROR: Invalid compset name, IHistClm60BgcNOSUCHCOMPSET, all stub components generated
FAIL SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcNOSUCHCOMPSET.container_gnu (phase CREATE_NEWCASE)
```

That is a genuine post-launch failure reaching the wrapper through `--wait`, confirming section 3's nonzero half by a second route.

The early-abort run exited **1**, on exactly the predicted path:

```
File "/ctsm/python/ctsm/run_sys_tests.py", line 273, in run_sys_tests
    _check_py_env(testname_list)
File "/ctsm/python/ctsm/run_sys_tests.py", line 806, in _check_py_env
    raise ModuleNotFoundError("modify_fsurdat" + err_msg) from err
ModuleNotFoundError: modify_fsurdat can't be loaded. Do you need to activate the ctsm_pylib conda environment?
```

Its testroot `tests_0930-221901ct` holds only `SRCROOT_GIT_STATUS` and `cs.status.fails` -- no case directory, no `STDOUT.`/`STDERR.` -- confirming `create_test` never ran. The two statuses being 100 and 1 is itself the evidence that they are distinct paths rather than one generic failure.

An earlier attempt at this second run returned 2, because the command here was missing `--skip-compare --skip-generate` and argparse rejected it before `_check_py_env`.

## 5. A full suite end to end

**What this tests:** the whole suite path, `-s aux_clm_mpi_serial`, including `cs.status` generation and the multi-test testroot layout.

Judge this run by the suite's failures, **not** by the wrapper's exit code: a nonzero exit is *expected* on a first full run, for reasons unrelated to this change -- `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` needs python modules the image lacks (same as section 4's second run), and the suite's NEON / `CLM_USRDAT` and FATES entries need user datasets and FATES build support that this change does not touch.

Do **not** use `-s clm_short` as a substitute "quick suite" -- it has exactly two derecho/gnu entries, `ERP_D_P64x2_Ld3.f10_f10_mg37.I1850Clm50BgcCrop` and `ERS_D_Ld3.f10_f10_mg37.I1850Clm50BgcCrop`, neither mpi-serial, and `P64x2` wants 64 MPI tasks against a machine config with `MAX_MPITASKS_PER_NODE=4` inside an 8-cpu PBS reservation; it will fail for reasons that have nothing to do with this change.

**Run it:**

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=12:00:00

# inside the session; podman's storage is node-local and does not survive it.
# The tarball restores as localhost/ctsm-ci-derecho-gnu:dev, the wrappers'
# default IMAGE_TAG, so nothing needs re-tagging.
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows

docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial --skip-compare --skip-generate
echo "wrapper exit status: $?"
```

The 12-hour walltime above matches the wrapper's own PBS header; this section is hours, not minutes. `qsub`-ing the script instead does not work: it takes no arguments that way, so it would run with neither `-s` nor `-t`.

**Check:** ignore the exit status, for the reasons above, and judge by the failures:

```bash
cime/CIME/Tools/cs.status --fails-only $SCRATCH/cases_devcontainer/tests_<MMDD-HHMMSS>ct/*/TestStatus
```

Expect `FSURDATMODIFYCTSM_D_Mmpi-serial_Ld1.5x5_amazon` plus the NEON / `CLM_USRDAT` and FATES entries to appear there; anything else is a finding.

**Result: not yet run.**

## 6. The replacement test, in the other two wrappers

**What this tests:** that `IHistClm60BgcQianRsGs` is a working substitute everywhere `IHistClm60Bgc` was used. It replaced it in `run-test-in-container.sh`'s `default_test` and in the `run-case-in-container.sh` example, for the DGLC reason in section 3, and neither has been run since. Running `run-test-in-container.sh` bare is the path that has been broken since the rebase.

These share a session with sections 2 and 3, so the block below does all four -- one `podman load` covers everything, and section 2 costs seconds.

**Run it** -- fastest first, so a failure stops you before the slow ones:

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=04:00:00

# inside the session; podman's storage is node-local and does not survive it.
# The tarball restores as localhost/ctsm-ci-derecho-gnu:dev, the wrappers'
# default IMAGE_TAG, so nothing needs re-tagging.
module load podman
export TMPDIR=/var/tmp/$USER      # rootless podman needs node-local scratch
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows

# section 2: costs seconds, nothing is built
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -s aux_clm_mpi_serial --dry-run -v --skip-compare --skip-generate
echo "section 2 exit status: $?"

# section 3: the exit-0 path, with the replacement test
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs \
    --skip-compare --skip-generate
echo "section 3 exit status: $?"

# the same test through the create_test wrapper, via its new default_test
docker/ctsm-ci-derecho-gnu/run-test-in-container.sh
echo "run-test exit status: $?"

# the create_newcase example from README "Running cases and tests"
docker/ctsm-ci-derecho-gnu/run-case-in-container.sh \
    --case brazil_test --compset IHistClm60BgcQianRsGs --res 1x1_brazil \
    --mpilib mpi-serial --run-unsupported
echo "run-case exit status: $?"
```

If `run-case-in-container.sh` is re-run, give it a fresh `--case` name or remove `$HOME/cases_devcontainer/brazil_test` first; `create_newcase` will not overwrite an existing case directory.

**Check:** all four exit 0. The failure to watch for is the one section 3 describes -- `single column mode for DGLC is not currently allowed` in SHAREDLIB_BUILD -- which would mean `IHistClm60BgcQianRsGs` is not on SGLC after all and the replacement is wrong. Anything else is a fault in that wrapper rather than in the test choice.

**Result: partly done, 2026-09-30 -- two of the four.** The first two commands in that block are sections 2 and 3, both of which passed, and section 3's PASS through SHAREDLIB_BUILD is the evidence that `IHistClm60BgcQianRsGs` is a working substitute. What remains is the other two wrappers: `run-test-in-container.sh` bare, which reaches the same test through `create_test` rather than `run_sys_tests`, and the `run-case-in-container.sh` example.

## What to watch for

Highest-risk failure modes first:
- `_record_git_status` can abort `run_sys_tests` up front if git's dubious-ownership / `safe.directory` check trips on the bind-mounted repo. Mitigation: pass `--skip-git-status`.
- The silent job log described in README "Running run_sys_tests" (Test output does not appear in the job log): with `--wait`, nothing streams to the PBS log between "Running: <create_test ...>" and the final exit code, so a long quiet job is expected, not a hang -- watch `<testroot>/STDOUT.<testid>` / `STDERR.<testid>` instead.
- `MAX_MPITASKS_PER_NODE=4` together with `GMAKE_J=4` means CIME may build and run up to 4 tests at once inside the 8-cpu PBS reservation these wrappers request -- expect concurrency, not one test running at a time.
- The `podman load` OOM (exit 137, prints nothing, `podman tag`/`podman run` then fail with the misleading "image not known") documented in `NEXT_STEPS.md` under "Resolved 2026-08-31: mpi-serial runs work" applies here too: load the image from a session with real memory before running any of the above.
