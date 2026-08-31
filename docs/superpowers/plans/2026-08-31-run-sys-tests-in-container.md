# run_sys_tests in the container — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make CTSM's `run_sys_tests` usable inside the `ctsm-ci-derecho-gnu` container, joining the three wrappers that already cover `create_newcase`, `create_test` and the pFUnit unit tests.

**Architecture:** Three small, independently useful additions to CTSM's python — a `container` entry in `MACHINE_DEFAULTS`, an `--xml-machine` option that separates "what machine am I" from "whose testlist do I read", and an opt-in `--wait` that makes `run_sys_tests` block and propagate `create_test`'s exit status — plus a fourth wrapper script that supplies the container plumbing and injects those flags.

**Tech Stack:** Python 3 (CTSM's `python/ctsm`, tested with `unittest` and formatted with `black`), bash, podman on NCAR HPC (Casper), PBS.

**Spec:** `docs/superpowers/specs/2026-08-31-run-sys-tests-in-container-design.md`

## Global Constraints

- **Python environment: `ctsm_pylib`.** All python unit tests run there. Do NOT create, modify, or survey any other environment; if `ctsm_pylib` is unavailable, stop and ask.
- **No image rebuild.** `run_sys_tests`, `python/ctsm` and `cime` all come from the repo mounted at `/ctsm`. Nothing here triggers the rebuild → revalidate → republish → tag-bump cycle. Do not edit `Dockerfile` or the tag in `.github/workflows/cirrus-testing.yml`.
- **Container-absolute paths, fixed by the wrapper's mounts:** repo `/ctsm`, cases `/cases`, scratch `/scratch`, baselines `/scratch/baselines`, inputdata `/opt/cesmdata/inputdata` (read-only), symlink root `/glade/campaign` (read-only).
- **Host-side equivalents:** `$HOME/cases_devcontainer` ↔ `/cases`, `$SCRATCH/cases_devcontainer` ↔ `/scratch`.
- **Formatting:** CTSM enforces `black`. After any python edit, run from `python/`: `black --config pyproject.toml ctsm` and verify with `black --check --config pyproject.toml .`
- **No `ccs_config` changes.** That submodule stays untouched; this is a settled project decision.
- **No `testlist_clm.xml` changes.** Out of scope per the spec.
- **Backward compatibility:** `--xml-machine` and `--wait` both default to today's behavior. No existing `run_sys_tests` invocation may change behavior.

---

### Task 1: `container` machine defaults

Adds the `MACHINE_DEFAULTS` entry whose absence makes `create_machine` log "machine %s not recognized; using generic no-batch settings" and leave `scratch_dir`/`baseline_dir` as `None`.

**Files:**
- Modify: `python/ctsm/machine_defaults.py`
- Test: `python/ctsm/test/test_unit_machine.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `MACHINE_DEFAULTS["container"]`, a `MachineDefaults` namedtuple with `job_launcher_type=JOB_LAUNCHER_NOBATCH`, `scratch_dir="/scratch"`, `baseline_dir="/scratch/baselines"`, `account_required=False`, `create_test_retry=0`, `create_test_queue=CREATE_TEST_QUEUE_UNSPECIFIED`, `job_launcher_defaults={}`.

- [x] **Step 1: Write the failing test**

In `python/ctsm/test/test_unit_machine.py`, add `mock` to the imports at the top of the file (it is not currently imported):

```python
from unittest import mock
```

Then add this method to `class TestCreateMachine`, after `test_unknownMachine_defaults`:

```python
    def test_containerMachine_defaults(self):
        """Tests the container machine

        Unlike the other machines, its paths are container-absolute rather than
        user-specific, because the container's mounts are always at the same
        place regardless of who ran it. It also requires no account, so
        create_machine must succeed with account=None; _get_account is mocked
        out so this doesn't depend on the environment the test runs in.
        """
        with mock.patch("ctsm.machine._get_account", return_value=None):
            machine = create_machine("container", MACHINE_DEFAULTS)
        self.assertMachineInfo(
            machine=machine,
            name="container",
            scratch_dir=os.path.join(os.path.sep, "scratch"),
            baseline_dir=os.path.join(os.path.sep, "scratch", "baselines"),
            account=None,
        )
        self.assertNoBatchInfo(machine)
```

- [x] **Step 2: Run test to verify it fails**

Run:
```bash
cd python && python -m unittest ctsm.test.test_unit_machine.TestCreateMachine.test_containerMachine_defaults -v
```
Expected: FAIL with `AssertionError: None != '/scratch'`. (`defaults.get("container")` returns `None`, so `create_machine` takes its unknown-machine branch rather than raising.)

- [x] **Step 3: Write minimal implementation**

In `python/ctsm/machine_defaults.py`, extend the import at the top — it currently imports only `JOB_LAUNCHER_QSUB`:

```python
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_QSUB, JOB_LAUNCHER_NOBATCH
```

Then add this entry to the `MACHINE_DEFAULTS` dict, after the `"cheyenne"` entry so the dict stays roughly alphabetical:

```python
    "container": MachineDefaults(
        # The ctsm-ci-derecho-gnu container (docker/ctsm-ci-derecho-gnu/).
        # NOTE: unlike every other machine here, these paths are NOT
        # get_user()-derived. Inside the container the mounts are always at
        # /scratch and /cases no matter who ran it; the host-side directories
        # they map to are chosen by the wrapper scripts.
        job_launcher_type=JOB_LAUNCHER_NOBATCH,
        scratch_dir=os.path.join(os.path.sep, "scratch"),
        baseline_dir=os.path.join(os.path.sep, "scratch", "baselines"),
        account_required=False,
        create_test_retry=0,
        create_test_queue=CREATE_TEST_QUEUE_UNSPECIFIED,
        # The container has BATCH_SYSTEM=none, so there is no queue to
        # configure and no non-default launcher worth having defaults for.
        job_launcher_defaults={},
    ),
```

- [x] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m unittest ctsm.test.test_unit_machine -v
```
Expected: all PASS, including the eight pre-existing tests.

- [x] **Step 5: Format and commit**

```bash
cd python && black --config pyproject.toml ctsm && cd ..
git add python/ctsm/machine_defaults.py python/ctsm/test/test_unit_machine.py
git commit -m "Add container to MACHINE_DEFAULTS

Without this, run_sys_tests falls through to its unknown-machine branch:
scratch_dir and baseline_dir come back None, making --testroot-base
mandatory and --compare/--generate unusable.

Unlike the other entries the paths are container-absolute rather than
get_user()-derived, because the container's mounts are always at the same
place regardless of who ran it."
```

---

### Task 2: `--xml-machine` — read another machine's testlist

`run_sys_tests` hardcodes `--xml-machine <machine.name>`, so suite mode inside the container would query a `container` testlist that does not exist (verified: zero `container` entries in `testlist_clm.xml`). This splits the knob so the container can run derecho's tests without claiming to be derecho.

**Files:**
- Modify: `python/ctsm/run_sys_tests.py`
- Test: `python/ctsm/test/test_unit_run_sys_tests.py`

**Interfaces:**
- Consumes: `MACHINE_DEFAULTS["container"]` from Task 1 (not directly exercised here).
- Produces: `run_sys_tests(..., xml_machine=None)` keyword argument, and the `--xml-machine` command-line option (`args.xml_machine`). When `None`, behavior is identical to today.

- [x] **Step 1: Write the failing test**

Add to `class TestRunSysTests` in `python/ctsm/test/test_unit_run_sys_tests.py`, after `test_createTestCommands_testsuiteSpecifiedCompilers`:

```python
    def test_createTestCommands_testsuiteXmlMachine(self):
        """With xml_machine given, the suite is queried for THAT machine

        The machine we are actually running on is unchanged: the testid still
        carries this machine's prefix. This is what lets the container run
        derecho's testlist without claiming to be derecho.
        """
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.datetime") as mock_date, mock.patch(
            "ctsm.run_sys_tests.get_tests_from_xml"
        ) as mock_get_tests:
            mock_date.now.side_effect = self._fake_now
            mock_get_tests.return_value = [{"compiler": "gnu"}]
            run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                suite_name="my_suite",
                xml_machine="derecho",
            )

        mock_get_tests.assert_called_once_with(xml_machine="derecho", xml_category="my_suite")
        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        self.assertRegex(all_commands[0].cmd, r"--xml-machine +derecho(\s|$)")
        self.assertNotRegex(all_commands[0].cmd, r"--xml-machine +{}(\s|$)".format(self._MACHINE_NAME))
        # The machine we are running on is unaffected: testid prefix is still "fa"
        self.assertRegex(
            all_commands[0].cmd,
            r"--test-id +{}_gnu(\s|$)".format(self._expected_testid()),
        )
```

- [x] **Step 2: Run test to verify it fails**

```bash
cd python && python -m unittest ctsm.test.test_unit_run_sys_tests.TestRunSysTests.test_createTestCommands_testsuiteXmlMachine -v
```
Expected: FAIL with `TypeError: run_sys_tests() got an unexpected keyword argument 'xml_machine'`.

- [x] **Step 3: Write the implementation**

Three edits in `python/ctsm/run_sys_tests.py`.

**(a)** In `_commandline_args()`, next to the existing `--machine-name` argument, add:

```python
    parser.add_argument(
        "--xml-machine",
        default=None,
        help="Machine name to use when looking up which tests are in a suite.\n"
        "This is separate from --machine-name, which says what machine we are\n"
        "actually running on. Separating them lets one machine run another\n"
        "machine's tests - e.g. a container that replicates derecho running\n"
        "derecho's suites, since the container has no testlist entries of its\n"
        "own.\n"
        "Only used together with --suite-name.\n"
        "Default: the value of --machine-name.",
    )
```

**(b)** In `main()`, add to the `run_sys_tests(...)` call, after `suite_compilers=args.suite_compiler,`:

```python
        xml_machine=args.xml_machine,
```

**(c)** In `run_sys_tests()`, add `xml_machine=None` to the keyword-only arguments, document it in the docstring alongside `suite_compilers`:

```python
    xml_machine (str or None): machine name to use when querying which tests are in
        the given suite (only used with suite_name); if not specified, use
        machine.name
```

and pass it through in the `_run_test_suite(...)` call:

```python
            xml_machine=xml_machine,
```

Then in `_run_test_suite()`, add `xml_machine` to its keyword-only parameters and replace the top of its body:

```python
    xml_machine_final = xml_machine if xml_machine else machine.name
    if not suite_compilers:
        suite_compilers = _get_compilers_for_suite(
            suite_name, xml_machine_final, running_ctsm_py_tests
        )
    for compiler in suite_compilers:
        test_args = [
            "--xml-category",
            suite_name,
            "--xml-machine",
            xml_machine_final,
            "--xml-compiler",
            compiler,
        ]
```

Note `_get_compilers_for_suite` already takes the machine name as its second parameter and uses it both for the `get_tests_from_xml` query and for its "No tests found for suite {} on machine {}" error, so passing `xml_machine_final` makes that error name the machine actually queried. No change is needed inside that function.

- [x] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m unittest ctsm.test.test_unit_run_sys_tests -v
```
Expected: all PASS. `test_createTestCommands_testsuite` is the no-regression check — it asserts `--xml-machine fake_machine` when `xml_machine` is not given.

- [x] **Step 5: Format and commit**

```bash
cd python && black --config pyproject.toml ctsm && cd ..
git add python/ctsm/run_sys_tests.py python/ctsm/test/test_unit_run_sys_tests.py
git commit -m "run_sys_tests: add --xml-machine

Separates 'what machine am I' from 'whose testlist do I read'. Previously
--xml-machine was hardcoded to the machine name, so a machine with no
testlist_clm.xml entries could not use suite mode at all.

The motivating case is the ctsm-ci-derecho-gnu container, which replicates
derecho's gnu stack and so wants to run derecho's suites -- notably
aux_clm_mpi_serial -- without pretending to be derecho.

Defaults to the machine name, so no existing invocation changes."
```

---

### Task 3: `--wait` — block and propagate the exit status

`JobLauncherNoBatch.run_command_impl` `Popen`s `create_test` and `run_sys_tests` returns without waiting. In a container the wrapper then exits, podman tears down the PID namespace, and `create_test` dies with it.

**Files:**
- Modify: `python/ctsm/joblauncher/job_launcher_base.py`
- Modify: `python/ctsm/joblauncher/job_launcher_no_batch.py`
- Modify: `python/ctsm/joblauncher/job_launcher_fake.py`
- Modify: `python/ctsm/run_sys_tests.py`
- Test: `python/ctsm/test/test_unit_run_sys_tests.py`

**Interfaces:**
- Consumes: `run_sys_tests(...)` from Task 2.
- Produces: `run_sys_tests(..., wait=False)` keyword argument; `run_sys_tests()` now **returns an int** (0 unless `wait=True`, in which case `create_test`'s return code); the `--wait` command-line option; `JobLauncherBase.wait_for_last_process_to_complete()` raising `RuntimeError`; `JobLauncherFake.set_return_code(int)`.

- [x] **Step 1: Write the failing tests**

Add to `class TestRunSysTests` in `python/ctsm/test/test_unit_run_sys_tests.py`:

```python
    def test_wait_returnsCreateTestStatus(self):
        """With wait=True, run_sys_tests returns create_test's exit status

        A failing test therefore surfaces as a nonzero exit, which is what CI
        needs and what makes the container wrapper report failures.
        """
        machine = self._make_machine()
        machine.job_launcher.set_return_code(3)
        return_code = run_sys_tests(
            machine=machine, cime_path=self._cime_path(), testlist=["foo"], wait=True
        )
        self.assertEqual(return_code, 3)

    def test_noWait_returnsZero(self):
        """Without wait, run_sys_tests returns 0 even if create_test would fail

        This is the pre-existing behavior and must not change: run_sys_tests
        normally returns as soon as the job is launched.
        """
        machine = self._make_machine()
        machine.job_launcher.set_return_code(3)
        return_code = run_sys_tests(
            machine=machine, cime_path=self._cime_path(), testlist=["foo"]
        )
        self.assertEqual(return_code, 0)

    def test_wait_withDryRun_returnsZero(self):
        """With dry_run there is no process to wait for, so the status is 0"""
        machine = self._make_machine()
        machine.job_launcher.set_return_code(3)
        return_code = run_sys_tests(
            machine=machine,
            cime_path=self._cime_path(),
            testlist=["foo"],
            wait=True,
            dry_run=True,
        )
        self.assertEqual(return_code, 0)
```

- [x] **Step 2: Run tests to verify they fail**

```bash
cd python && python -m unittest ctsm.test.test_unit_run_sys_tests.TestRunSysTests.test_wait_returnsCreateTestStatus -v
```
Expected: FAIL with `AttributeError: 'JobLauncherFake' object has no attribute 'set_return_code'`.

- [x] **Step 3: Write the implementation**

**(a)** In `python/ctsm/joblauncher/job_launcher_base.py`, add this method to `JobLauncherBase`, after `run_command_logger_message`:

```python
    def wait_for_last_process_to_complete(self):
        """Wait for the last process started by run_command to complete

        Returns that process's return code.

        Only launchers that run the job synchronously under this process can
        support this. For launchers whose job outlives us (e.g. qsub) there is
        nothing to wait for, so this raises.
        """
        raise RuntimeError(
            "Waiting for the launched job is only supported with the no-batch "
            "job launcher; this is a {}".format(type(self).__name__)
        )
```

**(b)** In `python/ctsm/joblauncher/job_launcher_no_batch.py`, replace the existing `wait_for_last_process_to_complete` with:

```python
    def wait_for_last_process_to_complete(self):
        """Waits for the last process started by run_command_impl (if any) to complete

        Returns that process's return code, or 0 if no process was started -
        which is the case after a dry run.
        """
        if self._process is not None:
            return self._process.wait()
        return 0
```

**(c)** In `python/ctsm/joblauncher/job_launcher_fake.py`, add a canned return code so the plumbing is testable. Replace `__init__` and add two methods:

```python
    def __init__(self):
        JobLauncherBase.__init__(self)
        self._commands = []
        self._return_code = 0

    def set_return_code(self, return_code):
        """Set the value returned by wait_for_last_process_to_complete"""
        self._return_code = return_code

    def wait_for_last_process_to_complete(self):
        """Pretend to wait; returns the canned return code set by set_return_code"""
        return self._return_code
```

**(d)** In `python/ctsm/run_sys_tests.py`, add the option in `_commandline_args()`, next to `--dry-run`:

```python
    parser.add_argument(
        "--wait",
        action="store_true",
        help="Wait for the launched create_test to finish, and exit with its\n"
        "status, rather than returning as soon as it has been launched.\n"
        "This makes a failing test produce a nonzero exit status.\n"
        "Only supported with the no-batch job launcher.\n"
        "Note that this waits for the LAST create_test launched, so if the\n"
        "suite spans multiple compilers it returns when the final one is done.\n"
        "Needed when running inside a container, which would otherwise be torn\n"
        "down while create_test was still running.",
    )
```

Add `wait=False` to the keyword-only arguments of `run_sys_tests()`, and document it:

```python
    wait (bool): if True, wait for the launched create_test to complete and return its
        exit status (rather than returning as soon as it is launched)
```

Add to the end of `run_sys_tests()`, after the `if suite_name: ... else: ...` block, at the same indentation as that `if`:

```python
    if wait and not dry_run:
        return machine.job_launcher.wait_for_last_process_to_complete()
    return 0
```

Finally in `main()`, pass the flag and exit with the result. Replace the bare `run_sys_tests(` call with an assignment, add `wait=args.wait,` to its arguments, and follow it with the exit:

```python
    return_code = run_sys_tests(
        machine=machine,
        cime_path=cime_path,
        ...
        xml_machine=args.xml_machine,
        wait=args.wait,
    )
    sys.exit(return_code)
```

`sys` is already imported at the top of the module.

- [x] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m unittest ctsm.test.test_unit_run_sys_tests ctsm.test.joblauncher.test_unit_job_launcher_no_batch -v
```
Expected: all PASS. The two pre-existing no-batch launcher tests already call `wait_for_last_process_to_complete()`, including once after a dry run, so they cover the `return 0` path.

- [x] **Step 5: Run the whole unit suite for regressions**

```bash
cd python && ./run_ctsm_py_tests --unit
```
Expected: no new failures against the pre-change baseline. If any test failed before your changes, note it and confirm it is unrelated rather than assuming.

- [x] **Step 6: Format and commit**

```bash
cd python && black --config pyproject.toml ctsm && cd ..
git add python/ctsm/run_sys_tests.py python/ctsm/test/test_unit_run_sys_tests.py \
        python/ctsm/joblauncher/job_launcher_base.py \
        python/ctsm/joblauncher/job_launcher_no_batch.py \
        python/ctsm/joblauncher/job_launcher_fake.py
git commit -m "run_sys_tests: add --wait

The no-batch launcher Popens create_test and run_sys_tests returns without
waiting. That is right on a login node, where the job should outlive the
shell, but fatal in a container: the wrapper exits, podman tears down the
PID namespace, and create_test dies with it.

--wait blocks on the launched process and exits with its status, so a
failing test surfaces as a nonzero exit. wait_for_last_process_to_complete()
already existed on the no-batch launcher with no production caller; it now
returns the process's code instead of discarding it.

Opt-in, so no existing invocation changes."
```

---

### Task 4: The wrapper script and docs

**Files:**
- Create: `docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh`
- Modify: `docker/ctsm-ci-derecho-gnu/README.md`
- Modify: `docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md`

**Interfaces:**
- Consumes: `MACHINE_DEFAULTS["container"]` (Task 1), `--xml-machine` (Task 2), `--wait` (Task 3). From `container-common.sh`: `ctsm_host_setup`, `ctsm_has_arg <flag> "$@"`, `ctsm_podman_run <cmd...>`, `ctsm_container_setup`, and the variables `CTSM_CASES_DIR`, `CTSM_SCRATCH_DIR`, `CTSM_INPUTDATA`.
- Produces: the executable wrapper.

- [x] **Step 1: Write the wrapper**

Create `docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh`:

```bash
#!/usr/bin/env bash
#PBS -N ctsm-devcontainer-sys-tests
#PBS -q casper
#PBS -l select=1:ncpus=8:mem=96GB
#PBS -l walltime=12:00:00
#PBS -j oe
#
# Run CTSM's run_sys_tests inside the ctsm-ci-derecho-gnu container.
#
# A thin wrapper around ./run_sys_tests: every argument is passed straight
# through to it, and this script only adds the container plumbing (mounts, the
# $HOME/.cime macro copy) plus the flags that tell run_sys_tests it is in a
# container.
#
# Why this rather than run-test-in-container.sh, which already wraps
# create_test: run_sys_tests adds the bookkeeping -- a dated testroot,
# cs.status / cs.status.fails aggregation, a recorded SRCROOT_GIT_STATUS,
# baseline compare/generate, and retry. It is also what a CTSM developer
# actually types.
#
# The larger PBS reservation than the other wrappers is deliberate: a suite
# runs its tests serially in one container (BATCH_SYSTEM=none), and podman
# load alone peaks around 54-56 GB because the vfs driver duplicates layers.
#
# Usage:
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh [run_sys_tests args]
# Examples:
#   # one test by name
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
#       -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc
#   # derecho's mpi-serial suite, filtered to this image's one compiler
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial
#   # see what would run, without running it
#   docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial --dry-run
#
# These defaults are injected only if you did not pass them yourself:
#   --machine-name container  picks up MACHINE_DEFAULTS["container"]: the
#                             no-batch launcher, /scratch as the testroot base,
#                             /scratch/baselines, and no account requirement
#   --wait                    run_sys_tests otherwise backgrounds create_test and
#                             returns, and podman would tear the container down
#                             while it was still running
#   --extra-create-test-args "--machine container --compiler gnu"
#                             run_sys_tests never passes --machine to
#                             create_test itself, only --xml-machine, so
#                             without this CIME tries to auto-detect the
#                             machine from the container hostname and fails.
#                             A value you pass is MERGED with this, not
#                             replaced.
# and, only when you asked for a suite with -s/--suite-name:
#   --xml-machine derecho     whose testlist to read. The container has no
#                             entries in cime_config/testdefs/testlist_clm.xml,
#                             and replicating derecho is the image's purpose.
#                             Override with $XML_MACHINE.
#   --suite-compiler gnu      the only compiler in the image
#
# Overridable via environment (see container-common.sh for the full list):
#   IMAGE_TAG  INPUTDATA  CASES_DIR  SCRATCH_DIR  MAKE_J  MOUNT_OPTS  DRY_RUN
#   XML_MACHINE  the machine whose testlist a suite is read from
set -eo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=container-common.sh
source "${here}/container-common.sh"

set -u

ctsm_host_setup

# run_sys_tests takes --extra-create-test-args as a SINGLE string, which it
# later splits. So a value the user passed has to be merged with ours rather
# than either one clobbering the other. Pull it out of the argument list here
# and re-add the combined value below.
injected_extra="--machine container --compiler gnu"
user_extra=""
args=()
skip_next=0
for arg in "$@"; do
    if [ "${skip_next}" -eq 1 ]; then
        user_extra="${arg}"
        skip_next=0
        continue
    fi
    case "${arg}" in
        --extra-create-test-args=*) user_extra="${arg#--extra-create-test-args=}" ;;
        --extra-create-test-args)   skip_next=1 ;;
        *)                          args+=("${arg}") ;;
    esac
done

cime_args=()
ctsm_has_arg --machine-name "${args[@]+"${args[@]}"}" \
    || cime_args+=(--machine-name container)
ctsm_has_arg --wait "${args[@]+"${args[@]}"}" \
    || cime_args+=(--wait)
cime_args+=(--extra-create-test-args "${injected_extra}${user_extra:+ ${user_extra}}")

# Suite-only injections. --xml-machine and --suite-compiler are meaningless
# without -s, and passing them anyway would just be confusing output.
if ctsm_has_arg -s "${args[@]+"${args[@]}"}" \
   || ctsm_has_arg --suite-name "${args[@]+"${args[@]}"}"; then
    ctsm_has_arg --xml-machine "${args[@]+"${args[@]}"}" \
        || cime_args+=(--xml-machine "${XML_MACHINE:-derecho}")
    ctsm_has_arg --suite-compiler "${args[@]+"${args[@]}"}" \
        || ctsm_has_arg --suite-compilers "${args[@]+"${args[@]}"}" \
        || cime_args+=(--suite-compiler gnu)
fi

cime_args+=("${args[@]+"${args[@]}"}")

# run_sys_tests names the testroot tests_<MMDD-HHMMSS><first 2 chars of the
# machine name>, so "co" here. Unlike the other wrappers there is no separate
# --test-root / --output-root: run_sys_tests passes --output-root <testroot> to
# create_test, so case dirs, bld/ and run/ all land inside the testroot.
echo "host-side directories:"
echo "  testroot   ${CTSM_SCRATCH_DIR}/tests_<MMDD-HHMMSS>co"
echo "             (cases, bld/ and run/ all land inside it)"
echo "  baselines  ${CTSM_SCRATCH_DIR}/baselines"
echo "note: run_sys_tests also links the testroot into ${CTSM_CASES_DIR},"
echo "      but that link points at the CONTAINER path and so is dangling"
echo "      when read from the host. Use the testroot path above."

log="${CTSM_SCRATCH_DIR}/run_sys_tests.$(date +%Y%m%d%H%M%S).log"
echo "log        ${log}"
echo "running    run_sys_tests ${cime_args[*]}"

set +e
ctsm_podman_run bash -s -- "${cime_args[@]}" <<'INSIDE' 2>&1 | tee "${log}"
set -eu

source /ctsm/docker/ctsm-ci-derecho-gnu/container-common.sh
ctsm_container_setup

echo "HOME=$HOME"

# run_sys_tests makes a symlink to the testroot in the working directory.
cd /cases

exec /ctsm/run_sys_tests "$@"
INSIDE
rc=${PIPESTATUS[0]}
set -e

# Same failure mode as run-test-in-container.sh: CIME cannot fetch missing
# inputdata through a read-only mount, and the first suspect is dangling
# symlinks rather than genuinely absent data. See that script for the full
# explanation of the wording matched here.
if [ "${rc}" -ne 0 ] \
   && grep -qE 'Could not find all inputdata on any server|wget failed|Read-only file system|Errno 30' "${log}"; then
    cat >&2 <<MSG

=====================================================================
ERROR: CIME could not find input data and failed trying to download it.

  host:      ${CTSM_INPUTDATA}
  container: /opt/cesmdata/inputdata (read-only)

FIRST SUSPECT: dangling symlinks, not genuinely missing data. Much of
the inputdata tree is absolute symlinks into sibling campaign
collections, which resolve only if that tree is mounted too. Check the
"symlinks" line in this script's output: if it is absent, set
INPUTDATA_SYMLINK_ROOT (default /glade/campaign) to the tree the
symlinks point into.

Full log: ${log}
=====================================================================
MSG
fi

exit "${rc}"
```

- [x] **Step 2: Verify the script parses and is executable**

```bash
chmod +x docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh
bash -n docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh && echo "syntax OK"
```
Expected: `syntax OK`.

- [x] **Step 3: Verify argument assembly without a container**

`DRY_RUN` is honored by `ctsm_podman_run`, so this exercises the injection logic on the host. Run all four and check the `running run_sys_tests ...` line:

```bash
cd /glade/work/samrabin/ctsm_cirrus-runner-workflows
DRY_RUN=1 docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60Bgc 2>&1 | grep '^running'
DRY_RUN=1 docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s aux_clm_mpi_serial 2>&1 | grep '^running'
DRY_RUN=1 docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh \
    -s aux_clm_mpi_serial --extra-create-test-args "--debug" 2>&1 | grep '^running'
DRY_RUN=1 XML_MACHINE=izumi docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -s fates 2>&1 | grep '^running'
```

Expected, in order:
1. `--machine-name container --wait --extra-create-test-args --machine container --compiler gnu -t SMS_...` — note **no** `--xml-machine` or `--suite-compiler`, because this is not a suite.
2. the same plus `--xml-machine derecho --suite-compiler gnu -s aux_clm_mpi_serial`.
3. the extra args **merged**: `--extra-create-test-args --machine container --compiler gnu --debug` — `--debug` present and the injected pair still there.
4. `--xml-machine izumi`, showing `XML_MACHINE` overrides.

If test 3 shows only `--debug`, or drops it, the merge loop is wrong — fix before continuing.

- [x] **Step 4: Commit the wrapper**

```bash
git add docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh
git commit -m "Add run-sys-tests-in-container.sh

The fourth wrapper, alongside run-case / run-test / run-unit-tests. It adds
what run_sys_tests gives over the bare create_test wrapper: a dated testroot,
cs.status aggregation, recorded git status, baselines and retry.

Injects --machine-name container, --wait, and the --machine/--compiler that
run_sys_tests never passes to create_test itself; and for suites,
--xml-machine derecho (the container has no testlist entries of its own) and
--suite-compiler gnu."
```

- [ ] **Step 5: Validate on Casper, in order** -- DEFERRED TO THE USER.
  Needs an interactive compute-node session; not done as part of this plan.
  Tracked in `docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md` under Remaining steps.

The wrapper reads `"$@"`, like the other three, so it takes no arguments
through `qsub -v`. Validate from an **interactive** session, which is also how
the other wrappers were validated:

```bash
execcasper -A <PROJECT> -l select=1:ncpus=8:mem=96GB -l walltime=04:00:00
module load podman
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_20260831.tar
```

If `podman load` is killed with exit 137, it was OOM-killed and the next
command will misleadingly say `image not known`; ask for more memory.

**The canonical five-step sequence lives in
`docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md`, "Remaining steps" item 3** --
follow it there rather than duplicating it here, since it also records what
each step does and does *not* prove, and the highest-risk failure modes.
In outline: (1) `query_testlists.py` on the host for real suite resolution;
(2) the wrapper with `-s aux_clm_mpi_serial --dry-run`; (3) one known-good
test with `-t`, checking `echo $?`; (4) a deliberately failing test, checking
`echo $?` is nonzero; (5) the full `-s aux_clm_mpi_serial`, judged by
`cs.status.fails` rather than the exit code.

Two validation steps written into an earlier draft of this plan were wrong
and were corrected there:

- **`-s ... --dry-run` does not prove suite resolution.** The wrapper always
  injects `--suite-compiler gnu`, which makes `run_sys_tests` skip
  `_get_compilers_for_suite` -- the only caller of `get_tests_from_xml`. Only
  `query_testlists.py` exercises that path.
- **`-s clm_short` is not a usable "quick suite" here.** It has exactly two
  `derecho`/`gnu` entries, neither mpi-serial, and one of them (`P64x2`)
  wants 64 MPI tasks against `MAX_MPITASKS_PER_NODE=4`.

Confirm across the sequence: the testroot appears under
`$SCRATCH/cases_devcontainer/` named `tests_<MMDD-HHMMSS>ct`, `cs.status.fails`
exists inside it, and a run containing a failing test exits nonzero.

- [x] **Step 6: Document in the container README**

In `docker/ctsm-ci-derecho-gnu/README.md`, extend the "Running cases and tests" section with a `run_sys_tests` subsection covering: the two example invocations from the script header; the injected-defaults table; that `-s` reads **derecho's** testlist via `--xml-machine`, overridable with `XML_MACHINE`; that the testroot holds cases, `bld/` and `run/` together, unlike the other wrappers; that the `/cases` symlink dangles when read from the host; and the ctsm_pylib caveat below.

Add this note verbatim, since it is the one surprising behavior:

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

- [x] **Step 7: Record status in NEXT_STEPS.md**

Add an "Added 2026-08-31: run_sys_tests wrapper" section recording: the four CTSM-side changes and why each was needed; that `container` now exists in `MACHINE_DEFAULTS`; that `--xml-machine` and `--wait` are upstreamable on their own merits and default to prior behavior; and that **no image rebuild was needed**, because everything runs from the repo mounted at `/ctsm`.

Also update the "Remaining steps" list: this work does not close any of the four items, but step 2 ("Wire runs into CI") should now note that `run_sys_tests --wait` gives CI the nonzero exit status it needs.

- [x] **Step 8: Commit the docs**

```bash
git add docker/ctsm-ci-derecho-gnu/README.md docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md
git commit -m "Document run_sys_tests in the container

Covers the injected defaults, that -s reads derecho's testlist, the
testroot layout, the host-dangling symlink, and the one surprise: injecting
--suite-compiler skips run_sys_tests' ctsm_pylib check."
```

---

## Notes for the implementer

- **Do not create or survey python environments.** Use `ctsm_pylib` and nothing else. If it is missing or broken, stop and ask rather than investigating alternatives.
- **Tasks 1–3 are testable without a container or a compute node.** Only Task 4 Step 5 needs Casper.
- **If a `run_ctsm_py_tests --unit` failure predates your change**, say so explicitly with the evidence rather than fixing it silently or assuming it is unrelated.
- **`git status` should be clean between tasks.** Each task commits its own files and nothing else.
