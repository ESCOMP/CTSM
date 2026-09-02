# Framework specifics

`SKILL.md` states the rules; this file gives the exact commands and the exact failure modes.
The pFUnit and CTSM rows were verified in this checkout. The pytest rows are the general shape
of the same rule and are worth confirming against your own suite before relying on them.

## The anti-retrofit check: which glob, and what it misses

| | Test files | Test-side code the glob does **not** match |
|---|---|---|
| pFUnit (CTSM) | `':(top)*.pf'` | `src/unit_test_shr/*.F90` — the shared fixture harnesses; a test directory's `CMakeLists.txt` |
| pytest | `':(top)test_*.py'`, `':(top)*_test.py'` | `conftest.py`, fixture and factory modules, test data files |

**The cwd trap, observed here:** run from inside one test directory, `git diff <base>..HEAD --
'*.pf'` reported on three changed files and came back clean about a fourth that had changed in a
sibling test directory. Quoting does not help — the shell was never the problem. Only `:(top)`,
`:/`, or running from the repository root does.

## Where a first assertion failure stops the test

| | Mechanism |
|---|---|
| pFUnit | The preprocessor emits `if (anyExceptions()) return` after every `@assert` it rewrites (`<PFUNIT_PATH>/PFUNIT-4.8/bin/funit/pFUnitParser.py:318`, and five sibling sites for the other assert forms), so the first failure returns from the test subroutine |
| pytest | A failing `assert` raises `AssertionError`, which ends that test function |

Either way, later assertions in the same test body never run — which is why a mutation has to be
narrowed until it reaches the assertion you are actually trying to exercise.

## What a test that cannot run costs you

| | Failure | Blast radius |
|---|---|---|
| pFUnit (CTSM) | A compile error in any one test directory | `run_tests.py` completes every build stage before it starts the first ctest stage (`cime/scripts/fortran_unit_testing/run_tests.py:440-492`), so **no test in the entire suite runs** |
| pytest | A collection error in a test module | That module's tests do not run; an error in `conftest.py`, or `-x`, can end the session |

## Confirming a new test actually ran

`SKILL.md` requires the writer to confirm its tests ran, not merely that the suite was green.
In CTSM this is not paranoia: adding a new `.pf` **file** to an existing test directory leaves
the generated pFUnit driver stale, and the run reports every test passing while executing the
previous set of tests. A red-first test can therefore appear to pass at its own test commit,
which destroys the evidence the commit exists to carry. The rebuild that fixes it, the narrow
boundary of the trap, and the per-directory executable whose reported count is the real
confirmation are in the **ctsm-unit-tests** skill.

## Mutations that abort instead of failing an assertion

CTSM's unit tests build under `CESM_DEBUG`, where `-check bounds` turns a broken index into an
abort and `-fpe0` turns a signalling NaN or a zero divisor into `forrtl` errors 65 and 73. So a
mutation you expected to fail an assertion may stop the binary instead, with no assertion
reported. That still counts as the test catching the mutation — say which form you saw. What
that build does and does not check is in the **ctsm-unit-tests** skill.
