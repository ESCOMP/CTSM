---
name: ctsm-unit-tests
description: Use when writing, building, running, or debugging a CTSM Fortran unit test - the pFUnit .pf files under src/*/test/ that run_tests.py builds into unit_tests.temp - or when a test build reports "Unrecognized token '@' skipped" or "Global name too long", or when a green run leaves you unsure whether your new or renamed test actually ran.
---

# CTSM unit tests

CTSM's Fortran unit tests are pFUnit 4 tests in `.pf` files under
`src/<area>/test/<Name>_test/`, built and run by `run_tests.py`.

## Start here

For the basics, read these rather than re-deriving them — everything below is what they
do not say.

- `src/README.unit_testing` — the command that runs the whole suite.
- `cime/scripts/fortran_unit_testing/README` — how the framework fits together.
- `doc/design/pfunit_testing.rst` — error handling and ESMF initialization in tests.

## Symptom → cause → fix

| Symptom | Cause | Fix |
|---|---|---|
| `error #5078: Unrecognized token '@' skipped`, then `#5082` on the same line | An `@assert` is split across lines with `&`, or carries a trailing comment | Put the whole `@assert` on one line; move the comment above it |
| Build is green and reports every test passing, but your new test never ran | You added a `.pf` **file** to an existing test directory; the driver object is stale | `rm -rf unit_tests.temp`, or delete only `<Name>_driver.F90.o` |
| `undefined reference to <module>_suite_` at link time | You *removed* a `.pf` from an existing test directory; same stale driver object | Same |
| You renamed a `@Test` subroutine and the run still reports the old name, or reports nothing different at all | Stale build | `rm -rf unit_tests.temp`; then confirm the executable's test count is what it was *before* the rename |
| `warning #5462: Global name too long, shortened from:` | `<module>_mp_<PROCEDURE>` exceeds 90 characters | Shorten the module name or the test name |
| `forrtl: error (65): floating invalid`, with no assertion reported | `-fpe0` trapped arithmetic on a signalling NaN — often a fixture component nothing ever set | Treat the abort as the failure and find the unset value; it is not a crash to debug |
| `forrtl: error (73): floating divide by zero` | `-fpe0` trapped a division a fixture left with a zero divisor | Same — the abort *is* the result, so read it as one |
| A compiler warning you expected is nowhere in the build log | `run_tests.py` prints compiler output only when the build fails | Run `make` yourself in the test's build directory |

## The incremental rebuild will lie to you

The dangerous case is a **green run that proves nothing**. `run_tests.py` reuses
`unit_tests.temp`, and a `.pf` file you have just added to an existing test directory is
compiled, linked, and then never called. The run exits 0 and prints
`100% tests passed, 0 tests failed out of <N>`. Your tests did not run.

Only the generated driver goes stale. pFUnit's driver does `#include _TEST_SUITES`
(`<PFUNIT_PATH>/PFUNIT-*/include/driver.F90`), and nothing tracks the generated
`<Name>.inc` as a dependency of `<Name>_driver.F90` — so the `.inc` is regenerated and
does list your new suite, while the object that would call it is not rebuilt. That cause
explains most of the boundary, which is narrower than it looks:

| Change to an existing, populated `unit_tests.temp` | Result |
|---|---|
| New `@Test` subroutine in an existing `.pf` | Picked up |
| **Renaming** a `@Test` subroutine in an existing `.pf` | **Needs a full rebuild** |
| New `.pf` file in an **existing** test directory, even added to `CMakeLists.txt` | **Silently dropped** |
| Removing a `.pf` from an existing test directory | **Link failure** |
| New test **directory** with its own `add_pfunit_ctest` | Picked up |

A new directory is safe because it has no stale driver yet. Editing
`CMakeLists.txt` does not rescue you.

The rename row is the one the stale-driver account does not predict — a rename changes
neither the suite list nor the suite's name, so on that mechanism it ought to behave like
adding a subroutine. It does not; it is recorded here from experience rather than from the
mechanism, so treat a rename as needing the wipe and do not reason from the cause above to
the contrary. Its check is different too, because the count is *supposed* to stay put: the
rename took if the executable still reports the same number of tests, and broke pFUnit's
discovery if the count fell.

The fix is `rm -rf unit_tests.temp` and a full rebuild. If you would rather not pay for
that, deleting the one stale object is enough:

```
unit_tests.temp/__command_line_test__/__command_line_test__/<lib>_test/<Dir>_test/CMakeFiles/<Name>.dir/<Name>_driver.F90.o
```

### Two things that are not evidence

**The reported test count counts executables, not tests.** The `out of N` in ctest's
summary is the number of ctest entries — one per `add_pfunit_ctest`, i.e. one per test
directory. It cannot move when you add a test, or a whole file, to a directory that
already exists, so watching it is not a check.

**A green build prints nothing.** `run_tests.py` captures cmake and make output and shows
it only on failure, so silence tells you nothing about warnings or about what compiled.

### What is evidence

Run the test directory's own executable and read its `(N tests)` line:

```
unit_tests.temp/__command_line_test__/__command_line_test__/<lib>_test/<Dir>_test/<Name>
```

Add `-v` to list the individual test names. Check that the count moved by the number of
tests you added. If something fails, ctest dumps the executable's verbose output and a
`Tests run: N, Failures: N, Errors: N` line — but only then.

## What the test build actually does

`run_tests.py` builds `CMAKE_BUILD_TYPE=CESM_DEBUG` by default; `--build-optimized`
selects `CESM` instead. The Fortran flags on Derecho with Intel are, in order:

```
-O0 -no-fma -g -check uninit -check bounds -check pointers -fpe0
-check noarg_temp_created -qno-opt-dynamic-align -convert big_endian
-assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source
-qopt-report -march=core-avx2 -check nouninit
```

Four consequences worth knowing before you write a fixture:

- **Uninitialized-variable checking is off**, despite `-check uninit` appearing in the
  list. `-check nouninit` comes *later* and last wins. It is appended by
  `ccs_config/machines/<machine>/intel_<machine>.cmake`, whose comment explains why. It
  cannot simply be turned back on: `ifx` implements `-check uninit` as MemorySanitizer
  instrumentation, and with `nouninit` removed the binary aborts inside glibc startup
  before any of your code runs — a backtrace through glibc, not a finding about your test.
  Re-check the flag order in a built `flags.make` when the machine files change.
- **`NDEBUG` is not defined**, so `#ifndef NDEBUG` blocks are live and count as coverage.
- **`-fpe0` turns a bad value into an abort, not a NaN.** Arithmetic on a signalling NaN
  stops the binary with `forrtl: error (65): floating invalid`, and a zero divisor with
  `error (73): floating divide by zero`. What traps is the *arithmetic*, not the read: a
  plain copy of a signalling NaN is silent, so a bad value can travel some way from the
  fixture that made it before the abort names a line. Either way, a fixture that leaves a
  divisor at zero does not produce a NaN result you can assert on.
- **`-check bounds` turns a broken index into an abort** rather than a failed assertion.

There are no `-warn` options at all.

## pFUnit authoring limits

The preprocessor rewrites `@` lines it can parse whole, and silently passes through the
ones it cannot. It reports success either way — `... Done.  Results in <name>.F90` — and
the compiler is what fails, so the error arrives one step removed from its cause.

**An `@assert` cannot be continued with `&`.**

```fortran
! Rejected
@assertEqual(expected_value(c), &
     actual_value(c), tolerance=tol)

! Fine - pull the operand into a local so the assert fits on one line
expected = expected_value(c)
@assertEqual(expected, actual_value(c), tolerance=tol)
```

**An `@assert` cannot carry a trailing comment.** Put the comment on the line above:

```fortran
! Rejected
@assertEqual(0._r8, h2osoi_liq(c))  ! the layer should be dry

! Fine
! The layer should be dry
@assertEqual(0._r8, h2osoi_liq(c))
```

Either mistake produces `error #5078: Unrecognized token '@' skipped` with a caret under
the `@`, then `error #5082: Syntax error, found '='`. The `.pf` is named through a doubled
path — the test's build directory, then the absolute path of the file. Everything after
the first two errors per file is cascade. Fix the first `#5078` and ignore the rest.

Use `message="..."` inside the assert when you want the explanation to appear in the
failure output rather than only in the source.

### Not every assert exists at every rank

In 4.8, `assertIsNaN` and `assertIsFinite` are **scalar only** — check an array element by
element, or assert on a reduction of it. The comparison asserts (`assertEqual`,
`assertLessThan`, `assertGreaterThan` and friends) exist for equal ranks, `(Nd, Nd)` up to
5d, and for a scalar against any rank, `(0d, Nd)` — but there is no `(1d, 0d)` form of any
assert, so a scalar has to be the expected value, on the left:

```fortran
! No specific procedure - a rank-1 expected against a scalar actual
@assertLessThan(bounds_array, upper_limit)

! Fine
@assertGreaterThan(upper_limit, bounds_array)
```

### Keep names under the symbol limit

`ifx` shortens any global name over 90 characters and warns
`#5462: Global name too long`. The symbol is `<module>_mp_<PROCEDURE>`, so your budget is
`len(module name) + 4 + len(test subroutine name) <= 90`, and truncation is from the
front. No CTSM test currently exceeds the limit, but the longest are in the high 80s, so
the margin is a few characters, not tens. The warning never appears in a green build log.

## Fixtures CTSM already provides

In `src/unit_test_shr/`:

| Module | What it gives you |
|---|---|
| `unittestSubgridMod` | `unittest_subgrid_setup_start` / `_setup_end` / `_teardown`, and `unittest_add_gridcell` / `_landunit` / `_column` / `_patch` for building a subgrid by hand |
| `unittestSimpleSubgridSetupsMod` | Ready-made subgrids: `setup_single_veg_patch`, `setup_n_veg_patches`, `setup_ncells_single_veg_patch`, `setup_landunit_ncols` |
| `unittestFilterBuilderMod` | `filter_from_range`, and `filter_empty` for the no-points case |
| `unittestWaterTypeFactory` | A complete `water_type`, without hand-rolling one. Its calls have a required order and one hidden side effect — see **`water-type-factory.md`** next to this file |
| `unittestArrayMod` | `grc_array` / `col_array` / `patch_array`, which build a subgrid-level 1-D array for you: called with no argument each returns an `r8` array prefilled with NaN, called with a value an array of that value's type filled with it. Also `logical_array_to_int` |
| `unittestTimeManagerMod` | `unittest_timemgr_setup`, `_set_curr_date`, `_set_nstep`, `_teardown` |

Worked examples, by file: `test_filter_col.pf` for filters and a simple subgrid,
`test_partition_precip.pf` for `setup_ncells_single_veg_patch`.

### Types with no factory

Hand-allocate only the components the routine under test actually touches, and:

```fortran
! in setUp - the pointer's association status is undefined until you do this,
! which makes an associated() guard in tearDown meaningless
nullify(this%temperature_inst%t_soisno_col)

! where you allocate it
allocate(this%temperature_inst%t_soisno_col(bounds%begc:bounds%endc, -nlevsno+1:nlevgrnd))
this%temperature_inst%t_soisno_col(:,:) = nan

! in tearDown
if (associated(this%temperature_inst%t_soisno_col)) then
   deallocate(this%temperature_inst%t_soisno_col)
end if
```

with `use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)`. For a plain 1-D
subgrid-level array, `col_array()` and its siblings above already do all of this.

**The `nan` prefill is the part people drop, and it is the part that matters.**
Uninitialized checking is off and cannot be turned on, so the only thing protecting you
is that CTSM's own init paths fill arrays with a *signalling* NaN, which `-fpe0` traps as
soon as anything computes with it. A bare `allocate` opts out of that, and then fails
**open**: fresh heap pages read as exactly `0.0`, so `@assertEqual(0._r8, ...)` against a
component nobody ever wrote will pass and tell you nothing.

For the allocate/deallocate pair alone, `test_init_columns.pf` is the plainest example.
`test_dyn_cons_biogeophys.pf` guards every deallocation with `associated()` — which is
exactly the case the `nullify` above exists to make meaningful.

## Setup traps

**Module state outlives a test.** A module variable set by one test is still set for the
next test in the same executable. Set what you depend on in `setUp` rather than inheriting
whatever the previous test left behind.

**Look for a testing entry point before writing to module state directly.** The naming is
not consistent — `InitForTesting`, `SetNMLForTesting`, `setParamsForTesting` and
`SnowHydrologySetControlForTesting` all exist — so grep the module for `ForTesting` rather
than guessing at a name.

**Setting a module's controls is not the same as initialising it.** A testing entry point
may set control variables without allocating the module's arrays, which the module's real
init routine does. If you skip that init, the arrays stay unallocated and the module's
matching cleanup routine fails trying to deallocate them. Check what actually allocates
the state you depend on.

## Running one test over several inputs

When the same behaviour has to hold across several values of a module-level configuration
variable, use pFUnit's parameterized test case rather than a loop inside one test body or
three copy-pasted subroutines. CTSM has no example of it, so the complete pattern, the two
pitfalls, and what the output looks like are in **`parameterized-tests.md`** next to this
file.
