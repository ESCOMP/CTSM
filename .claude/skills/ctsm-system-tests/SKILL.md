---
name: ctsm-system-tests
description: Use when adding or changing a CTSM system test - an entry in cime_config/testdefs/testlist_clm.xml, a testmod under testmods_dirs, or an ExpectedTestFails.xml entry - when running or re-running existing system tests, whether or not you changed them, when choosing between a system test and a unit test, when setting a test's wallclock or its machine and compiler coverage, and when verifying a CTSM source change or deciding how far that verification goes before handing off to the human who runs the suites.
---

# CTSM system tests

A CTSM system test is a whole-model CIME case, registered as an entry in
`cime_config/testdefs/testlist_clm.xml`, configured by testmods under
`cime_config/testdefs/testmods_dirs/`, and run as part of one or more named suites.

## Start here

For the basics, read these rather than re-deriving them — everything below is what they
do not say.

- `doc/source/users_guide/testing/testing.rst` — CTSM's own testing overview, and the
  links it carries to the CIME testing chapter and the CTSM wiki's System Testing Guide.
- The header comment of `cime_config/testdefs/testlist_clm.xml` — what every suite is for.
- The header comment of `cime_config/testdefs/ExpectedTestFails.xml` — that file's schema.
- `./cime/scripts/create_test --help` — test type names, phases, and `--xml-*` selection.

## Running tests is the human's call

**`run_sys_tests` is the human's to launch, foreground or background — and so is any
`create_test` that names a whole category. This holds just as much for running existing
tests unmodified as for tests you added.** A suite may submit hundreds of jobs — the large
ones do — to a shared queue, billed to the allocation of whoever launched it, running for
hours on a schedule that person owns. Nothing about a code change makes that call yours.

Assistant-side verification stops at the case build and the unit tests (below). Past that,
say what you would run and why, and let the human run it. **Never report or characterize a
suite result you were not given** — not from a partial log, not from a job that is still
queued, and not by inference from something else that passed.

Where this repository's committed hook is installed
(`.claude/hooks/ctsm-suite-launch-guard.py`, wired in `.claude/settings.json`), a suite
launch turns into a permission prompt rather than running on an allow rule. That prompt
exists so a human sees the command before it runs; it is not a lane for getting a suite
launched on your own judgment.

## Symptom → cause → fix

| Symptom | Cause | Fix |
|---|---|---|
| `ERS_Ly2` (or any `ER*` with `STOP_N < 3`) fails immediately | `REST_N` is computed as roughly half of `STOP_N` plus one and then asserted strictly less than `STOP_N` | Respell the duration so `STOP_N >= 3` — `ERS_Ld731` rather than `ERS_Ly2` |
| A testmod's `hist_fincl1` is missing variables another testmod set | Assignment instead of `+=`, and testmods apply in order with later ones winning | `hist_fincl1 += 'VAR'`, unless the replacement was the point |
| A namelist path that looks right reaches `lnd_in` verbatim and the run cannot open it | Relative-to-absolute conversion applies to namelist *defaults*, never to a user override | `$SRCROOT/...` in the testmod |
| An input file was never fetched, and case setup never complained it was missing | A `landroot`-typed override does not enter `ctsm.input_data_list`, so CIME neither fetches nor existence-checks it | Not yours to paper over: tell the human, and offer either to pre-stage the file everywhere the test runs or to re-type the override so CIME does check it |
| An `ExpectedTestFails` entry did not apply to a test that failed | Wrong phase named, or a testid suffix was appended to the test name | Name the phase the test actually fails in; never add a testid |

## The verdict has to come from the phase statuses

**A completed test must never require a human to do anything beyond reading PASS/FAIL for
each of its phases.** You may not ask someone to diff two runs by hand, open a history
file, or confirm that two tests agreed with each other. If a comparison is the
requirement, it belongs *inside* a test type that performs it — and if you find yourself
writing "then check that X matches Y" in a testing note, that is the signal you picked the
wrong test type.

The types differ in what they assess, so choose deliberately (`create_test --help` lists
them):

- **`SMS`** reports that the run completed — plus, where a baseline exists, that the fields
  it compares match. In CTSM, completion also means every armed balance check stayed closed
  while it ran, since those call `endrun`.
- **`ER*`, `REP`, `PEA` and the other multi-run types** run the model more than once and do
  the comparison for you: exact restart, reproducibility, a different PE layout.
- **CTSM's own types** — seven of them, in `cime_config/SystemTests/`, all built on
  `SystemTestsCompareTwo` — compare two *configurations*. That base class is what to reach
  for when the requirement is "these two must agree", and one of the existing seven may
  already be it.

**Prefer a unit test where the requirement is about internal detail.** Which layer received
what, whether a weight was applied identically on both sides, whether a guard fired: a
stock `SMS` cannot see any of it, and "the run exercises it" is not an assessment. A unit
test that runs in seconds is both cheaper and stricter — see the `ctsm-unit-tests` skill.

**A CIME test that aborts is a FAIL.** So a requirement of the form "this configuration
must abort cleanly" cannot be a suite entry at all. It stays a manual check, named as such
wherever the change records its testing, rather than left looking covered.

## Registering a test

An entry names the test, grid, compset and testmods, then lists machine/compiler/category
rows and options:

```xml
<test name="SMS_D_Ld5" grid="f10_f10_mg37" compset="I2000Clm60Bgc" testmods="clm/default--clm/myfeature">
  <machines>
    <machine name="derecho" compiler="intel" category="aux_clm"/>
  </machines>
  <options>
    <option name="wallclock">00:20:00</option>
    <option name="comment">What this test is, and why it earns its place</option>
  </options>
</test>
```

**The full test name is the key** for both baselines and `ExpectedTestFails.xml`:
`TESTNAME.GRID.COMPSET.machine_compiler.testmods`, with `/` becoming `-` inside each
testmod name and `--` joining them. Duration is part of it, so a `Ly2` twin of an `Ld5`
test is a *different* test with no baseline until one is generated.

**Suite categories need no registration** — a category exists because a
`<machine ... category="..."/>` line names it. Confirm a new one took:

```
./cime/scripts/query_testlists --list categories
```

**But every category belongs in the header comment** at the top of `testlist_clm.xml`,
with a line saying what it is for. Nothing enforces that, which is exactly why a new
category is easy to leave undocumented and impossible for the next person to interpret.

**`ERS` and `ERP` need `STOP_N >= 3`** whatever the `STOP_OPTION` — see the table above for
the mechanism, in `cime/CIME/SystemTests/system_tests_common.py`. The two fail differently
and it is worth knowing which you are buying: `ERP` fails at SETUP, while `ERS` fails in the
RUN phase, after a full build has already been paid for. `SMS` is unaffected.

**Do not duplicate an `SMS` test exactly with an `ER*` test.** Vary the duration or the
configuration, so the second entry buys something.

**`<option name="comment">` says what the test is and why it exists**: type and duration,
site and climate, configuration, then the reason this entry earns its place — especially
if it is the only cover for something.

**Machine and compiler coverage.** derecho intel, derecho gnu and izumi nag catch
different classes of bug: intel debug traps NaN where gnu does not; nag's runtime pointer
checking catches illegal associations neither of the others sees; and gfortran rejects
integer-as-logical that Intel accepts as a DEC extension. **Not every test needs all
three** — aim for a *set* of new entries that reaches each of them, and put a given entry
on the combination that can actually catch what it is there for.

**Wallclock** is derived from existing entries, never guessed, and scales by resolution
rather than by whether a test is "global". The method and the anchor values are in
**`wallclock.md`** next to this file.

## Landing a test before the capability exists

Add the test early, with an `ExpectedTestFails.xml` entry, and give the change that makes
it pass an explicit step *removing* the entry:

```xml
<test name="SMS_D_Ld5.f10_f10_mg37.I2000Clm60Bgc.derecho_intel.clm-myfeature">
  <phase name="SETUP">
    <status>FAIL</status>
    <issue>#1234</issue>
  </phase>
</test>
```

- **Name the phase the test actually fails in.** A namelist-build guard aborts at SETUP,
  and a RUN entry would not cover it. A phase block may be repeated for a test that fails
  in more than one.
- **Do not append a testid suffix to the name** — that narrows the entry to a single run.
- **Never put an `ExpectedTestFails` entry on a test that is also serving as a baseline
  instrument.** A test recorded as an expected RUN fail generates no baseline, which
  silently removes bit-for-bit coverage for every later change — the entry does its job and
  costs you the thing you were relying on.

## Testmods, and four ways to get one silently wrong

Testmods apply in order and later ones win. Pick a composition order convention — feature
first, then site is one — and apply it consistently, so that a shared testmod is never
half-overridden depending on which entry composed it.

- **Extend history output rather than assigning it**: `hist_fincl1 += 'VAR'`. A plain
  assignment wipes what earlier testmods and the defaults set, invisibly, on every test
  that does not itself compose the testmod that set it. Assignment is the right tool when
  replacing the list *is* the point — a testmod that exists to output one narrow set of
  fields, say. Then it is a deliberate choice, and the testmod says so.
- **Point at an in-repo file with `$SRCROOT/...`.** XML variables expand in `user_nl_clm`,
  but the relative-to-absolute conversion applies only to namelist *defaults*, never to a
  user override, so a relative path reaches `lnd_in` verbatim.
- **A `landroot`-typed override never enters `ctsm.input_data_list`**, so CIME will not
  fetch it and will not check that it exists at case-setup time. The file must already be
  present on every machine the test runs on.
- **A feature that adds a required restart field cannot start from an existing `finidat`.**
  If your restart read aborts when its new field is absent, set `finidat = ' '` in the
  feature's testmod so its tests cold-start, and check whether
  `finidat_interp_source` is refused in that configuration too — interpolation is not a way
  out.

## Baselines

Whenever any test is added or changed, a **complete** new baseline is generated — including
tests the change did not touch. A `--generate` name that identifies the branch and commit
(`<branch-basename>.<shorthash>`) keeps them tellable apart; `--compare` names the previous
baseline. `run_sys_tests` requires one of `--compare`/`--skip-compare` and one of
`--generate`/`--skip-generate`, so the choice is always explicit.

The human launching the suite decides this, which means it is something to **prompt them
about** when a change has touched the test list — not something to arrange yourself.

## Verifying a CTSM change without a suite

These checks are yours to run:

- **Build a case.** A dedicated case built against the checkout catches everything the
  compiler can catch across the whole model, which the unit test build does not cover.
  **If your change added a new `.F90` and you are reusing a build directory, refresh the
  cache first:** CIME derives `Srcfiles` from `Filepath`, and adding a file to an existing
  source directory does not touch `Filepath`, so the build fails with
  `error #7002: Error in opening the compiled module file`. `./case.build --clean lnd`, or
  `touch` the `Filepath` under the case's `bld/.../clm/obj/`. This is a hazard of a
  long-lived verification case only — a CIME test builds in a fresh directory, and test
  directories should not be reused.
- **Run the unit tests** (`cime/scripts/fortran_unit_testing/run_tests.py`; the
  `ctsm-unit-tests` skill covers the traps, including an incremental rebuild that reports
  success while running the previous set of tests).
- **A change that touches namelist code additionally runs**
  `bld/unit_testers/build-namelist_test.pl` — see `testing.rst` for how to read its log.

Then stop, and say plainly what was run and what was not.

## Amending commits

**Amend freely until a human says a system test has been run against a commit; never
after.** From that point a fix is a new commit, because the commit someone tested has to
keep existing for the result to mean anything. If it is ambiguous whether anything has been
run against it, ask rather than assuming it is still yours to rewrite.

Note that commits which add or change *unit* tests carry a stricter rule of their own —
they are never amended once run. See the `writing-tests-before-the-implementer` skill.
