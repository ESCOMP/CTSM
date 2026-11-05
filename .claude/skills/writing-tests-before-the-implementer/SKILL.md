---
name: writing-tests-before-the-implementer
description: Use when a task's unit tests are written by a different agent than the one implementing it, when deciding whether a test is committed before or after the code it covers, when a committed test fails and the implementation does not exist yet, when tempted to edit a test file to make it pass, when a test passed the first time it ran, when amending a commit that adds or changes a test, or when checking whether a task's tests were reshaped to match the implementation that landed.
---

# Writing tests before the implementer

A dedicated test-writing agent runs **before** the implementing agent. It writes the task's
tests, runs them, and commits them on its own. The implementer that follows may not edit a
test file: making the tests pass without touching them is the job.

**Core principle:** a test written after the implementation encodes the bug as readily as the
fix. The guard is not intent, it is evidence — and evidence has to be produced before the code
exists, because afterwards there is nothing left to produce it against.

**Violating the letter of the rules is violating the spirit of the rules.**

This matters most where it looks least necessary. On a one-slot index rewrite, what the code
does and what the code *should* do look identical to a reader. A test written afterwards agrees
with whichever one it was written against.

Exact commands and per-framework failure modes are in **`test-frameworks.md`** next to this
file; everything below is framework-independent.

## What this adds

This skill assumes two practices it does not teach: ordinary test-first discipline, and a
workflow that dispatches a fresh agent per task with a review after it. Where the superpowers
plugin is installed, `superpowers:test-driven-development` and
`superpowers:subagent-driven-development` are those two in full, and this skill is an addition
to that family. Where it is not, **`assumed-background.md`** next to this file states the part
this skill relies on. Nothing below needs the plugin.

Neither of those covers a *separate* test author. Test-driven development assumes one agent
alternating between test and code; a per-task dispatch gives each task one implementer that
"implements, tests, commits". Four things are added here:

| Added | Because |
|---|---|
| A test author who is not the implementer | An implementer with a failing test in front of it can make it pass by editing it |
| The red-evidence commit | A commit anyone can check out and watch fail — verifiable with `git`, not on trust |
| A mutation duty for tests that are never red | Red/green says nothing about a test that is green at both ends |
| The anti-retrofit diff check | The only mechanical proof that no test was reshaped to match the code |

## The arrangement

1. **The test writer is dispatched first.** It receives the task's requirement, the spec, and
   the constraints — **never a description of how the implementation will be shaped**, so what
   it writes encodes the requirement rather than someone's guess at the code.
2. **It runs the tests itself.** A test file that has not been executed must not be committed:
   an unrun test's failure is indistinguishable from a typo. Confirm the tests really ran —
   a suite that reports success without having picked them up is worse than a failure, because
   it looks like evidence.
3. **It labels every test** red-first or green-throughout (below) in its hand-off, and produces
   the mutation evidence for the green ones.
4. **It commits the tests standalone.**
5. **The implementer's dispatch states that it must not edit a test file.** Say it explicitly.
   Failing tests in front of an implementer are an invitation to adjust the tests, which
   silently converts this back into writing them afterwards.
6. **The reviewer sees the union of the test commit and the implementation commit.** A reviewer
   shown only the implementation commit cannot tell whether the tests were written to the
   requirement or to the code. The reviewer also confirms the test files were not touched after
   they were run.

## Two tracks — label every test as one

**Red-first — behaviour the change alters.** It must **fail at the test commit and pass at the
end of the task**. That pair of commits is its evidence: check out the test commit, run the
suite, watch it fail. This is stronger than a mutation, because the mutation is the pre-change
code itself rather than an edit chosen after the fact by whoever wrote the test.

**Green-throughout — behaviour the change must preserve.** Green before and green after, so
red/green proves nothing about it and it **owes a mutation instead** (below). Committing it
early still buys what a test written afterwards cannot: green, then an untouched test file,
then green again is an invariance proof, where a test written at the end only pins the end
state.

**If no evidence can be produced for a test — neither a red run nor a mutation — it is not
pinning anything. Cut it.**

## The test commit

**Standalone, and never amended once it has been run.** Amending destroys the red evidence: the
commit a reviewer would check out to watch the test fail stops existing, and an amendable test
commit makes the diff check below meaningless. This **overrides any amend-the-commit review
loop** for these commits specifically — a review finding against a test becomes a new commit.

**It has to actually run.** A test that calls something which does not exist yet is not a red
test: it does not fail, it fails to build or to be collected, and a test that never ran
produced no evidence. The blast radius is usually wider than the one test — in a compiled suite
a single bad test file can abort the whole run before anything executes, taking every other
test's result with it. So where a change adds a new callable, its tests commit at the earliest
point whose code makes them compile, and the task text says which point and why. The common
case is unaffected: a test that sets up the new configuration, calls an existing routine
through its existing interface, and fails on the assertion.

**It is exempt from "all tests pass", and only from that.** Its red cases fail there by design.
It must still build, and every test that existed before it must still pass. Every later commit
in the task is held to the full gate.

## Mutation evidence

For every green-throughout assertion: break the thing under test, confirm **that specific
assertion** fails, restore it, and report which assertion caught which break. A test written
against working code proves nothing until you have seen it fail.

- **Reverting the line under test to its pre-change form beats inventing a mutation.** An
  invented mutation is a guess at what could go wrong; the reversion is the thing that actually
  would have.
- **Two-sided.** Red on the new case, green on its unchanged mirror. A mutation that reddens
  both has not isolated anything.
- **A mirror test is blind by design** to a mutation of the new-path-only term — that is what
  makes it a mirror. Mutate it the other way instead: make the new-path term unconditional.
- **Narrow the mutation until it reaches your assertion.** Assertion frameworks stop a test at
  its first failure, so a mutation that trips an *earlier* assertion means the one you meant to
  exercise never ran. A test that stopped at assertion 2 tells you nothing about assertion 7.
- **Claiming a test adds coverage? Run the control.** The mutation must fail with the new test
  present and pass with it absent.
- **A crash counts.** A mutation may abort or error the test rather than fail an assertion — a
  bounds check, a trapped arithmetic exception, an exception where an assertion failure was
  expected. The test still caught it; record which form you saw.

## Retrofitting coverage to code that already landed

There is no red-first track at all: the whole task is green-throughout and the writer expects
green. The mutation duty then **splits**. The writer mutates every assertion it lands, **and
the reviewer independently mutates too**, choosing its own mutations rather than re-running the
writer's. A test that detects only the mutation its author had in mind reads as coverage
without being coverage.

## The anti-retrofit check

At the end of the task, from the **repository root**:

```bash
git diff <test-commit>..HEAD -- ':(top)<test-file glob>'
```

It must come back empty. Two ways this check quietly lies:

- **Pathspecs are relative to the current directory.** Run from inside a test directory, the
  same command reports on that directory alone and comes back empty for every test file outside
  it. Anchor it with `:(top)` — or `:/` — or run it from the repository root. Quoting protects
  you from the shell, not from this.
- **A test-file glob does not cover all test-side code.** Shared fixture harnesses, test
  helpers, per-directory test configuration and a test directory's build list do not match it, and a red
  test can be turned green by editing any of them. Check the implementer's whole diff for
  test-side files, not just the glob.

**A non-empty diff is a finding either way:** either a test was retrofitted to the
implementation, or the requirement was misread before the tests were written.

## A test that has to change is a finding, not a formality

It means the requirement was misread before implementation started. Allowed — never silent:
**its own commit**, the reason in the hand-off, and the requirement text corrected to match.
Never folded into the implementation commit.

## Rationalizations

| Excuse | Reality |
|---|---|
| "The test has an obvious typo, I'll just fix it" | Then no one can tell your fix from a retrofit. Report it; it gets its own commit. |
| "The test asserts the wrong thing, so editing it *is* the work" | Maybe. That is a finding about the requirement, not an edit to make in passing. |
| "I'll write the tests once I can see the shape of the code" | That is tests-after with extra steps. The shape of the code is exactly the bias the order exists to exclude. |
| "The test can't be red — the thing it calls doesn't exist yet" | Correct, and that means it is not a red test but a build or collection error, which runs nothing. Commit it where it runs. |
| "It passes, so the behaviour is covered" | A test that has never failed is not evidence that it can fail. |
| "My mutation didn't fail, so the code must be right" | Or the mutation never reached your assertion. Check where the test stopped. |
| "I'll amend the test commit — it's tidier, and nothing is pushed" | Amending deletes the one commit the evidence lives in. New commit. |
| "I'll fold the test correction into the implementation commit" | Then the test change is invisible in the log, which is the whole failure mode. |
| "The diff check came back empty" | From which directory, with which pathspec, over which files? An empty result from a subdirectory means nothing. |
| "I made the test pass by fixing the fixture module, not the test" | The check is blind to that. It is the same violation with a different file extension. |
| "It's red at the test commit, so no mutation is needed" | True for that test. Every green-throughout test in the same commit still owes one. |
| "The suite came back green, so my new test passes" | Or it never ran. Confirm the count moved by the number of tests you added. |
| "Both reviewers would just re-run the same mutation" | On retrofitted coverage that is the point of failure — the reviewer picks its own. |

## Red flags — STOP

- Editing a test file while implementing
- A test that passed the first time it ran, reported with no mutation
- `git commit --amend` anywhere near a commit that added or changed a test
- A test committed without having been run, or without confirming it ran
- A red test made green by editing a fixture harness, a build list, or a test helper
- Reviewing the implementation commit without the test commit
- An empty diff check run from anywhere but the repository root
- "I'll add the tests after this builds"
