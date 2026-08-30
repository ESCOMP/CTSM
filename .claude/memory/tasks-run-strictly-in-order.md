---
name: tasks-run-strictly-in-order
description: "Every plan's tasks are strictly sequential — never start Task N+1 before Task N is finished, even when they look independent"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 5449dcdc-ad21-4faf-ba68-fa93f308aaa0
  modified: 2026-08-30T16:45:27.591Z
---

Treat every task in a plan as order-critical. Task N+1 does not begin until Task N is
finished. Stated 2026-08-30, in answer to a narrower question about whether two particular
tasks had to run in sequence — the answer generalized to all of them, so do not read task
independence into a plan just because two tasks touch different files.

**Why:** the reason he gave for the general rule is the same reason the specific case had:
finishing a task can change what a later task is even measuring. In the skills work, running a
baseline after the fix it is meant to baseline destroys the measurement with no way to recover
it. Overlapping tasks also make a review gate meaningless, since the reviewer can no longer
tell which task produced what.

**How to apply:** write plans so each task ends with an independently reviewable deliverable,
and do not batch or interleave, even to save a dispatch. Where an ordering is load-bearing
rather than merely conventional, say so *in the plan* and say what breaks if it is violated —
a plan that only implies its ordering will get reordered by someone optimizing for
parallelism. Related: [[scope-plan-sections-to-own-task]], [[unit-tests-are-test-first]].
