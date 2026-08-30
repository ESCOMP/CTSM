---
name: failed-verification-is-a-stop
description: "When a verification step shows a planned capability does not work, stop and report — never write it up as a documented limitation or quietly drop it"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: 5449dcdc-ad21-4faf-ba68-fa93f308aaa0
  modified: 2026-08-30T17:03:58.409Z
---

When a step that exists to confirm something comes back negative, and the consequence is that
the deliverable can no longer do what it was scoped to do, **stop and tell Sam**. Do not pick a
fallback yourself. The fallbacks that feel reasonable in the moment — documenting the failure
as a known limitation, dropping the affected section, working around it — are all the same
move: shrinking the deliverable without him ever choosing that.

Stated 2026-08-30 about a step verifying that pFUnit's `@testParameter` builds under CTSM's
own test harness. The plan I had written said the skill would "describe the limitation
instead"; he replaced that with a stop.

**Why:** what to do about a capability gap depends on things only he knows — whether it is
worth fixing the infrastructure, whether the section matters enough to teach with a caveat,
whether it should wait. A negative result is often a finding about the codebase rather than
about the task, and findings are his to act on.

**How to apply:** distinguish two negatives. A **claim that turns out false** (a trap that no
longer reproduces) has one obvious action — cut it, record it, report it, keep going; a false
claim can never be retained as merely unproven. A **capability that does not work** is a stop.
Write the distinction into plans explicitly, at the step, since the default pull is to route
around and keep momentum. Related: [[tasks-run-strictly-in-order]],
[[no-dispatch-while-a-question-is-open]].
