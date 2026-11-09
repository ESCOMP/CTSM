# Assumed background

`SKILL.md` builds on two practices it does not teach. Where the superpowers plugin is
installed, `superpowers:test-driven-development` and `superpowers:subagent-driven-development`
are the full treatments and are worth reading instead of this file. This file exists so the
skill stands on its own where the plugin is not installed: it is the part `SKILL.md` actually
relies on, and no more.

## Test-first discipline

The order is the point, not just the coverage.

- Write the test before the code it covers. Run it. Watch it fail, and check that it failed for
  the reason you intended rather than because of a typo or a missing import.
- Then write the smallest code that makes it pass, and run it again.
- A test that has never been seen to fail is not known to be able to fail. Passing the first
  time it runs proves nothing about the test.
- Tests written afterwards answer "what does this do?". Tests written first answer "what should
  this do?" — which is the question the requirement asked.
- The iron law: no production code without a failing test first. "I'll test after", "I already
  checked it by hand", and "deleting what I already wrote would waste those hours" turn up every
  time, and none of them restores the missing evidence. Sunk cost is spent either way; the
  choice is only between code you can trust and code you cannot.

## A fresh agent per task, with review

- Work is decomposed into tasks. Each task gets a **fresh** agent that begins with no memory of
  the session: it knows only what its dispatch hands it, plus whatever it can discover for
  itself as a skill. Anything it must know has to be in one of those two places.
- A dispatch carries the task's requirement, the interfaces it touches, and the constraints that
  bind it — not the session's history and not the other tasks.
- Each task ends in a commit, then a review by another fresh agent that sees the task's
  requirement and the task's diff. Two stages are usual: does the diff meet the requirement and
  only the requirement, then is it correct and well made.
- Findings go back to the implementing agent, which fixes them and is re-reviewed.

That is the whole shape `SKILL.md` assumes. It adds one agent to it — a test author ahead of the
implementer — and the evidence rules that make the addition mean something.
