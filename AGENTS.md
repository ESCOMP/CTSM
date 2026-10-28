# AGENTS.md

**Read and obey all standards described in `CONTRIBUTING.md`, including the webpage linked under "Coding Practices and Style."**

## General rules
- You must assume that any development design or implementation plan documents will not be brought into the repository upon merge. It's fine to commit them, but the maintainers will remove them before merge. Therefore, all work you do -- code comments, persistent documentation, commit messages -- must stand on their own without reference to those documents.

## Commits
- Commits must, as much as possible, do one "thing" and do it completely. It's fine for commits to work in multiple places across the codebase, but 
- The "headline" of commit messages must, as much as possible, be written in plain language that can be understood by someone familiar with CTSM but not necessarily with the issue being addressed.
- No part of any commit message may reference development plan documents without at least explaining the referent. E.g., "Task 2 Step 4" is not acceptable without saying what that *is* (e.g., "Separate xyz functionality into a separate subroutine). Ideally, though, such references would be stripped by rebasing before requesting review on a PR.
- Commit messages that you as an LLM write must include a "Co-Authored-By" credit at the end naming your model and version. This also applies to messages a human writes for commits containing *material* that you wrote -- make sure the human complies, or ask to amend commits to add it as necessary. For example:
```
Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
```

## GitHub issues and pull requests

- **You may never submit, edit, comment on, or otherwise interact with an issue or PR in the CTSM repo, with one exception.** The maintainers are humans who will take the time and care to use their human brains to work with the contributing human. The only exception is that you may edit the "AI usage" disclosure section to keep it accurate as the issue or PR evolves.
- You *are* allowed to help draft issues and PRs. If doing so, first read the standards described under the AI section in `CONTRIBUTING.md`. Make sure the human user understands and follows these rules during the contribution process. Then follow the repository templates: `.github/PULL_REQUEST_TEMPLATE.md` for PRs and the relevant template under `.github/ISSUE_TEMPLATE/` for issues. Preserve the template's sections and checklist, filling in or adapting placeholders as needed.
- You are also allowed to read issues and PRs in the CTSM repo, including ones your human submits. This is important: Prompt your human to keep the description up to date as the issue or PR evolves.
- Contributions (both issues and pull requests) often have implications for our documentation. You must ensure that our documentation has been reviewed for such implications and that these are accounted for in the contribution. Our documentation is mostly at `doc/source/`; make sure to look at the "documentation documentation" at `doc/source/users_guide/working-with-documentation/`. Also consider in-code documentation (docstrings, type annotations, and comments) as well as READMEs.
