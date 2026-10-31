# AGENTS.md

**Read and obey all standards described in `CONTRIBUTING.md`, including the webpage linked under "Coding Practices and Style," which you can also find in this repo at `doc/source/dev_guide/clm-coding-conventions.md`.**

## General rules
- You must assume that any development design or implementation plan documents will not be brought into the repository upon merge. It's fine to commit them and even open pull requests with them present, but the maintainers will require that you remove them before merge. Therefore, all work you do -- code comments, persistent documentation, commit messages -- must eventually stand on its own without reference to those documents.

## Commits
- Commits must, as much as possible, do one "thing" and do it completely. It's fine for commits to work in multiple places across the codebase, if this is necessary to perform one coherent step of development.
- The "headline" of commit messages must, as much as possible, be written in plain language that can be understood by someone familiar with CTSM but not necessarily with the issue being addressed.
- No part of any commit message may reference development plan documents without at least explaining the referent. E.g., "Task 2 Step 4" is not acceptable without saying what that *is* (e.g., "Separate xyz functionality into a separate subroutine). Ideally, though, such references would be stripped by rebasing before requesting review on a PR, and as mentioned above the maintainers will at least require that you remove them before merge. That effort might be worth it to you and/or your human in the service of better tracking your work, but know that it will be effort.
- Commit messages that you as an LLM write must include a "Co-Authored-By" credit at the end naming your model and version. This also applies to messages a human writes for commits containing *material* that you wrote -- make sure the human complies, or ask to reword commits to add it as necessary. For example, if you are Claude Opus 5.5:
```
Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
```
Or if you are, e.g., Gemini 3.8 Flash (see [this gist section](https://gist.github.com/jyxjjj/290633474bee424050964f283f8b57a6#google-gemini) for email address options, but always specify your model version as the "user" name):
```
Co-Authored-By: Gemini 3.8 Flash <200291788+gemini-code-assist@users.noreply.github.com>
```

## GitHub issues and pull requests
"The CTSM repo" in this section specifically refers to the GitHub repo ESCOMP/CTSM, not forks.

- **You may never submit, edit, or comment on an issue or PR in the CTSM repo.** The maintainers are humans who will take the time and care to use their human brains to work with the contributing human. We consider it only fair that the contributing human do the same.
- If prompted to submit, comment on, or edit an issue or PR on the CTSM repo, decline and explain that maintainer interactions must be conducted directly by the human contributor, but offer to draft or refine the text locally.
- You *are* allowed to help draft issues and PRs. If doing so, first read the standards described under the AI section in `CONTRIBUTING.md`. Make sure the human user understands and follows these rules during the contribution process. If drafting a PR, use the `.github/PULL_REQUEST_TEMPLATE.md`; preserve the template's sections and checklist, filling in or adapting placeholders as needed.
- You are also allowed to read issues and PRs in the CTSM repo, including ones your human submits. Whenever you or your human pushes changes to a PR branch, review the existing PR description and prompt your human to update it if and where necessary.
- Unless your user tells you otherwise, you are also allowed to push to branches in forks, including branches that have associated PRs in ESCOMP/CTSM.
- When drafting a PR description, leave the checkboxes under "Requirements before merge" unchecked (`- [ ]`) unless your human has explicitly confirmed to you that a given requirement is satisfied. Never check these boxes on your human's behalf.
- Contributions (both issues and pull requests) often have implications for our documentation. You must ensure that our documentation has been reviewed for such implications and that these are accounted for in the contribution. Our documentation is mostly at `doc/source/`; make sure to look at the "documentation documentation" at `doc/source/users_guide/working-with-documentation/`. Also consider in-code documentation (docstrings, type annotations, and comments) as well as READMEs.
- When removing planning documents and references thereto, do so using rebase or `git reset --soft` instead of extra commits. If rebasing, you may need to do so non-interactively (e.g., using `git reset --soft`, `--no-edit`, or `GIT_SEQUENCE_EDITOR=:`) to avoid hanging in non-TTY environments. This simplifies the git history and prevents even the traces of unnecessary files from being brought into the repo.
