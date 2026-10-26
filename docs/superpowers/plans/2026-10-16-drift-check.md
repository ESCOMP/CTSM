# Drift check (NEXT_STEPS item 5) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One checker script behind two workflows: a PR gate that blocks merges when the recipe and `ccs_config` disagree, and a scheduled monitor that compares the published image and the recipe against derecho's live stack and opens an issue on drift.

**Architecture:** The container side is split by what can change it. The executable test (does a case build link the image's mpi-serial and PIO, not CTSM's submodules) runs where a build already happens — as a post-build assertion in `cirrus-testing.yml`'s `simple-build-create_test` job — because only a repo change can break it. The scheduled monitor, which watches the side that moves without commits, runs the pinned image to read its as-built facts, reads derecho's module tree and `spack.lock` from `/glade`, and compares both against the Dockerfile ARGs. A hosted-runner heartbeat job notices when the self-hosted runner has gone quiet.

**Tech Stack:** Python 3 stdlib (the existing checker's style, `unittest`), bash, GitHub Actions (`gha-runner-ctsm` self-hosted + `ubuntu-latest`), `JasonEtco/create-an-issue` (already pinned in `docs-update-ctsm_pylib.yml`).

**Spec:** `docker/ctsm-ci-derecho-gnu/DRIFT_CHECK_DESIGN.md` plus `docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md` item 5. Where they disagree with each other, the decisions in this plan's Context section win, and Task 1 writes those decisions back into both.

**Status: blocked on Jira CCPP-502** (the read-only `/glade/u/apps` mount on `gha-runner-ctsm`). The user's decision on 2026-10-16: no implementation on item 5 until that ticket is answered. Tasks 2–4 and 7 do not technically depend on the ticket; Tasks 5 and 6 do (they need `/glade/u/apps` to be validated), and Appendix A exists only if the ticket is refused. This file is `docs/superpowers/specs/2026-10-16-drift-check-plan.md`; NEXT_STEPS item 5 points here.

## Context

Item 5 asks for "one scheduled thing that compares built executables". The documents disagree on what that means for the container side: the design's premise (`DRIFT_CHECK_DESIGN.md:15-16`) is executables, its Phase 2 (`:207`) and NEXT_STEPS 4b (`:429`) say item 5 reads `stack-facts.txt` with `podman run … cat`, and the "one script, two workflows" table (`:127-130`) plus NEXT_STEPS item 5 make the derecho side metadata (module tree + `spack.lock`), never a derecho-built binary. The mpi-serial gap that motivated "executables" was closed by mechanism in item 4b, which made the ARGs load-bearing and left the binary with one job: confirming the mechanism via the `libraries/mpi-serial` strings test.

Decisions (user, 2026-10-16):

- **B+.** The strings assertion runs after every CI case build (`cirrus-testing.yml`). The scheduled monitor reads as-built facts from inside the pinned image (no case build) and compares them with the Dockerfile ARGs and with derecho's live tree. Derecho's binary is never built in the loop; the design's "Tier 3" derecho binary is a one-off manual validation that the metadata read is a faithful proxy.
- **Heartbeat.** A hosted-runner job in the scheduled workflow opens an issue when no scheduled run has *completed* (any conclusion) in 10 days — one missed week is caught at the next.
- **Derecho side: wait for CCPP-502.** The live reader is built either way (Task 5). If the mount is refused, Appendix A adds the cron-written snapshot and the stale-report failure.
- **All conflicting documentation gets reconciled** (Task 1 first; each later task updates the README where it changes what the README says).

Verified on 2026-10-16 from a glade login node, so the implementer does not have to re-derive it: `/glade/u/apps/derecho/modules/25.10/gcc/14.3.0/netcdf/4.9.3.lua` line 37 is `setenv("NETCDF", "/glade/u/apps/derecho/25.10/spack/opt/spack/netcdf/4.9.3/gcc/14.3.0/s6hn")` (also `NCAR_ROOT_NETCDF`, line 39); `/glade/u/apps/derecho/25.10/envs/public/spack.lock` has top-level keys `_meta, spack, roots, concrete_specs`; `concrete_specs` is keyed by 32-char hash, exactly one key starts with `s6hn`, that spec is `name=netcdf version=4.9.3`, its `dependencies` is a list of `{"name", "hash", "parameters"}` and the transitive closure holds exactly one version each of `netcdf-c 4.9.3`, `netcdf-fortran 4.6.2`, `hdf5 1.14.6`.

## Global Constraints

- **No image rebuild.** Everything the monitor reads from the image must be readable from `ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261015` as published, from inside a `container:` job (no `podman inspect`, no labels). Dockerfile edits are comment-only.
- **No `ccs_config` changes.** Settled project decision.
- **Python: stdlib only** in `check-derecho-versions.py` and every new script (they run under `actions/setup-python` 3.x on hosted runners, on the self-hosted runner, and inside the image's `python3`). Tests are `unittest`, invoked as `python docker/ctsm-ci-derecho-gnu/test_<name>.py -v`. On the login node use the `ctsm_pylib` env's python (carried over from the previous plan; do not create or survey environments). Match the checker's existing black-formatted style.
- **Exit codes:** the checker keeps `0` = all checks pass, `1` = a `❌` line, uncaught exception = setup problem. New scripts use `0`/`1`/`2` as stated per task.
- **Report lines** keep the checker's conventions: one `✅`/`❌`/`⚠️` line per result, failures say what to do, summary line last.
- **Only run the workflow being worked on.** `derecho-drift.yml` may be exercised from the branch (via its own `pull_request` self-trigger or `workflow_dispatch` once on the default branch). Never dispatch other workflows to answer a question; ask the user.
- **The user runs Casper.** Anything needing a compute allocation or a built `cesm.exe` is a user-run verification step, written out in full.
- **Never write under `/glade/campaign`.** That write (Appendix A's cron) is the user's to make.
- **Pinned actions:** copy SHAs from existing workflows (`actions/checkout@11bd719…`, `actions/setup-python@a26af69…`, `JasonEtco/create-an-issue@1b14a70…`; upload/download-artifact from `docs.yml`).
- **Docs:** never hard-wrap Markdown prose; one line per paragraph or list item. Each task updates `README.md`/`NEXT_STEPS.md` wherever it changes what they describe.
- **Every task opens with Step 0, done by the orchestrator itself, never an implementer:** (a) review what earlier tasks landed, (b) re-read the task's own text for anything that work or new measurements made stale, and amend it, (c) collect outstanding questions. Then **stop for the user's approval, no matter what** — even when nothing changed and there are no questions. The message is brief: what was changed or added to the task (or "nothing"), and numbered questions at the bottom. No implementer is dispatched until the user approves.
- **Every task ends by committing, then stopping** for the user to review. Commit before the stop.
- **Commits:** body ≤ 200 words, wrapped at 72, honest verification status ("not run in CI" where true), the attribution trailer the session's instructions specify (no model name is hardcoded here: the plan may be implemented with a different model than wrote it).

## Review Focus

1. A PR from a fork: `github.event.pull_request.head.sha` is a commit the base repo only has via the merge ref; with `fetch-depth: 0` it is reachable (second parent of the merge commit). → Task 2 test `test_head_sha_reachable_only_via_merge_commit`.
2. A PR whose base is a release branch, not `master`. → Task 2 test `test_base_ref_is_not_master`.
3. Two `concrete_specs` keys sharing the 4-char prefix. → Task 5 test `test_ambiguous_hash_prefix_is_a_finding`.
4. The `[snapshot]` ini drifting from live while every ARG still matches live — the PR gate would keep passing on a stale snapshot. → Task 5 test `test_stale_snapshot_reported_against_live`.
5. A tool missing inside the image job (`nc-config` off `PATH`): a partial facts file must not read as "matches". → Task 4 test `test_missing_tool_aborts_without_writing`.

---

### Task 1: Decision record

**Files:**
- Modify: `docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md` (item 5, `:442-566`; the 4b sentence at `:429`; header date `:3`)
- Modify: `docker/ctsm-ci-derecho-gnu/DRIFT_CHECK_DESIGN.md` (`:15-16`, `:141-191` intro, `:239-273` Phase 5, `:305`)

**Interfaces:** none. Prose only.

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what earlier tasks landed; re-read this task for anything they or new measurements made stale, and amend it; collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Rewrite NEXT_STEPS item 5's lead paragraph** (`:442-447`) so it states the B+ shape in the Context section's words: executable assertion per CI build; scheduled monitor = image facts + recipe vs derecho live metadata; derecho binary never built in the loop. Keep the pointer to `DRIFT_CHECK_DESIGN.md` and to this plan.
- [ ] **Step 2: Add a dated "Decisions 2026-10-16" block** under item 5 with the four decisions from Context, the heartbeat threshold (10 days, `--status completed`), the issue mechanism (`JasonEtco/create-an-issue`, `update_existing: true`, `search_existing: open`, one stable title per failure kind), the `workflow_dispatch` input `open_issue`, and the blocked-on-CCPP-502 status with the two derecho-side paths. Replace the `[snapshot]`/`spack.lock` bullet's "verified 2026-10-09" paragraph with the exact `setenv("NCAR_ROOT_NETCDF", …)` → hash-prefix → closure method from Context (verified 2026-10-16). Change the "cron fallback" bullet's first sentence to say it is built only if CCPP-502 is refused (Appendix A of this plan).
- [ ] **Step 3: Fix the 4b sentence** at `:429` ("item 5 reads them with `podman run ... cat`") to "item 5 reads them from inside a `container:` job on the Cirrus runner (`image_facts.py`)".
- [ ] **Step 4: Reconcile DRIFT_CHECK_DESIGN.md.** Under "two things changed" (`:15-16`) add one dated sentence: the executable comparison lives in the per-build strings assertion; the schedule compares image facts and recipe against derecho's live metadata. In the instrument-map intro (`:147`) say the derecho column is the one-off manual validation of the metadata proxy, not something the monitor runs. Replace Phase 5's `derecho-drift.yml` bullet (`:250-252`) with the four jobs (`image-facts`, `check`, `report`, `heartbeat`) and note `schedule`/`workflow_dispatch` leave the PR gate. At `:305` replace "Tier 3" with "the manual derecho-binary validation (instrument map, derecho column)".
- [ ] **Step 5: Bump the header date** in NEXT_STEPS (`:3`).
- [ ] **Step 6: Commit** — `docs: record item 5 decisions in NEXT_STEPS and the design doc`. Then stop and ask the user to review.

---

### Task 2: PR gate hygiene — `derecho-version-check.yml` + staleness guard

**Files:**
- Create: `docker/ctsm-ci-derecho-gnu/pr_gate_staleness.py`
- Create: `docker/ctsm-ci-derecho-gnu/test_pr_gate_staleness.py`
- Modify: `.github/workflows/derecho-version-check.yml`
- Modify: `docker/ctsm-ci-derecho-gnu/README.md` ("Bumping versions", `:721-765`: the PR gate no longer has a scheduled half; describe the guard)

**Interfaces:**
- Consumes: nothing.
- Produces: `read_pull_request_paths(workflow_path: str) -> list[str]` — the `on: pull_request: paths:` items, in order; raises `ValueError` if the block is absent/empty or an entry contains `*`, `?` or `[` (the guard uses literal pathspecs; the current list is all literals). `changed_watched_paths(base_ref: str, head_sha: str, paths: list[str], cwd: str) -> list[str]` — `git merge-base <head_sha> origin/<base_ref>`, then `git diff --name-only <merge_base> origin/<base_ref> -- <paths…>`, returned as a sorted list. CLI: `pr_gate_staleness.py --workflow PATH --base-ref REF --head-sha SHA` → exit 0 and `✅ base <ref> has not moved in watched paths since <merge-base[:7]>`; exit 1 and `❌ stale: base <ref> moved in watched paths since this branch forked (<comma-separated paths>). This run evaluated the merge commit, so the result is a stale green rather than a wrong answer; merge or rebase <ref> into this branch and push to re-run.`; exit 2 on `ValueError`/git failure, message on stderr.

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what earlier tasks landed; re-read this task for anything they or new measurements made stale, and amend it; collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Write the failing tests** in `test_pr_gate_staleness.py`. A `setUp` builds a temp repo with `git init -b main`, one commit holding `Dockerfile`, `ccs_config` (a plain file standing in for the gitlink), `unrelated.txt`, and a fixture workflow whose `on.pull_request.paths` lists `Dockerfile` and `ccs_config`; creates branch `feature` with one commit; sets `refs/remotes/origin/main` with `git update-ref`. Helper `run_guard(base_ref, head_sha)` runs the CLI via `subprocess` and returns `(returncode, stdout, stderr)`. Tests:
  - `test_irrelevant_base_move_passes`: commit to `main` touching `unrelated.txt`, update the origin ref → rc 0, stdout starts with `✅`.
  - `test_watched_base_move_fails`: commit to `main` touching `ccs_config` → rc 1, stdout contains `❌ stale` and `ccs_config`.
  - `test_head_sha_reachable_only_via_merge_commit`: after a watched base move, `git merge main` into `feature` (HEAD is now the merge commit) and pass `--head-sha` = the pre-merge `feature` tip → rc 1; additionally assert `git merge-base --is-ancestor origin/main HEAD` succeeds, documenting the trap the guard avoids.
  - `test_base_ref_is_not_master`: repeat the watched-move case with the base branch named `release-clm5.0` → rc 1.
  - `test_real_workflow_paths`: `read_pull_request_paths` on the real `.github/workflows/derecho-version-check.yml` returns a list containing `docker/ctsm-ci-derecho-gnu/Dockerfile` and `ccs_config` and has no duplicates.
  - `test_no_paths_is_setup_error`: a fixture workflow without `paths:` → rc 2, `ValueError` text on stderr.
  - `test_glob_entry_rejected`: an entry `docker/**` → `ValueError`.
- [ ] **Step 2: Run** `python docker/ctsm-ci-derecho-gnu/test_pr_gate_staleness.py -v` → FAIL (module not found).
- [ ] **Step 3: Implement `pr_gate_staleness.py`.** `read_pull_request_paths` is a line-based parser: find the line matching `^\s*pull_request:\s*$`, then within its indentation block the line `^\s*paths:\s*$`, then collect `^\s*-\s*['"]?(.+?)['"]?\s*$` lines until the indentation drops back; no YAML library. `changed_watched_paths` uses `subprocess.run(..., check=True, capture_output=True, text=True)`.
- [ ] **Step 4: Run the tests** → all PASS.
- [ ] **Step 5: Edit `derecho-version-check.yml`.** `on:` becomes `pull_request:` only (keep `paths:`; add `docker/ctsm-ci-derecho-gnu/pr_gate_staleness.py` and `docker/ctsm-ci-derecho-gnu/test_pr_gate_staleness.py`); delete `push`, `schedule`, `workflow_dispatch` and the per-event `runs-on` expression (`runs-on: ubuntu-latest`); add `permissions: contents: read`; checkout gets `fetch-depth: 0`; after the existing test step add `python docker/ctsm-ci-derecho-gnu/test_pr_gate_staleness.py -v`; then a step `Staleness guard` running the CLI with `--base-ref "${{ github.event.pull_request.base.ref }}" --head-sha "${{ github.event.pull_request.head.sha }}"`; the checker step always passes `--skip-module-tree` (drop the conditional). Update the file's header comment: the scheduled half now lives in `derecho-drift.yml` (Task 6).
- [ ] **Step 6: Update README "Bumping versions"** to describe the gate as PR-only with the staleness guard, and point the scheduled description at Task 6's workflow (write it as "`derecho-drift.yml`, see below" — Task 6 fills it in).
- [ ] **Step 7: Commit** — `ci: make derecho-version-check a PR-only gate with a staleness guard`. State in the body that the workflow change is not CI-verified until a PR touching a watched path runs it. Then stop and ask the user to review.

---

### Task 3: Executable assertion after every CI case build

**Files:**
- Create: `docker/ctsm-ci-derecho-gnu/assert-links-image-stack.sh`
- Create: `docker/ctsm-ci-derecho-gnu/test_assert_links_image_stack.py`
- Modify: `.github/workflows/cirrus-testing.yml` (`simple-build-create_test`, after the `create_test` step)
- Modify: `docker/ctsm-ci-derecho-gnu/README.md` (mpi-serial section `:489-523`: the acceptance test now runs in CI; name the script)
- Modify: `docker/ctsm-ci-derecho-gnu/cime-macros/gnu_container.cmake` (the "ACCEPTANCE TEST" comment: name the script)

**Interfaces:**
- Consumes: nothing.
- Produces: `assert-links-image-stack.sh <cesm.exe>` — runs `strings -a` once into `$TMPDIR`; exit 1 with `❌ cesm.exe compiled CTSM's submodule: <match>` if any line matches `libraries/mpi-serial/` or `libraries/parallelio/` (path-prefix-independent: the repo is `/ctsm` on Casper but `/__w/CTSM/CTSM` in CI); exit 1 with `❌ no compiled-source path for /usr/local/mpi-serial found; a non-DEBUG or stripped binary makes this test vacuous` if no line contains `/usr/local/mpi-serial`; otherwise exit 0 and `✅ cesm.exe links the image's mpi-serial and PIO (<n> image-path hits, 0 submodule hits)`. Exit 2 if the file is missing or `strings` is not on `PATH`.

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what earlier tasks landed; re-read this task for anything they or new measurements made stale, and amend it; collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Write the failing tests.** `strings -a` on a text file yields its printable lines, so fixtures are text files. Tests: `test_clean_binary_passes` (lines include `/usr/local/mpi-serial/include/mpif.h`), `test_submodule_path_fails` (adds `/__w/CTSM/CTSM/libraries/mpi-serial/mpi.c` → rc 1, message names it), `test_parallelio_submodule_fails`, `test_no_positive_control_fails` (no `/usr/local/mpi-serial` line → rc 1, "vacuous"), `test_missing_file_is_setup_error` (rc 2).
- [ ] **Step 2: Run** `python docker/ctsm-ci-derecho-gnu/test_assert_links_image_stack.py -v` → FAIL.
- [ ] **Step 3: Implement the script** (`set -euo pipefail`, `mktemp` under `${TMPDIR:-/tmp}`).
- [ ] **Step 4: Run the tests** → PASS.
- [ ] **Step 5: Add the CI step** to `simple-build-create_test`, after `Call create_test`: find the executable with `find "$RUNNER_TEMP" -path '*/bld/cesm.exe' -type f`, fail if the count is not exactly 1, then run the script on it. Keep it in its own step named `Assert cesm.exe links the image's mpi-serial and PIO`.
- [ ] **Step 6: Update README and the macro comment** as listed.
- [ ] **Step 7: Commit** — `ci: assert each CI case build links the image's mpi-serial and PIO`. Body: the CI step is unverified because `simple-build-create_test` is red (NEXT_STEPS item 13); the script itself was tested on fixtures only. Then stop and ask the user to review; Step 8 is theirs to run.
- [ ] **Step 8 (user, Casper): validate against a real binary.** From a Casper allocation with the image loaded, after any `run-test-in-container.sh` run: `docker/ctsm-ci-derecho-gnu/assert-links-image-stack.sh "$(find $SCRATCH/cases_devcontainer -path '*/bld/cesm.exe' | head -1)"` → expect `✅`. Record the date in README's mpi-serial section.

---

### Task 4: Image facts from inside the image, and the checker's image-vs-recipe comparison

**Files:**
- Create: `docker/ctsm-ci-derecho-gnu/image_facts.py`
- Create: `docker/ctsm-ci-derecho-gnu/test_image_facts.py`
- Modify: `docker/ctsm-ci-derecho-gnu/check-derecho-versions.py` (`run_checks`, `main`, module docstring)
- Modify: `docker/ctsm-ci-derecho-gnu/test_check_derecho_versions.py`
- Modify: `docker/ctsm-ci-derecho-gnu/Dockerfile` (comment only, the `stack-facts` block near `:836`)
- Modify: `docker/ctsm-ci-derecho-gnu/README.md` ("Bumping versions": what the image facts are and where they come from)

**Interfaces:**
- Consumes: nothing.
- Produces (`image_facts.py`): `gcc_version(run) -> str` from `gfortran -dumpfullversion`; `netcdf_c_version(run) -> str` from `nc-config --version` (`netCDF 4.9.3` → `4.9.3`); `netcdf_fortran_version(run) -> str` from `nf-config --version` (`netCDF-Fortran 4.6.2` → `4.6.2`); `hdf5_version(h5pubconf_text: str) -> str` from `#define H5_VERSION "1.14.6"`; `esmf_facts(esmf_mk_text: str) -> dict` with keys `esmf` (from `ESMF_VERSION_STRING=`) and `esmf_comm` (from `ESMF_COMM=`); `stack_facts(text: str) -> dict` parsing `key=value` lines. `run` is `Callable[[list[str]], str]` returning stdout, injected for tests (default: `subprocess.run(..., check=True, capture_output=True, text=True).stdout`). `collect(run, h5pubconf_path="/usr/local/include/H5pubconf.h", esmf_mk_path=os.environ["ESMFMKFILE"], stack_facts_path="/opt/ctsm-container/stack-facts.txt") -> dict` with exactly these keys: `gcc, netcdf_c, netcdf_fortran, hdf5, esmf, esmf_comm, netcdf_fortran_stackfacts, netcdf_c_stdfilters`. `write_facts(facts: dict, out_path: str)` writes `key=value` lines sorted by key. CLI `image_facts.py --out PATH`: exit 1 and nothing written if any source is missing or a parser finds no match.
- Produces (checker): `load_image_facts(path) -> dict` (`FileNotFoundError` if absent — a setup problem); `check_image_facts(args: dict, facts: dict) -> (ok, lines)` with `IMAGE_FACT_ARGS = [("GCC_VERSION","gcc"),("NETCDF_C_VERSION","netcdf_c"),("NETCDF_FORTRAN_VERSION","netcdf_fortran"),("HDF5_VERSION","hdf5"),("ESMF_VERSION","esmf")]`; each line `✅ image <key> <v> matches ARG <ARG>` or `❌ image <key> is <v> but the Dockerfile ARG <ARG> is <a>: the published image was not built from this recipe`; plus `❌ image esmf_comm is <x>, expected mpiuni` and `❌ image netcdf_fortran_stackfacts <x> ≠ nf-config <y>`; a missing key is `❌ image facts: no "<key>" line`. `run_checks(..., image_facts=None)` runs it when a path is given; CLI flag `--image-facts PATH`.

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what earlier tasks landed; re-read this task for anything they or new measurements made stale, and amend it; collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Write the failing tests** in `test_image_facts.py`: one parser test per function on literal fixture strings from the Interfaces block; `test_collect_keys_exact` (a fake `run` returning canned outputs; assert `set(facts) == {…8 keys…}`); `test_missing_tool_aborts_without_writing` (fake `run` raising `FileNotFoundError` for `nc-config` → CLI rc 1, out file absent); `test_write_facts_format` (sorted `key=value`, newline-terminated).
- [ ] **Step 2: Run** `python docker/ctsm-ci-derecho-gnu/test_image_facts.py -v` → FAIL.
- [ ] **Step 3: Implement `image_facts.py`.**
- [ ] **Step 4: Run** → PASS.
- [ ] **Step 5: Write the failing checker tests** in `test_check_derecho_versions.py`: add `BASE_IMAGE_FACTS` text matching `BASE_DOCKERFILE`'s ARGs, a `make_image_facts` switch in the fixture, and class `TestImageFacts`: baseline adds zero failures; each of the five version keys perturbed → exactly one `❌` naming that key; `esmf_comm=mpich` → one `❌`; stackfacts ≠ nf-config → one `❌`; a deleted line → one `❌` with `no "hdf5" line`; a nonexistent path → `FileNotFoundError` (in `TestSetupProblemsRaise`). `TestExitCodes` gains a case running the CLI with `--image-facts` on a mismatching file → rc 1 and `FAILED`.
- [ ] **Step 6: Run** `python docker/ctsm-ci-derecho-gnu/test_check_derecho_versions.py -v` → new tests FAIL.
- [ ] **Step 7: Implement** in the checker; append the image lines after the `CHECKS` loop and before the module tree. Update the module docstring's subject paragraph (`:4-11`) to say the image's as-built facts are compared too when `--image-facts` is given.
- [ ] **Step 8: Run** the checker tests → PASS. Run the real checker on the login node with `--skip-module-tree` (no facts file) → unchanged `OK:` output.
- [ ] **Step 9: Dockerfile comment** near `:836`: the reader is `image_facts.py` inside a `container:` job. README "Bumping versions": add the image-facts paragraph.
- [ ] **Step 10: Commit** — `drift check: read the image's as-built facts and compare them to the recipe`. Then stop and ask the user to review.

---

### Task 5: Live netCDF dependency versions from derecho's `spack.lock`

**Files:**
- Modify: `docker/ctsm-ci-derecho-gnu/check-derecho-versions.py` (`check_module_tree`, new functions, `run_checks`, constants `:85-108`, docstring `:33-36`)
- Modify: `docker/ctsm-ci-derecho-gnu/test_check_derecho_versions.py`
- Modify: `docker/ctsm-ci-derecho-gnu/derecho-versions.ini` (header `:5-29`: the scheduled run now re-measures live; the ini is for the PR gate and is checked against live)
- Modify: `docker/ctsm-ci-derecho-gnu/README.md` ("Bumping versions": HDF5/netCDF-Fortran are read live on the schedule)

**Interfaces:**
- Consumes: nothing from Task 4 (independent).
- Produces: `check_module_tree(config, roots) -> (ok, lines, resolved)` where `resolved: dict[str, str]` maps `"<name>/<version>"` to the modulefile path found (only resolved pins; the existing two call sites and tests are updated for the third element). `DERECHO_APPS_ROOT = "/glade/u/apps/derecho"`; `SYSTEM_MODULE_ROOT = os.path.join(DERECHO_APPS_ROOT, "modules")` (same value as today). `spack_lock_path(apps_root, env_version) -> str` = `<apps_root>/<env_version>/envs/public/spack.lock`. `netcdf_spec_hash_prefix(lua_path) -> str` — basename of the path in `setenv("NCAR_ROOT_NETCDF", "…")`; raises `ValueError` naming the file if the line is absent. `spack_closure(lock: dict, hash_prefix: str) -> (root: tuple[str,str], versions: dict[str, set[str]])` — unique prefix match over `lock["concrete_specs"]` (`ValueError` on 0 or >1 matches, message naming the prefix and the count), depth-first over `dependencies[].hash`, `versions[name]` collects `version` strings. `LIVE_NETCDF_DEPS = [("hdf5", "HDF5_VERSION", ("snapshot", "hdf5")), ("netcdf-fortran", "NETCDF_FORTRAN_VERSION", ("snapshot", "netcdf_fortran"))]`. `check_live_netcdf_deps(args, snap, config, resolved, env_version, apps_root) -> (ok, lines)`: resolves `netcdf/<pinned>` from `resolved` (⚠️ not checked if the module did not resolve — the tree check already printed the ❌), reads the lock, asserts the root spec is `("netcdf", <pinned>)` and that `versions["netcdf-c"] == {<pinned>}`, then per entry: exactly one version in the closure (else ❌ listing them), equal to the ARG (`✅ live hdf5 <v> matches ARG HDF5_VERSION` / `❌ …`), and equal to the ini snapshot (`❌ [snapshot] hdf5 = <s> is stale: netcdf/<pinned> on derecho today pairs with hdf5 <v>. Re-measure and update derecho-versions.ini`). `run_checks(..., apps_root=DERECHO_APPS_ROOT)`; runs only when `skip_module_tree` is false.

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what earlier tasks landed; re-read this task for anything they or new measurements made stale, and amend it; collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Write the failing tests.** Fixture additions: the fake system tree's `netcdf/4.9.3.lua` gets the `setenv("NCAR_ROOT_NETCDF", "<tmp>/25.10/spack/opt/spack/netcdf/4.9.3/gcc/14.3.0/s6hn")` line; a fake `<tmp>/25.10/envs/public/spack.lock` written from a dict with specs `s6hn…` (netcdf 4.9.3; deps hdf5, netcdf-c, netcdf-fortran, gcc-runtime), `netcdf-c` depending on `hdf5` again (closure revisits), and an unrelated `zzzz…` spec. `run_checks` is called with `apps_root=<tmp>`. Class `TestLiveNetcdfDeps`: baseline → two new `✅ live` lines, zero `❌`; `test_live_hdf5_mismatch_fails` (lock hdf5 → 1.14.5 → exactly one `❌` naming `HDF5_VERSION`… plus the stale-snapshot line, so assert both and that they are distinct); `test_live_netcdf_fortran_mismatch_fails`; `test_stale_snapshot_reported_against_live` (ini `hdf5 = 1.12.2`, ARG `HDF5_VERSION=1.14.6`, lock 1.14.6 → exactly one `❌` containing `[snapshot] hdf5` and `stale`); `test_ambiguous_hash_prefix_is_a_finding` (second key starting `s6hn` → `❌` mentioning `2 specs`); `test_unknown_hash_prefix_is_a_finding`; `test_lua_without_ncar_root_is_a_finding`; `test_two_hdf5_versions_in_closure_fails`; `test_root_spec_is_not_pinned_netcdf_fails` (lock root version 4.9.2); `test_skip_module_tree_skips_live_deps` (no live lines, no lock read — point `apps_root` at a nonexistent dir). Existing `TestModuleTree` tests unpack the third return value.
- [ ] **Step 2: Run** → FAIL.
- [ ] **Step 3: Implement.** Findings inside `check_live_netcdf_deps` are `❌` lines, not exceptions (`ValueError` from the helpers is caught there and reported); a missing lock file is also a `❌` naming the path (an unreadable live source is drift-check failure by design, like a missing module root).
- [ ] **Step 4: Run** the checker tests → PASS.
- [ ] **Step 5: Run the real checker on the login node without `--skip-module-tree`** (the login node sees `/glade/u/apps`): expect `✅ live hdf5 1.14.6 …`, `✅ live netcdf-fortran 4.6.2 …`, and `OK:`. Paste the two lines into the commit body.
- [ ] **Step 6: Docs.** Checker docstring `:33-36` (the hdf5-repointing gap is now closed on the scheduled path); ini header; README.
- [ ] **Step 7: Commit** — `drift check: read HDF5 and netCDF-Fortran live from derecho's spack.lock`. Then stop and ask the user to review.

---

### Task 6: `derecho-drift.yml` — the scheduled monitor, the issue, the heartbeat

**Files:**
- Create: `.github/workflows/derecho-drift.yml`
- Create: `.github/workflows/derecho-drift.issue_template.md`
- Create: `.github/workflows/derecho-drift-silent.issue_template.md`
- Create: `docker/ctsm-ci-derecho-gnu/drift_heartbeat.py`
- Create: `docker/ctsm-ci-derecho-gnu/test_drift_heartbeat.py`
- Modify: `docker/ctsm-ci-derecho-gnu/README.md` ("Bumping versions": the monitor, its image pin, the issue titles, the heartbeat) and `NEXT_STEPS.md` item 5 (tick what landed)

**Interfaces:**
- Consumes: `image_facts.py --out` (Task 4); checker `--image-facts` (Task 4) and live reads (Task 5).
- Produces: `drift_heartbeat.py --runs-json PATH --max-age-days N` — `PATH` is the output of `gh run list --workflow derecho-drift.yml --event schedule --status completed --limit 1 --json createdAt,conclusion,url`; exit 0 and `✅ last completed scheduled run: <url> (<age> days ago, <conclusion>)` when `createdAt` is within N days; exit 0 and `⚠️ no completed scheduled run yet (first run)` on an empty list; exit 1 and `❌ no scheduled run has completed in <age> days (last: <url>). gha-runner-ctsm is probably offline or the job is queued; check the runner before trusting "no drift".` otherwise; exit 2 on unparseable JSON. `newest_completed_age_days(runs: list[dict], now: datetime) -> float | None` is the tested function.

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what earlier tasks landed; re-read this task for anything they or new measurements made stale, and amend it; collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Write the failing heartbeat tests**: empty list → `None`; fresh run → `< 1`; 12-day-old run → `≈ 12`; CLI rc 0/0/1 for those; malformed JSON → rc 2.
- [ ] **Step 2: Run** → FAIL. **Step 3: Implement.** **Step 4: Run** → PASS.
- [ ] **Step 5: Write the workflow.** `name: Derecho drift check`. Triggers: `schedule: - cron: '23 13 * * 1'` (moved from the PR gate); `workflow_dispatch: inputs: open_issue: {type: boolean, default: false, description: "Open/update the drift issue even though this is not a scheduled run"}`; `pull_request: paths: ['.github/workflows/derecho-drift.yml']` as the branch-test trigger, with the same explanatory comment as `probe-derecho-modules.yml:37-46`. `permissions: {contents: read, issues: write, actions: read}`. Jobs:
  - `image-facts`: `runs-on: gha-runner-ctsm`; `container: image: ghcr.io/escomp/ctsm/ctsm-ci-derecho-gnu:20261015` with the pin comment copied from `cirrus-testing.yml:50-61` plus "bump both files together"; `timeout-minutes: 30`; steps: checkout; `python3 docker/ctsm-ci-derecho-gnu/image_facts.py --out "$RUNNER_TEMP/image-facts.txt"`; `cat` it; upload-artifact `image-facts`.
  - `check`: `runs-on: gha-runner-ctsm` (no container); `needs: image-facts`; `timeout-minutes: 30`; checkout with `fetch-tags: true`; `git describe --tags --always` printed under the heading `CTSM ref`; setup-python 3.x; `bin/git-fleximod update ccs_config`; download-artifact; run `test_check_derecho_versions.py -v`, `test_image_facts.py -v`; then `set -o pipefail; python docker/ctsm-ci-derecho-gnu/check-derecho-versions.py --image-facts "$RUNNER_TEMP/image-facts.txt" | tee "$RUNNER_TEMP/drift-report.txt"`; upload-artifact `drift-report` with `if: always()`.
  - `report`: `runs-on: ubuntu-latest`; `needs: check`; `if: failure() && (github.event_name == 'schedule' || inputs.open_issue)`; checkout; download-artifact `drift-report` with `continue-on-error: true`; a step that writes `DRIFT_REPORT` to `$GITHUB_ENV` via a heredoc (`cat` of the report, or `"(no report file; the image-facts or check job failed before producing one — see the run)"`) and `RUN_URL=${{ github.server_url }}/${{ github.repository }}/actions/runs/${{ github.run_id }}`; `JasonEtco/create-an-issue@1b14a70e4d8dc185e5cc76d3bec9eab20257b2c5` with `filename: .github/workflows/derecho-drift.issue_template.md`, `update_existing: true`, `search_existing: open`, `env: GITHUB_TOKEN`.
  - `heartbeat`: `runs-on: ubuntu-latest`; `if: github.event_name == 'schedule'`; no `needs`; checkout; `gh run list --workflow derecho-drift.yml --event schedule --status completed --limit 1 --json createdAt,conclusion,url > "$RUNNER_TEMP/runs.json"` with `GH_TOKEN: ${{ secrets.GITHUB_TOKEN }}`; `python docker/ctsm-ci-derecho-gnu/drift_heartbeat.py --runs-json "$RUNNER_TEMP/runs.json" --max-age-days 10`; on failure (`if: failure()` step) the same action with `derecho-drift-silent.issue_template.md`.
- [ ] **Step 6: Templates.** `derecho-drift.issue_template.md`: frontmatter `title: Drift check: container vs derecho mismatch`, `labels: bug`; body: `{{ date }} — the scheduled drift check failed. Run: {{ env.RUN_URL }}` then the report in a fenced block via `{{ env.DRIFT_REPORT }}`, then one line pointing at `docker/ctsm-ci-derecho-gnu/README.md` "Bumping versions". `derecho-drift-silent.issue_template.md`: `title: Drift check: no scheduled run has completed in 10+ days`, body naming the runner and the run URL. Confirm from the action's README at the pinned SHA that `{{ env.X }}` is supported; if not, replace the two issue steps with `gh issue list --search "<title> in:title" --state open --json number` + `gh issue create`/`gh issue edit --body-file` (same semantics: one open issue per title, body replaced).
- [ ] **Step 7: README and NEXT_STEPS** as listed; README's "Bumping versions" now has the second image pin to bump.
- [ ] **Step 8: Commit** — `ci: add derecho-drift.yml, the scheduled drift monitor with issue and heartbeat`. Body states which parts ran (see Step 9) and which only ran as unit tests.
- [ ] **Step 9: Exercise from the branch.** Open/update the PR so the `pull_request` self-trigger fires: `image-facts` and `check` must both run on `gha-runner-ctsm`; `check` must print the two `✅ live` lines and `✅ image …` lines; `report` and `heartbeat` are skipped by their `if:`. Then, once on the default branch, a `workflow_dispatch` with `open_issue: true` and a deliberately wrong ARG on a scratch branch is the way to see the issue get created and then updated; the heartbeat can only be seen on the first real Monday run. Record outcomes in NEXT_STEPS. Commit that record, then stop and ask the user to review.

---

### Task 7: Gate the image's couplings to ccs_config's container machine and cime's cprnc

Added 2026-10-26. Depends on Task 2 (the gate's shape, and the staleness guard reads the same `paths:` list); independent of CCPP-502. Can run any time after Task 2.

**Why:** since 2026-10-26 (NEXT_STEPS item 13) the image (a) bakes `docker/ctsm-ci-derecho-gnu/cime-machines/container/config_machines.xml`, a stopgap copy of `ccs_config/machines/container/config_machines.xml` plus one `<CCSM_CPRNC>` line, which CIME reads *instead of* ccs_config's file; and (b) builds cprnc at the Dockerfile's `CPRNC_VERSION`, which is meant to be the cprnc commit cime records for `CIME/non_py/cprnc`. Both are recipe-vs-checkout couplings, which is this gate's job. Today (a) is guarded only by a `diff` step in `cirrus-testing.yml`, which item 10 will make on-demand, and (b) by nothing. Neither drift breaks a case build, so neither would be noticed: a stale (a) masks ccs_config changes to the container machine; a stale (b) leaves the image's cprnc behind cime's.

**Files:**
- Modify: `docker/ctsm-ci-derecho-gnu/check-derecho-versions.py`, `docker/ctsm-ci-derecho-gnu/test_check_derecho_versions.py`
- Modify: `.github/workflows/derecho-version-check.yml`
- Modify: `.github/workflows/cirrus-testing.yml` (remove the `Check container machine override against ccs_config` step)
- Modify: `docker/ctsm-ci-derecho-gnu/README.md` ("cprnc is built into the image", "Bumping versions"), `docker/ctsm-ci-derecho-gnu/NEXT_STEPS.md` (item 13), and the header comment of `cime-machines/container/config_machines.xml` (which names `cirrus-testing.yml` as the guard)

**Interfaces:**
- `override_differences(override_path: str, upstream_path: str) -> list[str]` — strip from the override the leading comment block (a line `<!--` through a line `-->`) and the single `<CCSM_CPRNC>` line, then return `difflib.unified_diff` lines against the upstream file; empty list = match. Raise `ValueError` if the override has zero or more than one `<CCSM_CPRNC>` line (that is a malformed override, i.e. a setup problem). Report lines: `✅ container machine override matches ccs_config's except CCSM_CPRNC` / `❌ container machine override differs from ccs_config/machines/container/config_machines.xml beyond CCSM_CPRNC (first lines of the diff …): re-copy ccs_config's file, re-add the CCSM_CPRNC line, rebuild and republish the image (the copy is baked in).`
- `cime_cprnc_commit(cime_dir: str) -> str` — `git -C <cime_dir> ls-tree HEAD CIME/non_py/cprnc`, the gitlink SHA. Needs cime checked out, not cprnc.
- `resolve_cprnc_tag(tag: str) -> str` — `git ls-remote https://github.com/ESMCI/cprnc refs/tags/<tag> refs/tags/<tag>^{}`, the peeled commit if present, else the plain one; `ValueError` if the tag does not exist. Network failure propagates (setup problem, per the exit-code constraint).
- `run_checks` gains `override=…`, `upstream_container_xml=…`, `cime_dir=…` and `resolve_tag=resolve_cprnc_tag` (injectable for tests). Report lines: `✅ CPRNC_VERSION <tag> is the cprnc commit cime records (<sha7>)` / `❌ CPRNC_VERSION <tag> is <sha7> but cime records cprnc <sha7>: set CPRNC_VERSION to cime's cprnc tag, rebuild and republish.`

- [ ] **Step 0 (orchestrator): Reconcile and stop.** Review what Tasks 1–6 landed; re-read this task against the current `cime-machines/` file, `CPRNC_VERSION`, cime's cprnc pin and `cirrus-testing.yml`, and amend anything stale (in particular: whether ccs_config has since gained `CCSM_CPRNC` upstream, in which case the override and half this task go away). Collect questions. **Stop for the user's approval before anything else in this task.**
- [ ] **Step 1: Write the failing tests** in `test_check_derecho_versions.py`, using temp files and a temp git repo:
  - `test_override_matches`: override = upstream text with a header comment and a `<CCSM_CPRNC>` line inserted → `[]`.
  - `test_override_extra_difference_fails`: additionally change `<GMAKE_J>` → non-empty, and the report line starts with `❌ container machine override`.
  - `test_override_without_cprnc_line_is_setup_error` and `test_override_with_two_cprnc_lines_is_setup_error` → `ValueError`.
  - `test_cprnc_pin_matches` / `test_cprnc_pin_mismatch_fails`: a temp repo whose index carries a gitlink (`git update-index --add --cacheinfo 160000,<sha>,CIME/non_py/cprnc`, then commit) and a fake `resolve_tag` returning the same / a different SHA.
  - `test_real_override_matches_checkout`: the real override against the real `ccs_config` file (skip if ccs_config is not populated).
- [ ] **Step 2: Run** `python docker/ctsm-ci-derecho-gnu/test_check_derecho_versions.py -v` → new tests FAIL.
- [ ] **Step 3: Implement** in the checker (stdlib only: `difflib`, `subprocess`). Update the module docstring's subject paragraph to say the gate also covers these two couplings.
- [ ] **Step 4: Run** the checker tests → PASS. Run the real checker on the login node with `--skip-module-tree` → the two new `✅` lines and `OK:` (needs `cime` populated, which a normal checkout is).
- [ ] **Step 5: Edit `derecho-version-check.yml`.** Add `docker/ctsm-ci-derecho-gnu/cime-machines/container/config_machines.xml` and `cime` to `pull_request: paths:` (so the staleness guard watches them too), with a comment in the style of the existing `ccs_config` one; add a step `bin/git-fleximod update cime` next to the `ccs_config` one.
- [ ] **Step 6: Remove the `diff` step from `cirrus-testing.yml`** and point the override file's header comment at the gate instead.
- [ ] **Step 7: Docs.** README "cprnc is built into the image" and "Bumping versions" (`CPRNC_VERSION` is now checked against cime's pin); NEXT_STEPS item 13.
- [ ] **Step 8: Commit** — `drift check: gate the container machine override and the cprnc pin`. Body states the workflow change is not CI-verified until a PR touching a watched path runs it. Then stop and ask the user to review.

---

### Appendix A (only if CCPP-502 is refused): cron-written derecho snapshot

Build this only after the ticket is refused; expand into full steps then. Shape, so the interfaces above are not designed against it:

- `check-derecho-versions.py --write-derecho-facts PATH` writes, atomically (`tmp` + `os.replace`), JSON `{"schema": 1, "generated_at": "<UTC ISO-8601>", "apps_root": …, "modulefiles": [<every .lua path relative to each root under environment/*/, <env>/Core/**, <env>/gcc/<gccver>/** for every env and gcc version present, both roots>], "netcdf": {"<relative lua path>": {"root": ["netcdf", "<ver>"], "closure": {"hdf5": […], "netcdf-c": […], "netcdf-fortran": […]}}}}` — pin-independent so the cron needs no CTSM checkout and never talks to GitHub.
- `--derecho-facts PATH --max-age-days N` makes `check_module_tree` test membership in `modulefiles` instead of `os.path.isfile`, and `check_live_netcdf_deps` read the `netcdf` map; a report older than N days is `❌ derecho snapshot is <age> days old (limit N): the cron on <host> has stopped`.
- The user installs the cron (daily, on a Casper or derecho login node — their choice) writing to a path they choose under `/glade/campaign/cesm/cesmdata/cseg/`, not inside `inputdata`; `derecho-drift.yml`'s `check` job passes `--derecho-facts <that path> --max-age-days 3` instead of reading live. Tests mirror Task 5's with the snapshot as the source, plus staleness and schema-version cases.

---

## Verification (end to end)

- Unit: `for t in docker/ctsm-ci-derecho-gnu/test_*.py; do python "$t" -v || exit 1; done` — all green on the login node.
- Checker live: `python docker/ctsm-ci-derecho-gnu/check-derecho-versions.py` (no flags) on the login node prints the two `✅ live` lines and `OK:`.
- PR gate: a PR touching `docker/ctsm-ci-derecho-gnu/Dockerfile` runs `derecho-version-check.yml` on `ubuntu-latest` and passes the staleness guard; after a watched file lands on the base branch and before the PR is updated, a re-run of the gate fails with the `❌ stale` message.
- Monitor: the `pull_request` self-trigger run of `derecho-drift.yml` (Task 6 Step 9) shows `image-facts` reading all eight keys from `:20261015` and `check` reporting `OK:` with image and live lines.
- Executable assertion (user, Casper): Task 3 Step 8.
- Override and cprnc pin (Task 7): a PR touching `cime-machines/container/config_machines.xml` runs the gate and prints both new `✅` lines; a scratch commit that edits a second line of the override, or sets a wrong `CPRNC_VERSION`, turns each into `❌`.
- Manual derecho-binary validation (user, derecho; once): build a current `_Mmpi-serial` `derecho_gnu` case and confirm its `readelf -d cesm.exe` RPATH names `mpi-serial/2.5.3`, `parallelio/2.6.8`, `esmf/8.9.1`, and that `nf-config --version` under the netcdf prefix prints `4.6.2` — i.e. that the `spack.lock` read in Task 5 reports what a derecho binary actually links. Record the date in `DRIFT_CHECK_DESIGN.md`'s instrument map.
- Schedule itself: only observable after merge to the default branch; the first Monday run should show all four jobs, with `report` skipped and `heartbeat` printing `⚠️ no completed scheduled run yet`.

## Deliberately left out

- Comparing compression filters (instrument map's last row): the mapping between Spack dependency names in `spack.lock` and `nc-config --has-stdfilters` tokens is unverified; it needs a measurement of both sides before a pass/fail rule can be written. The container's set is already recorded in `netcdf_c_stdfilters` and now travels in the image-facts file, so the measurement is one `cat` away once the monitor runs.
- Item 10 (every-push triggers on `cirrus-testing.yml`) and item 13 (the red jobs): untouched, except that Task 3 adds a step to a job item 13 must fix.
