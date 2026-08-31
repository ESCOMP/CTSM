#!/usr/bin/env bash
# Shared plumbing for the ctsm-ci-derecho-gnu container wrappers,
# run-case-in-container.sh and run-test-in-container.sh.
#
# Nothing here runs at source time: this file only defines functions, one set
# for each side of the container boundary.
#
#   ctsm_host_setup       HOST side, before podman. Validates and creates the
#                         persistent host directories, then builds the podman
#                         argument array.
#   ctsm_podman_run       HOST side. Runs the container (or, under DRY_RUN,
#                         prints what it would have run).
#   ctsm_has_arg          HOST side. Did the user already pass this flag?
#   ctsm_container_setup  INSIDE the container, first thing. Puts
#                         gnu_container.cmake where CIME will look for it.
#
# The repo is mounted at /ctsm, so the in-container script sources this very
# same file. That is why it must stay free of top-level side effects -- a
# `module load podman` at file scope would fire inside the container too.
#
# WHY THE PERSISTENT DIRECTORIES ARE CONTAINER-SPECIFIC
# Cases built here are built by a different toolchain than a native derecho or
# casper build (see README), so their bld/ and run/ directories must not be
# mixed in with native ones. Hence the `_devcontainer` suffix on both host
# paths: $HOME/cases_devcontainer and $SCRATCH/cases_devcontainer.

# --------------------------------------------------------------------------
# Host side
# --------------------------------------------------------------------------
# Sets for the caller:
#   CTSM_REPO         the CTSM checkout this script lives in
#   CTSM_IMAGE        image to run
#   CTSM_PODMAN_ARGS  array of `podman run` flags (mounts, workdir, env)
#   CTSM_INPUTDATA    resolved host inputdata path
#   CTSM_CASES_DIR    resolved host case directory
#   CTSM_SCRATCH_DIR  resolved host bld/run/archive directory
# The three resolved paths are exported so callers can name host-side paths in
# error messages without re-deriving the defaults (and drifting from them).
# Honors, from the environment:
#   IMAGE_TAG    image to run          (default localhost/ctsm-ci-derecho-gnu:dev)
#   INPUTDATA    host inputdata        (default the CESM tree on campaign)
#   CASES_DIR    host case dirs        (default $HOME/cases_devcontainer)
#   SCRATCH_DIR  host bld/run/archive  (default $SCRATCH/cases_devcontainer)
#   MAKE_J       build parallelism     (default 16)
#   MOUNT_OPTS   extra comma-separated mount options, e.g. "z" if SELinux
#                relabeling turns out to be needed. No leading colon.
#   INPUTDATA_SYMLINK_ROOT
#                tree that inputdata's absolute symlinks point into, mounted
#                read-only at the SAME path inside the container so they
#                resolve (default /glade/campaign). Set to the empty string to
#                skip it. See the mount block below for why this is needed.
ctsm_host_setup() {
    local here repo inputdata cases scratch suffix symlink_root

    here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    repo="$(cd "${here}/../.." && pwd)"

    module load podman 2>/dev/null || true

    inputdata="${INPUTDATA:-/glade/campaign/cesm/cesmdata/cseg/inputdata}"
    if [ ! -d "${inputdata}" ]; then
        echo "ERROR: inputdata directory not found: ${inputdata}" >&2
        echo "       Set INPUTDATA to override. On NCAR HPC this is normally" >&2
        echo "       /glade/campaign/cesm/cesmdata/cseg/inputdata." >&2
        return 1
    fi

    cases="${CASES_DIR:-${HOME}/cases_devcontainer}"

    if [ -n "${SCRATCH_DIR:-}" ]; then
        scratch="${SCRATCH_DIR}"
    elif [ -n "${SCRATCH:-}" ]; then
        scratch="${SCRATCH}/cases_devcontainer"
    else
        echo "ERROR: neither SCRATCH_DIR nor SCRATCH is set, so there is no" >&2
        echo "       place to put bld/run/archive. On NCAR HPC \$SCRATCH is" >&2
        echo "       normally /glade/derecho/scratch/\$USER." >&2
        return 1
    fi

    mkdir -p "${cases}" "${scratch}"

    # MOUNT_OPTS is appended after :rw / :ro, so it is given without a leading
    # colon (MOUNT_OPTS=z -> ":rw,z"), unlike the older run-unit-tests script.
    suffix="${MOUNT_OPTS:+,${MOUNT_OPTS}}"

    CTSM_REPO="${repo}"
    CTSM_IMAGE="${IMAGE_TAG:-localhost/ctsm-ci-derecho-gnu:dev}"
    CTSM_INPUTDATA="${inputdata}"
    CTSM_CASES_DIR="${cases}"
    CTSM_SCRATCH_DIR="${scratch}"

    # /opt/cesmdata/inputdata is DIN_LOC_ROOT for CIME's container machine
    # ($CESMDATAROOT/inputdata, and the image sets CESMDATAROOT=/opt/cesmdata).
    # The Dockerfile creates it empty, so it is a ready-made mount point and
    # no xmlchange or CESMDATAROOT override is needed. Read-only on purpose:
    # nothing in a run should ever write to the shared inputdata tree.
    #
    # /cases and /scratch are fixed container paths rather than $HOME-relative
    # ones because HOME differs by context (/root under podman, /github/home
    # in GitHub Actions); fixed paths keep these scripts portable to CI.
    CTSM_PODMAN_ARGS=(
        --rm -i
        -v "${repo}:/ctsm:rw${suffix}"
        -v "${inputdata}:/opt/cesmdata/inputdata:ro${suffix}"
        -v "${cases}:/cases:rw${suffix}"
        -v "${scratch}:/scratch:rw${suffix}"
        -w /cases
        -e "MAKE_J=${MAKE_J:-16}"
    )

    # WHY A SECOND, WIDER READ-ONLY MOUNT
    # Much of the CESM inputdata tree is not files but absolute symlinks into
    # sibling campaign collections -- e.g. the GSWP3 datm forcing under
    #   .../inputdata/atm/datm7/atm_forcing.datm7.GSWP3.0.5d.v1.c170516/TPHWL/
    # points at /glade/campaign/collections/gdex/data/d651077/cesmdata/...
    # A bind mount does not rewrite symlink targets, so inside the container
    # those absolute paths do not exist and every such file looks MISSING.
    # CIME then tries to download it, which the read-only mount refuses, and
    # the run dies with a wall of "wget failed ... No such file or directory".
    # Mounting the enclosing tree at its own path makes them resolve. Measured
    # on 2026-08-31: all 744 files a single-point IHistClm60Bgc test reported
    # missing were symlinks, and every one pointed under /glade/campaign.
    symlink_root="${INPUTDATA_SYMLINK_ROOT-/glade/campaign}"
    if [ -n "${symlink_root}" ] && [ -d "${symlink_root}" ]; then
        CTSM_PODMAN_ARGS+=(-v "${symlink_root}:${symlink_root}:ro${suffix}")
    fi

    echo "image      ${CTSM_IMAGE}"
    echo "repo       ${repo} -> /ctsm"
    echo "inputdata  ${inputdata} -> /opt/cesmdata/inputdata (read-only)"
    echo "cases      ${cases} -> /cases"
    echo "scratch    ${scratch} -> /scratch (bld, run, archive)"
    if [ -n "${symlink_root}" ] && [ -d "${symlink_root}" ]; then
        echo "symlinks   ${symlink_root} -> ${symlink_root} (read-only; inputdata symlink targets)"
    fi
}

# Print where the case, build, run and archive directories land on the HOST,
# so a run can be found afterwards without reverse-engineering the mounts.
# CIME puts bld/ and run/ at $CIME_OUTPUT_ROOT/$CASE/{bld,run}, and this
# wrapper points DOUT_S_ROOT at /scratch/archive/$CASE.
#   $1  host-side case directory
#   $2  case name (CIME's $CASE, i.e. the case directory's basename)
# create_test generates its case names from the test id, so the test wrapper
# passes placeholders rather than real names.
ctsm_print_host_dirs() {
    local host_case="$1" case_name="$2"
    echo "host-side directories:"
    echo "  case     ${host_case}"
    echo "  build    ${CTSM_SCRATCH_DIR}/${case_name}/bld"
    echo "  run      ${CTSM_SCRATCH_DIR}/${case_name}/run"
    echo "  archive  ${CTSM_SCRATCH_DIR}/archive/${case_name}"
}

# Is `--foo` (or `--foo=bar`) already present in the caller's arguments?
# Used to make our injected defaults yield to anything the user passed.
ctsm_has_arg() {
    local needle="$1"
    shift
    local arg
    for arg in "$@"; do
        [ "${arg}" = "${needle}" ] && return 0
        case "${arg}" in
            "${needle}="*) return 0 ;;
        esac
    done
    return 1
}

# Run the container. With DRY_RUN set to anything non-empty, print the command
# instead -- useful for checking the mounts and the assembled CIME argument
# list without a container or an allocation.
ctsm_podman_run() {
    if [ -n "${DRY_RUN:-}" ]; then
        printf 'DRY_RUN: podman run'
        printf ' %q' "${CTSM_PODMAN_ARGS[@]}" "${CTSM_IMAGE}" "$@"
        printf '\n'
        return 0
    fi
    podman run "${CTSM_PODMAN_ARGS[@]}" "${CTSM_IMAGE}" "$@"
}

# --------------------------------------------------------------------------
# Container side
# --------------------------------------------------------------------------
# CIME's copy_local_macros_to_dir() copies $HOME/.cime/*.cmake into the case.
# The image bakes a copy into /root/.cime, which is enough for a plain
# `podman run`, but HOME is overridden in GitHub Actions (/github/home), so
# place the file wherever HOME actually points. See the README section
# "Running CTSM's unit tests" and the header of cime-macros/gnu_container.cmake.
# The REPO's copy wins over the image's baked one whenever the repo is
# mounted. The image bakes a snapshot at build time, so without this a macro
# change appears to have no effect until the image is rebuilt -- which cost
# real debugging time once already. The image copy stays as the fallback for
# running without a mounted checkout.
ctsm_container_setup() {
    local repo_macro=/ctsm/docker/ctsm-ci-derecho-gnu/cime-macros/gnu_container.cmake
    local image_macro=/opt/ctsm-container/cime-macros/gnu_container.cmake

    # CIME snapshots each case directory with `git commit`. With no identity
    # configured the container has none, and a case build fills its log with
    # non-fatal but alarming
    #   ERROR: Command: 'git -C /cases/... commit ...' failed with error
    #          'Author identity unknown'
    # Supply one, without overriding anything the caller already set.
    : "${GIT_AUTHOR_NAME:=ctsm-container}"
    : "${GIT_AUTHOR_EMAIL:=ctsm-container@localhost}"
    : "${GIT_COMMITTER_NAME:=${GIT_AUTHOR_NAME}}"
    : "${GIT_COMMITTER_EMAIL:=${GIT_AUTHOR_EMAIL}}"
    export GIT_AUTHOR_NAME GIT_AUTHOR_EMAIL GIT_COMMITTER_NAME GIT_COMMITTER_EMAIL

    mkdir -p "${HOME}/.cime"
    if [ -f "${repo_macro}" ]; then
        echo "installing gnu_container.cmake into ${HOME}/.cime (from the mounted repo)"
        cp "${repo_macro}" "${HOME}/.cime/"
    elif [ ! -f "${HOME}/.cime/gnu_container.cmake" ]; then
        echo "installing gnu_container.cmake into ${HOME}/.cime (baked into the image)"
        cp "${image_macro}" "${HOME}/.cime/"
    fi
}
