# Source this (don't execute) to set up an environment with nvfortran available.
# Intended use:   . ./env.sh && nvfortran ...   or   . ./env.sh && make
# Safe to source multiple times.

if [ -z "${LMOD_CMD:-}" ]; then
    # shellcheck disable=SC1091
    . /etc/profile.d/lmod.sh
fi

module load nvhpc >/dev/null 2>&1 || {
    echo "env.sh: 'module load nvhpc' failed" >&2
    return 1 2>/dev/null || exit 1
}

if ! command -v nvfortran >/dev/null 2>&1; then
    echo "env.sh: nvfortran still not on PATH after module load nvhpc" >&2
    return 1 2>/dev/null || exit 1
fi
