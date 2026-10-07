#!/usr/bin/env bash
# Smoke-test the built ctsm-ci-derecho-gnu image. Asserts the toolchain versions
# match derecho's gnu stack and that a small MPI + netCDF Fortran program
# compiles, links (-lnetcdff -lnetcdf -llapack -lblas) and runs under
# mpiexec. Run after build-on-casper.sh tags the image.
#
# Usage: docker/ctsm-ci-derecho-gnu/smoke-test.sh
#   IMAGE_TAG overrides the image (default localhost/ctsm-ci-derecho-gnu:dev)
set -eo pipefail

module load podman 2>/dev/null || true

image="${IMAGE_TAG:-localhost/ctsm-ci-derecho-gnu:dev}"
echo "Smoke-testing ${image}"

# Expected versions come from the image's OWN labels, which the Dockerfile sets
# from the version ARGs -- never from a copy kept here, which silently goes
# stale on every version bump and then fails a correct image. Comparing the
# labels against what the tools actually report is a real check: it catches a
# build that used a cached layer, or that quietly produced something other than
# the ARG asked for.
lbl() {
    local v
    v="$(podman inspect --format "{{ index .Labels \"$1\" }}" "${image}" 2>/dev/null)" \
        || v=""
    if [ -z "$v" ] || [ "$v" = "<no value>" ]; then
        v="$(podman inspect --format "{{ index .Config.Labels \"$1\" }}" "${image}" 2>/dev/null)" \
            || v=""
    fi
    [ "$v" = "<no value>" ] && v=""
    printf '%s' "$v"
}

for l in gcc mpich netcdf-c netcdf-fortran pnetcdf; do
    if [ -z "$(lbl "$l")" ]; then
        echo "SMOKE FAIL: image has no '$l' label to check against" >&2
        exit 1
    fi
done

# The checks run INSIDE a fresh container so we exercise exactly the baked-in
# environment (non-login shell, image ENV only) that GitHub Actions sees.
podman run --rm -i \
    -e "WANT_GCC=$(lbl gcc)" \
    -e "WANT_MPICH=$(lbl mpich)" \
    -e "WANT_NETCDF_C=$(lbl netcdf-c)" \
    -e "WANT_NETCDF_F=$(lbl netcdf-fortran)" \
    -e "WANT_PNETCDF=$(lbl pnetcdf)" \
    "${image}" bash -s <<'INSIDE'
set -u
fail() { echo "SMOKE FAIL: $*" >&2; exit 1; }

# Literal (-F) match of the expected version anywhere in the tool's own output.
check() {
    local what="$1" want="$2" got="$3"
    printf '%s' "$got" | grep -Fq -- "$want" \
        || fail "$what: expected $want, got: $got"
}

echo "### gcc";        gcc --version | head -1
check gcc "$WANT_GCC" "$(gcc --version | head -1)"
echo "### gfortran";   gfortran --version | head -1
check gfortran "$WANT_GCC" "$(gfortran --version | head -1)"

echo "### mpich";      mpichversion | head -1
check mpich "$WANT_MPICH" "$(mpichversion)"

echo "### netcdf-c";   nc-config --version; echo "prefix=$(nc-config --prefix)"
check netcdf-c "$WANT_NETCDF_C" "$(nc-config --version)"
[ "$(nc-config --prefix)" = /usr/local ]            || fail "netcdf-c prefix != /usr/local"

echo "### netcdf-fortran"; nf-config --version
check netcdf-fortran "$WANT_NETCDF_F" "$(nf-config --version)"

echo "### pnetcdf";    pnetcdf-config --version
check pnetcdf "$WANT_PNETCDF" "$(pnetcdf-config --version)"

echo "### esmf";       echo "ESMFMKFILE=${ESMFMKFILE:-<unset>}"
[ -n "${ESMFMKFILE:-}" ] && [ -f "$ESMFMKFILE" ]    || fail "ESMFMKFILE missing"

echo "### perl XML::LibXML"
perl -MXML::LibXML -e 'print "XML::LibXML OK\n"'    || fail "perl XML::LibXML"

echo "### MPI + netCDF Fortran build/run"
cat > /tmp/hello.f90 <<'EOF'
program hello
  use mpi
  use netcdf
  implicit none
  integer :: ierr, rank, nprocs
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  if (rank == 0) print '(A,I0,A,A)', 'ranks=', nprocs, ' netcdf=', trim(nf90_inq_libvers())
  call MPI_Finalize(ierr)
end program hello
EOF
mpifort -I/usr/local/include /tmp/hello.f90 -o /tmp/hello \
    -L/usr/local/lib -lnetcdff -lnetcdf -llapack -lblas || fail "compile/link"
mpiexec -n 2 /tmp/hello                                 || fail "mpiexec run"

echo "ALL SMOKE TESTS PASSED"
INSIDE
