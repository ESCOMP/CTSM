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

# The checks run INSIDE a fresh container so we exercise exactly the baked-in
# environment (non-login shell, image ENV only) that GitHub Actions sees.
podman run --rm -i "${image}" bash -s <<'INSIDE'
set -u
fail() { echo "SMOKE FAIL: $*" >&2; exit 1; }

echo "### gcc";        gcc --version | head -1
gcc --version | head -1 | grep -q '12\.2\.0'        || fail "gcc != 12.2.0"
echo "### gfortran";   gfortran --version | head -1
gfortran --version | head -1 | grep -q '12\.2\.0'   || fail "gfortran != 12.2.0"

echo "### mpich";      mpichversion | head -1
mpichversion | grep -q 'Version:[[:space:]]*3\.4\.3' || fail "mpich != 3.4.3"

echo "### netcdf-c";   nc-config --version; echo "prefix=$(nc-config --prefix)"
nc-config --version | grep -q '4\.9\.2'             || fail "netcdf-c != 4.9.2"
[ "$(nc-config --prefix)" = /usr/local ]            || fail "netcdf-c prefix != /usr/local"

echo "### netcdf-fortran"; nf-config --version
nf-config --version | grep -q '4\.6\.1'             || fail "netcdf-fortran != 4.6.1"

echo "### pnetcdf";    pnetcdf-config --version
pnetcdf-config --version | grep -q '1\.12\.3'       || fail "pnetcdf != 1.12.3"

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
