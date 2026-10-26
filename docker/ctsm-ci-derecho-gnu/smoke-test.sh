#!/usr/bin/env bash
# Smoke-test the built ctsm-ci-derecho-gnu image. Asserts the toolchain versions
# match derecho's gnu + mpi-serial stack, that the PIO a case build links is
# actually installed and advertised, and that a small netCDF Fortran program
# compiles, links (-lnetcdff -lnetcdf -llapack -lblas) and runs. Run after
# build-on-casper.sh tags the image.
#
# There is no MPI here to test: since NEXT_STEPS item 4b the image carries one
# serial stack and no MPICH.
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

# mpi-serial and pio are existence-only, here and below: there is no honest
# source in the image to compare their values against (see the Dockerfile's
# LABEL block). Their presence is still worth asserting -- item 5 reads the
# container's mpi-serial and PIO versions from exactly these two labels, so an
# image that lost them would leave the drift check with nothing to read.
for l in gcc netcdf-c netcdf-fortran hdf5 esmf mpi-serial pio; do
    if [ -z "$(lbl "$l")" ]; then
        echo "SMOKE FAIL: image has no '$l' label to check against" >&2
        exit 1
    fi
done

# The checks run INSIDE a fresh container so we exercise exactly the baked-in
# environment (non-login shell, image ENV only) that GitHub Actions sees.
podman run --rm -i \
    -e "WANT_GCC=$(lbl gcc)" \
    -e "WANT_NETCDF_C=$(lbl netcdf-c)" \
    -e "WANT_NETCDF_F=$(lbl netcdf-fortran)" \
    -e "WANT_HDF5=$(lbl hdf5)" \
    -e "WANT_ESMF=$(lbl esmf)" \
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

# CC/CXX/FC decide which toolchain a CMake build inside the container picks
# when it is not told: cmake reads these or falls back to searching for
# `cc`/`f95` -- names that exist only in /usr/bin, as AlmaLinux's bootstrap
# GCC. That is how cprnc v1.1.2 linked the image's GCC-built static HDF5 with
# the wrong libgcc and died on undefined __truncsfhf2. cprnc v1.1.5
# (cime6.5.5) takes its compilers from the case's Macros.cmake instead, so
# this now guards other CMake builds rather than cprnc's. Checking the version each one REPORTS, not just that the path
# exists, is what makes this a real check: the Dockerfile already asserts
# existence at build time, which says nothing about the image in front of us.
echo "### CC/CXX/FC (what a cmake build in here will use)"
for v in CC CXX FC; do
    eval "p=\${$v:-}"
    [ -n "$p" ]                                     || fail "$v is unset; cmake would fall back to /usr/bin"
    case "$p" in /usr/bin/*) fail "$v=$p is the distro bootstrap compiler" ;; esac
    echo "$v=$p"
    check "$v" "$WANT_GCC" "$("$p" --version | head -1)"
done

echo "### netcdf-c";   nc-config --version; echo "prefix=$(nc-config --prefix)"
check netcdf-c "$WANT_NETCDF_C" "$(nc-config --version)"
[ "$(nc-config --prefix)" = /usr/local ]            || fail "netcdf-c prefix != /usr/local"

echo "### netcdf-fortran"; nf-config --version
check netcdf-fortran "$WANT_NETCDF_F" "$(nf-config --version)"

# h5dump rather than h5cc: the tools are installed whatever the wrapper is
# called. Fall back to whichever wrapper exists if h5dump is not there, so this
# reports a real version or fails, never silently passes.
echo "### hdf5"
hdf5_out="$(h5dump --version 2>/dev/null)" || hdf5_out=""
[ -n "$hdf5_out" ] || hdf5_out="$(h5cc  -showconfig 2>/dev/null | grep -i 'HDF5 Version')"
[ -n "$hdf5_out" ] || fail "hdf5: no h5dump/h5cc available to report a version"
echo "$hdf5_out"
check hdf5 "$WANT_HDF5" "$hdf5_out"

echo "### esmf";       echo "ESMFMKFILE=${ESMFMKFILE:-<unset>}"
[ -n "${ESMFMKFILE:-}" ] && [ -f "$ESMFMKFILE" ]    || fail "ESMFMKFILE missing"
# esmf.mk records the version as a line-start assignment; the mk file is the
# thing CTSM actually consumes, so read the version from it rather than from a
# directory name that happens to contain digits.
esmf_out="$(grep '^ESMF_VERSION_STRING=' "$ESMFMKFILE")" \
    || fail "esmf: no ESMF_VERSION_STRING in $ESMFMKFILE"
echo "$esmf_out"
check esmf "$WANT_ESMF" "$esmf_out"
# Since item 4b the default is the mpiuni build, because it is the only build.
# Catching a tarball that still carries an ESMF_COMM=mpich default here is the
# cheapest way to notice that the image in front of us predates the switch.
grep -q 'ESMF_COMM=mpiuni' "$ESMFMKFILE" \
    || fail "$ESMFMKFILE is not an ESMF_COMM=mpiuni build; this image predates item 4b"

# The same file reached the other way: through the hardcoded path the CIME
# macro selects when MPILIB=mpi-serial. The Dockerfile asserts both at build
# time, but that says nothing about the image in front of us now -- a
# mis-tagged or older tarball passes an assertion that ran in some other build.
echo "### esmf, mpiuni flavor (from the baked-in CIME macro)"
macro=/opt/ctsm-container/cime-macros/gnu_container.cmake
[ -f "$macro" ] || fail "no CIME macro at $macro"
serial_mk="$(sed -n 's/^ *set(ESMFMKFILE "\([^"]*\)").*/\1/p' "$macro" | tail -1)"
[ -n "$serial_mk" ] || fail "no set(ESMFMKFILE ...) in $macro"
echo "macro ESMFMKFILE=$serial_mk"
[ -f "$serial_mk" ] || fail "macro points ESMFMKFILE at $serial_mk, which does not exist"
grep -q 'ESMF_COMM=mpiuni' "$serial_mk" \
    || fail "$serial_mk is not an ESMF_COMM=mpiuni build"
serial_out="$(grep '^ESMF_VERSION_STRING=' "$serial_mk")" \
    || fail "esmf mpiuni: no ESMF_VERSION_STRING in $serial_mk"
echo "$serial_out"
check "esmf mpiuni" "$WANT_ESMF" "$serial_out"

# The PIO a case build links. These three files are what makes the switch real:
# cime/CIME/Tools/Makefile:667 names libpiof.a and libpioc.a by path, and the
# CTSM sources USE pio. Without them a case build fails at the link, having
# already skipped CTSM's libraries/parallelio submodule.
#
# EXISTENCE ONLY, NO VERSION COMPARISON, for PIO or for mpi-serial. Neither the
# `pio` nor the `mpi-serial` label can be checked against what the image
# reports: nothing in the image reports an mpi-serial version at all, and PIO's
# own lib/libpio.settings reports 2.6.6 for every tag from pio2_6_6 through
# pio2_6_8 because upstream stopped bumping the in-tree version. Comparing
# either label to itself would prove nothing; comparing `pio` to
# libpio.settings would fail a correct image. See the LABEL block in the
# Dockerfile.
echo "### PIO install"
for f in /usr/local/lib/libpiof.a /usr/local/lib/libpioc.a /usr/local/include/pio.mod; do
    [ -f "$f" ] || fail "PIO: $f missing"
    echo "  $f"
done

# The environment half of that same switch, and it has to be the ENVIRONMENT:
# cime/CIME/build_scripts/buildlib.pio:82-91 reads PIO_VERSION_MAJOR from
# os.environ, which no macro file reaches. Unset PIO_TYPENAME_VALID_VALUES does
# not fall back to a default either -- line 90 hands None to
# _set_pio_valid_values and it raises TypeError.
echo "### PIO environment"
for v in PIO PIO_VERSION_MAJOR PIO_TYPENAME_VALID_VALUES; do
    eval "val=\${$v:-}"
    [ -n "$val" ] || fail "$v is not set in the image environment"
    echo "  $v=$val"
done
[ "${PIO}" = /usr/local ]                || fail "PIO=${PIO}, expected /usr/local"
[ "${PIO_VERSION_MAJOR}" = 2 ]           || fail "PIO_VERSION_MAJOR=${PIO_VERSION_MAJOR}, expected 2"
[ "${PIO_TYPENAME_VALID_VALUES}" = netcdf ] \
    || fail "PIO_TYPENAME_VALID_VALUES=${PIO_TYPENAME_VALID_VALUES}, expected netcdf"

# Two facts about the built stack that no label can carry and no container
# cesm.exe reports: see the Dockerfile for why they are a file. Item 5's drift
# check reads this; check here that it exists and is populated, so a missing or
# empty field is caught at build time rather than read as agreement.
echo "### stack facts"
facts=/opt/ctsm-container/stack-facts.txt
[ -f "$facts" ] || fail "no $facts in the image"
cat "$facts"
# Both KEYS must be present. Only netcdf_fortran_version must have a VALUE:
# an empty netcdf_c_stdfilters is a legitimate answer meaning "no standard
# filters", which is a real difference from derecho for item 5 to report, not a
# broken image. See the Dockerfile for the same distinction at build time.
for k in netcdf_fortran_version netcdf_c_stdfilters; do
    grep -q "^${k}=" "$facts" || fail "stack facts: no ${k}= line in $facts"
done
[ -n "$(sed -n 's/^netcdf_fortran_version=//p' "$facts")" ] \
    || fail "stack facts: netcdf_fortran_version is empty in $facts"
# The recorded netCDF-Fortran version must be the one the image advertises;
# otherwise the drift check and the version check describe different images.
check "stack facts netcdf-fortran" "$WANT_NETCDF_F" \
    "$(sed -n 's/^netcdf_fortran_version=//p' "$facts")"

echo "### perl XML::LibXML"
perl -MXML::LibXML -e 'print "XML::LibXML OK\n"'    || fail "perl XML::LibXML"

# No MPI here: gfortran, not mpifort, and no mpiexec. -lnetcdff -lnetcdf is the
# link ccs_config/machines/container/container.cmake appends to SLIBS.
echo "### netCDF Fortran build/run"
cat > /tmp/hello.f90 <<'EOF'
program hello
  use netcdf
  implicit none
  print '(A,A)', 'netcdf=', trim(nf90_inq_libvers())
end program hello
EOF
gfortran -I/usr/local/include /tmp/hello.f90 -o /tmp/hello \
    -L/usr/local/lib -lnetcdff -lnetcdf -llapack -lblas \
    $(nc-config --libs --static | tr ' ' '\n' | grep '^-l' | tr '\n' ' ') \
    || fail "compile/link"
/tmp/hello                                          || fail "run"

echo "ALL SMOKE TESTS PASSED"
INSIDE
