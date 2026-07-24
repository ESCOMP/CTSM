# ctsm-ci-derecho-gnu: CTSM build/test container for GitHub Actions

This image lets CIME's `container` machine (`--machine container`, compiler
`gnu`; defined in `ccs_config/machines/container/`) compile CTSM cases in
GitHub Actions, e.g. the `simple-build-create_test` job in
[.github/workflows/cirrus-testing.yml](../../.github/workflows/cirrus-testing.yml).
Despite the name, it is **not** meant to run on the derecho HPC machine
itself — it *replicates* derecho's gnu software stack inside a container so
CTSM can be built and tested in CI. An image for running on derecho proper
would be Cray-based and separate.

It is built **from scratch** on AlmaLinux 9 (no vendor base image): every
component below is compiled or installed here so it matches derecho's gnu
stack wherever an open equivalent exists.

## Design goals

1. **Match derecho's gnu software stack** (what
   `ccs_config/machines/derecho/config_machines.xml` loads for gnu builds) as
   closely as a public container allows.
2. **Satisfy CIME's `container` machine unmodified**: netCDF/PnetCDF under
   `/usr/local`, LAPACK/BLAS linkable as `-llapack -lblas`, `ESMFMKFILE`,
   `USER`, and `CESMDATAROOT` preset.
3. **Work in GitHub Actions**: GHA runs steps in *non-login* shells, so
   profile scripts never fire. Everything needed is baked in as Docker `ENV`.

## Software stack vs. derecho gnu

The toolchain (GCC, MPICH) lives under `/opt`; libraries install under
`/usr/local` (where the `container` machine hard-codes `NETCDF_PATH` /
`PNETCDF_PATH`); ESMF is located via `ESMFMKFILE`.

| Component | derecho gnu | this image | note |
|---|---|---|---|
| GCC | 12.2.0 | 12.2.0 | exact match, built from source under `/opt/gcc` |
| MPI | cray-mpich 8.1.27 | MPICH 3.4.3 (ch4:ofi) | cray-mpich 8.x is MPICH-3.4-ABI-derived; Cray code is proprietary |
| HDF5 | hdf5-mpi/1.12.2 | 1.12.2 (parallel) | confirmed on derecho: loaded by netcdf-mpi/4.9.2 under ncarenv/23.09 |
| netCDF-C | netcdf-mpi/4.9.2 | 4.9.2 | |
| netCDF-Fortran | bundled in netcdf/4.9.2 | 4.6.1 | confirmed on derecho: `nf-config --version` with ncarenv/23.09 + netcdf-mpi/4.9.2 |
| PnetCDF | parallel-netcdf/1.12.3 | 1.12.3 | |
| ESMF | esmf/8.6.0-debug and esmf/8.6.0 | 8.6.0, both flavors | see "ESMF flavors" below |
| BLAS/LAPACK | cray-libsci | reference `lapack`/`blas` (dnf) | libsci is proprietary |
| PIO | parallelio/2.6.2 module | none installed | CIME builds PIO from CTSM's pinned ParallelIO submodule during the case build |
| conda, nco | loaded on derecho | not included | not needed for build-only CI; revisit when tests actually run |

`git` is also built from source (it is not part of derecho's stack): CIME's
cprnc and PIO builds clone `genf90` via git during the case build. git can't
be `dnf install`ed here because that pulls in `openssh`, whose setuid
`ssh-keysign` fails to `chown` under NCAR's rootless single-UID podman — see
[NEXT_STEPS.md](NEXT_STEPS.md) for the full story.

## ESMF flavors

Two ESMF 8.6.0 trees are installed:

- `/usr/local/esmf-8.6.0-debug` (`ESMF_BOPT=g`) — mirrors the
  `esmf/8.6.0-debug` module derecho loads for gnu `DEBUG=TRUE` builds.
  **This is the default** (`ESMFMKFILE` points here) because the current CI
  test is `SMS_D...` (debug).
- `/usr/local/esmf-8.6.0` (`ESMF_BOPT=O`) — mirrors `esmf/8.6.0`, for
  non-debug tests.

To use the optimized flavor in a workflow step:

```yaml
env:
  ESMFMKFILE: /usr/local/esmf-8.6.0/lib/esmf.mk
```

## Building the image

**On NCAR HPC (Casper), rootless podman** — use the wrapper, which sets the
node-local `TMPDIR` that rootless podman/buildah requires and carries PBS
headers for a build allocation:

```bash
docker/ctsm-ci-derecho-gnu/build-on-casper.sh
```

Request generous memory (e.g. `mem=256GB`): the build compiles the whole
toolchain and the final image commit is memory-hungry. See
[NEXT_STEPS.md](NEXT_STEPS.md) for the full set of rootless-podman build
constraints.

**On a host with unrestricted Docker/root** (a workstation or CI runner), a
plain build works:

```bash
podman build -t ctsm-ci-derecho-gnu:dev docker/ctsm-ci-derecho-gnu/   # or: docker build ...
```

Validate a built image with `docker/ctsm-ci-derecho-gnu/smoke-test.sh`.

## Bumping versions

Component versions are `ARG`s near the top of the `Dockerfile`
(`GCC_VERSION`, `MPICH_VERSION`, `HDF5_VERSION`, `NETCDF_C_VERSION`,
`NETCDF_FORTRAN_VERSION`, `PNETCDF_VERSION`, `ESMF_VERSION`), plus
`GIT_VERSION` and `MAKE_JOBS` (build parallelism). Change one and rebuild.

## Baked-in environment (why each matters)

| Variable | Value | Why |
|---|---|---|
| `PATH` / `LD_LIBRARY_PATH` | `/usr/local` first, then `/opt` mpich + gcc | non-login GHA shells need the toolchain on `PATH` without profile scripts |
| `CC/CXX/FC/F77` | MPICH wrappers (`mpicc`, …) | build everything against MPICH |
| `ESMFMKFILE` | `/usr/local/esmf-8.6.0-debug/lib/esmf.mk` | CMEPS/CDEPS builds require it; the `container` machine config does not set it |
| `CESMDATAROOT` | `/opt/cesmdata` | `DIN_LOC_ROOT=$CESMDATAROOT/inputdata` must exist even for `--no-run` |
| `USER` | `root` | CIME requires `$USER`; GHA container jobs don't set it |
| `PKG_CONFIG_PATH` / `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS` | `/usr/local/lib/pkgconfig` / `1` | CIME's cprnc locates netCDF via pkg-config with the non-system `/opt` toolchain |
