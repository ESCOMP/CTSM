---
name: build-ctsm-ci-image
description: Use when building, rebuilding, or revalidating the ctsm-ci-derecho-gnu container image on Casper — including after any change to docker/ctsm-ci-derecho-gnu/Dockerfile, and when deciding between an interactive allocation and a batch job.
---

# Building and validating ctsm-ci-derecho-gnu

All commands run **from the repo root**. The image is built on Casper with rootless podman; `build-on-casper.sh` exists to handle the node-local `TMPDIR` that buildah requires, so always go through it rather than calling `podman build` yourself.

## Pick the path first

The build and the validation runs need the image in podman's storage, and that storage is **node-local** — it dies with the allocation. So the two halves either share one allocation or the image has to travel between them on disk.

| | Interactive | Batch + load |
|---|---|---|
| Build and validate in one allocation | yes | no |
| Ties up a session for ~1h | yes | no |
| Costs a queue wait and a ~3.3 GB write | no | yes |

### Interactive (fewer moving parts)

```
execcasper -A <account> -l select=1:ncpus=16:mem=256GB -l walltime=06:00:00
docker/ctsm-ci-derecho-gnu/build-on-casper.sh
```

Build is roughly 50 minutes, so 6h leaves room for validation afterward. Then run the validation section below in the same shell.

### Batch, then load in a second allocation

```
qsub -A <account> docker/ctsm-ci-derecho-gnu/build-on-casper.sh
```

**A `qsub`ed build cannot also validate** — the job ends when the script does. On success the script saves the image to GLADE and prints `Restore with: podman load -i ...`. In a later allocation:

```
podman load -i /glade/work/$USER/ctsm-ci-derecho-gnu_YYYYMMDD.tar
```

The tar is named for the build date, but a second build on the same day gets a timestamped name instead so it never clobbers a known-good tarball — read the script's output rather than assuming the plain date form. After loading, the image is `localhost/ctsm-ci-derecho-gnu:dev` again.

## Validate

Every wrapper defaults to `localhost/ctsm-ci-derecho-gnu:dev`, which is what the build produces, so no `IMAGE_TAG` is needed. Run all five:

```
docker/ctsm-ci-derecho-gnu/smoke-test.sh
docker/ctsm-ci-derecho-gnu/run-unit-tests-in-container.sh
docker/ctsm-ci-derecho-gnu/run-test-in-container.sh
docker/ctsm-ci-derecho-gnu/run-sys-tests-in-container.sh -t SMS_D_Ld1_Mmpi-serial.1x1_brazil.IHistClm60BgcQianRsGs
docker/ctsm-ci-derecho-gnu/run-case-in-container.sh --case brazil_test --compset IHistClm60BgcQianRsGs --res 1x1_brazil --run-unsupported
```

`--machine container`, `--compiler gnu` and `--mpilib mpi-serial` are injected for you; anything you pass yourself wins.

## The acceptance test is on the executable

A build that succeeds proves nothing on its own — a case build that quietly compiles CTSM's own submodules passes every version check while linking something different from derecho. That is the failure this image exists to end, so check the binary. Against the `cesm.exe` from `run-test-in-container.sh`:

```
strings -a <bld>/cesm.exe | grep -c /ctsm/libraries/mpi-serial   # must be 0
strings -a <bld>/cesm.exe | grep -c /ctsm/libraries/parallelio   # must be 0
strings -a <bld>/cesm.exe | grep -c /usr/local                   # must be nonzero
```

The first two being zero is the real result; the third is the positive control that stops a zero from meaning "grepped the wrong file".

## Don't

- Don't call `podman build` directly — the node-local `TMPDIR` handling lives in the wrapper, and without it the build dies creating the rootfs before any Dockerfile instruction runs.
- Don't set `CTSM_BUILD_NO_SAVE=1` on a real build. Podman's storage dies with the allocation and hours of compute go with it.
- Don't trust a version check that passes without the `strings` test above.
