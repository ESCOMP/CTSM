# Standalone perf testing for `p2c_2d_filter`

`p2c_2d_filter` is a small subroutine in
[../../src/main/subgridAveMod.F90](../../src/main/subgridAveMod.F90) that
performs patch → column averaging on multi-level patch arrays. This
directory contains a self-contained copy plus a timing harness so the
routine can be built, run, and instrumented without the rest of the CTSM
build.

## Files

- `p2c_2d_filter.F90` — the subroutine. Self-contained: only depends on
  intrinsic kinds (`selected_real_kind(12)` defines `r8` locally).
- `driver.F90` — synthetic timing harness; allocates inputs, times
  `niters` calls, prints results, writes `last_run.txt`.
- `baseline_checksum.txt` — committed reference output of a canonical run
  (default params). The driver compares against this when present.
- `Makefile` — builds the driver. Default `FC=nvfortran`, override-able.
- `env.sh` — `. ./env.sh` to load `nvhpc` so `nvfortran` is on PATH.
- `clean.sh` — removes build artifacts and `last_run.txt`.

## Build & run

```bash
. ./env.sh                 # makes nvfortran available
make                       # builds ./driver with -O3 -g
./driver                   # runs with default params (8000 16 10 8000 1)
```

Override params positionally — `ncol npatch_per_col nlev numfc niters`:

```bash
./driver 16000 32 15 16000 100
```

Override compiler / flags at make time:

```bash
make clean && make FC=gfortran FFLAGS="-O2 -g -fopenmp"
make clean && make FFLAGS="-O3 -g -gpu=cc80 -acc"   # for OpenACC variants
```

### Disabling the built-in timing

The driver's internal `system_clock` instrumentation (and the `elapsed (s)`
/ `per call (s)` lines) is gated by the cpp macro `PERF_TIMING`. The
Makefile defines it by default. To leave wall-clock measurement to an
external profiler, build with:

```bash
make clean && make TIMING=0
```

The `TIMING=0` binary still allocates inputs, runs the call loop, computes
the checksum, and writes `last_run.txt` / compares against the baseline —
just nothing in the call loop's surrounding region except the loop itself.

## Output

Each run prints config + elapsed time + per-call time + checksum, e.g.:

```
=== p2c_2d_filter standalone driver ===
  ncol            = 8000
  npatch_per_col  = 16
  nlev            = 10
  numfc           = 8000
  niters          = 1
  npatch_total    = 128000
  elapsed (s)     =   4.526100E-03
  per call (s)    =   4.526100E-03
  checksum        =   5.1204800000000000E+09
  baseline        = MATCH (|diff| =   0.000000E+00)
```

Every run also writes `last_run.txt` (gitignored) with the same
parameter set + checksum.

## Baseline checksum

`baseline_checksum.txt` is committed. It captures the checksum of a
canonical run (default parameters) and serves as a correctness reference
for future optimized variants of `p2c_2d_filter`. The driver:

- prints `MATCH` if the param set matches the baseline and the checksum
  agrees within `1e-10 * max(|baseline|, 1)`;
- prints `MISMATCH` (with the diff and tol) if the checksum has drifted —
  treat this as a correctness regression;
- skips the comparison when the param set doesn't match the baseline,
  since the checksum is parameter-dependent.

To **regenerate** the baseline (e.g. after deliberately changing the
algorithm or input fill pattern):

```bash
./clean.sh && make && ./driver
cp last_run.txt baseline_checksum.txt
git add baseline_checksum.txt
git commit -m "Regenerate p2c_2d_filter baseline_checksum.txt"
```
