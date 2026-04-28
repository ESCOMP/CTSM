# Standalone perf testing for `SoilBiogeochemCompetition`

`SoilBiogeochemCompetition` is the routine in
[../../src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90](../../src/soilbiogeochem/SoilBiogeochemCompetitionMod.F90)
that resolves plant / heterotroph / nitrifier / denitrifier competition
for mineral N. This directory contains a self-contained copy plus a
timing harness so the routine can be built, run, and instrumented
without the rest of the CTSM build.

The standalone has **FUN dropped** (`local_use_fun = .false.`) — the
CNFUN call path, the `cnveg_*` / water / temperature / soilstate /
canopystate types, and the FUN-only consumer code are gone, so the
arg list is dramatically simpler than the in-tree routine. FUN-removal
is a single isolated commit on the branch and can be `git revert`ed if
FUN turns out to matter for perf.

Both `use_nitrif_denitrif` branches are preserved and runtime-selectable
via a driver arg.

## Files

- `SoilBiogeochemCompetition.F90` — the routine. Self-contained: only
  depends on intrinsic kinds (`selected_real_kind(12)` defines `r8`
  locally).
- `driver.F90` — synthetic timing harness; allocates inputs (pointer
  arrays where the routine signature requires pointer), times `niters`
  calls, prints results, writes `last_run.txt`, compares against
  `baseline_checksum.txt`.
- `baseline_checksum.txt` — committed reference output of the canonical
  run (default params). Driver compares against this when present.
- `Makefile` — tiny wrapper that sets `OBJ` and includes
  [../Makefile.common](../Makefile.common) (which carries `FC`,
  `FFLAGS`, the `PERF_TIMING` macro plumbing, and the `clean` target).

Shared across all `perf_testing/` subdirs (one level up):

- [../env.sh](../env.sh) — `. ../env.sh` to load `nvhpc` so `nvfortran`
  is on PATH.
- [../.gitignore](../.gitignore) — recursive; covers `*.o`, `*.mod`,
  `driver`, `last_run.txt`.
- [../Makefile.common](../Makefile.common) — shared build rules.

## Build & run

```bash
. ../env.sh                # makes nvfortran available (shared)
make                       # builds ./driver with -O3 -g -DPERF_TIMING
./driver                   # canonical params (8000 10 8 -1 1 .true. .false.)
```

Override params positionally — `ncol nlevdecomp ndct numfc niters use_nitrif_denitrif carbon_only`:

```bash
./driver 16000 15 12 16000 100 .true.  .false.   # bigger run, both branches
./driver 8000  10  8 -1    1   .false. .false.   # exercise the .not. nitrif branch
./driver 8000  10  8 -1    1   .true.  .true.    # exercise carbon_only path
```

`numfc=-1` is a sentinel meaning "use ncol".

Override compiler / flags at make time:

```bash
make clean && make FC=gfortran FFLAGS="-O2 -g -fopenmp"
make clean && make FFLAGS="-O3 -g -gpu=cc80 -acc"   # for OpenACC variants
```

### Disabling the built-in timing

The driver's internal `system_clock` instrumentation (and the
`elapsed (s)` / `per call (s)` lines) is gated by the cpp macro
`PERF_TIMING`. The Makefile defines it by default. To leave wall-clock
measurement to an external profiler:

```bash
make clean && make TIMING=0
```

The `TIMING=0` binary still allocates inputs, runs the call loop,
computes the checksum, and writes `last_run.txt` / compares against the
baseline — just nothing in the call loop's surrounding region except
the loop itself.

`make clean` removes `driver`, `*.o`, `*.mod`, and `last_run.txt`. It
does not touch `baseline_checksum.txt`.

## Output

Each run prints config + (if `PERF_TIMING`) elapsed/per-call time +
checksum, e.g.:

```
=== SoilBiogeochemCompetition standalone driver ===
  ncol                 = 8000
  nlevdecomp           = 10
  ndct                 = 8
  numfc                = 8000
  niters               = 1
  use_nitrif_denitrif  = T
  carbon_only          = F
  decomp_method        = 2
  elapsed (s)          =   8.399000E-03
  per call (s)         =   8.399000E-03
  checksum             =   9.5970435393765438E+04
  baseline             = MATCH (|diff| =   0.000000E+00)
```

Every run also writes `last_run.txt` (gitignored) with the parameter
fingerprint + checksum.

## Baseline checksum

`baseline_checksum.txt` is committed. It captures the checksum of the
canonical run (default parameters) and serves as a correctness reference
for future optimized variants of `SoilBiogeochemCompetition`. The
driver:

- prints `MATCH` if the parameter fingerprint matches the baseline and
  the checksum agrees within `1e-10 * max(|baseline|, 1)`;
- prints `MISMATCH` (with the diff and tol) if the checksum has drifted
  — treat this as a correctness regression;
- skips the comparison when the parameter set doesn't match the
  baseline, since the checksum is parameter-dependent (e.g. the
  `use_nitrif_denitrif` flag changes which code paths execute).

To **regenerate** the baseline (e.g. after deliberately changing the
algorithm or input fill pattern):

```bash
make clean && make && ./driver
cp last_run.txt baseline_checksum.txt
git add baseline_checksum.txt
git commit -m "Regenerate SoilBiogeochemCompetition baseline_checksum.txt"
```

## Notes for future optimization stages

- The routine accumulates into `actual_immob`, `potential_immob`,
  `sminn_to_plant` in the non-nitrif branch (uses `+=` without
  re-zeroing). The driver zeros them once before the call loop;
  setting `niters > 1` makes the checksum reflect cumulative state and
  the canonical baseline (niters=1) won't match.
- Pointer attributes on the pointer args mirror the in-tree
  `soilbiogeochem_*_inst%*_col` declarations. Don't switch them to
  assumed-shape `intent(in/out)` during extraction; save signature
  changes for explicit optimization commits.
- `cascade_receiver_pool`, `landunit` are integer pointers in CTSM and
  are passed as such here. `dzsoi_decomp` is `allocatable` in CTSM
  (declared as assumed-shape `intent(in)` here, which accepts
  allocatable / pointer / regular contiguous arrays).
