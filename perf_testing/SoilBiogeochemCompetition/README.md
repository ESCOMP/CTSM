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

Both `use_nitrif_denitrif` branches are preserved, plus the
`carbon_only` and `decomp_method == mimics_decomp` switches. The driver
exercises all 8 combinations by default (see [Driver modes](#driver-modes)).

## Files

- `SoilBiogeochemCompetition.F90` — the routine. Self-contained: only
  depends on intrinsic kinds (`selected_real_kind(12)` defines `r8`
  locally).
- `driver.F90` — synthetic timing harness; allocates inputs (pointer
  arrays where the routine signature requires pointer), runs all 8
  config combinations per iter (or 1 with `--fast`), prints results,
  writes `last_run.txt`, compares against `baseline_checksum.txt`.
- `baseline_checksum.txt` — committed reference output of the canonical
  run (default params, `--all` mode). Driver compares against this when
  the fingerprint matches.
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
make                       # builds ./driver with -DPERF_TIMING
./driver                   # default: --all mode, 8 configs, default sizes
./driver --fast            # canonical config only (1 call instead of 8)
```

Override sizes positionally — `ncol nlevdecomp ndct numfc niters`:

```bash
./driver 16000 15 12 16000 100         # bigger --all run
./driver --fast 16000 15 12 16000 100  # bigger run, single config
```

`numfc=-1` (default) means "use ncol". `--fast` may appear before or
after the positional args.

## Driver modes

- **`--all` (default)** — runs each iter as 8 calls covering every
  combination of (`use_nitrif_denitrif`, `carbon_only`,
  `decomp_method == mimics_decomp`). The reported checksum is the sum
  of all 8 per-config checksums, so it locks correctness across every
  top-level branch in the routine. Per-call time = elapsed / (niters * 8).
- **`--fast`** — runs only the canonical config
  (`use_nitrif_denitrif=.true.`, `carbon_only=.false.`,
  `decomp_method=mimics_decomp`). Use it for tight perf-iteration loops
  where covering every branch every time is unnecessary. Has its own
  fingerprint, so an `--all` baseline doesn't `MATCH` a `--fast` run
  (the driver prints `param set differs; skipping compare`).

Within a single config, the synthetic inputs are rigged so per-cell
branches inside the routine fire on different cells (`sminn_vr` ranges
0.05..2 g/m3 so both `supply > demand` and `supply < demand` cells
exist; `cascade_receiver_pool` cycles 1..5 so the MIMICS pool-id test
both hits and misses).

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
checksum, e.g. (default `--all` mode):

```
=== SoilBiogeochemCompetition standalone driver ===
  mode                 = all
  ncol                 = 8000
  nlevdecomp           = 10
  ndct                 = 8
  numfc                = 8000
  niters               = 1
  configs / iter       = 8
  total calls          = 8
  elapsed (s)          =   4.957000E-02
  per call (s)         =   6.196250E-03
  checksum             =   7.6772246368780360E+05
  baseline             = MATCH (|diff| =   0.000000E+00)
```

Every run also writes `last_run.txt` (gitignored) with the parameter
fingerprint (`mode`, sizes, `niters`) + checksum.

## Baseline checksum

`baseline_checksum.txt` is committed. It captures the summed checksum
of the canonical default run (`--all` mode, default sizes) and serves
as a correctness reference for future optimized variants. The driver:

- prints `MATCH` if the parameter fingerprint matches the baseline and
  the checksum agrees within `1e-10 * max(|baseline|, 1)`;
- prints `MISMATCH` (with the diff and tol) if the checksum has drifted
  — treat this as a correctness regression;
- skips the comparison when the fingerprint doesn't match (e.g.
  different sizes, `niters`, or running `--fast` against an `--all`
  baseline).

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
  `sminn_to_plant` in the non-nitrif branch (`+=` without re-zeroing).
  The driver re-zeros all output / inout arrays before every call, so
  per-call results are independent of `niters` and config order.
- Pointer attributes on the pointer args mirror the in-tree
  `soilbiogeochem_*_inst%*_col` declarations. Don't switch them to
  assumed-shape `intent(in/out)` during extraction; save signature
  changes for explicit optimization commits.
- `cascade_receiver_pool`, `landunit` are integer pointers in CTSM and
  are passed as such here. `dzsoi_decomp` is `allocatable` in CTSM
  (declared as assumed-shape `intent(in)` here, which accepts
  allocatable / pointer / regular contiguous arrays).
