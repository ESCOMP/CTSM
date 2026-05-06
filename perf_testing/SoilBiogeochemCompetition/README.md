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

- `SoilBiogeochemCompetition.F90` — the routine plus a set of
  per-element helper procedures (sibling subroutines/functions in the
  same module) factored out of the canonical-path science loops so the
  upcoming OpenACC `do`-loop apparatus doesn't pollute the science
  code. Helper layout:
  - `accum_sminn_tot`, `compute_nuptake_prof`, `accum_dz_weighted`,
    `compute_fraction_or_one`, `compute_residual_smin_vr`,
    `distribute_residual_to_plant` — per-element math.
  - `compete_nh4`, `compete_no3`, `compute_n2o_emissions`,
    `apply_carbon_only_adjustment`, `compute_competition_summary` —
    sub-blocks of the main competition loop (Loop 17).
  - The non-canonical branches (`use_nitrif_denitrif=.false.` block,
    Loop 19 MIMICS overflow) are intentionally not refactored.
  - Depends on intrinsic kinds (`selected_real_kind(12)` defines `r8`
    locally) and on [`../perf_timers_mod.F90`](../perf_timers_mod.F90)
    (which is a no-op when `INNER_TIMING` is undefined).
- `driver.F90` — synthetic timing harness; allocates inputs (pointer
  arrays where the routine signature requires pointer), runs all 8
  config combinations per iter (or 1 with `--fast`), prints results,
  writes `last_run.txt`, compares against `baseline_checksum.txt`.
- `baseline_checksum.txt` — committed reference output of the `--all`
  run (default params). Driver compares against this when run without
  `--fast` and the fingerprint matches.
- `baseline_checksum_fast.txt` — committed reference output of the
  canonical `--fast` run (default params). Driver compares against
  this when run with `--fast` and the fingerprint matches.
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
  (`use_nitrif_denitrif=.true.`, `carbon_only=.false.`, MIMICS off —
  i.e. `decomp_method /= mimics_decomp`). Use it for tight
  perf-iteration loops where covering every branch every time is
  unnecessary. Compares against its own baseline file
  (`baseline_checksum_fast.txt`), not the `--all` baseline.

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

### Per-loop ("inner") timing

A separate cpp macro `INNER_TIMING` enables per-loop wall-clock
instrumentation for the canonical-path science loops in
`SoilBiogeochemCompetition` (init / accum / nuptake-prof / main
competition / etc.). Default off:

```bash
make clean && make INNER_TIMING=1
./driver --fast
```

When enabled, each canonical loop's elapsed time and call count
are accumulated by [`../perf_timers_mod.F90`](../perf_timers_mod.F90)
(intrinsic `system_clock`, no external library). At the end of a
run the driver:

- prints a `--- per-loop wall-clock timers ---` table to stdout, and
- writes one row per label to `last_run_timings.csv` (gitignored;
  see [`../.gitignore`](../.gitignore)).

`INNER_TIMING` is independent of `TIMING` / `PERF_TIMING` (which
gates the driver-level total-time block), so you can measure
per-loop times alone, total time alone, both, or neither. Use
`./verify.sh INNER_TIMING=1` to build, run, and confirm both
`--fast` and `--all` still MATCH with timers on.

`make clean` removes `driver`, `*.o`, `*.mod`, and `last_run.txt`. It
does not touch `baseline_checksum.txt` or `baseline_checksum_fast.txt`.

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

Two committed baseline files, picked by mode:

- `baseline_checksum.txt` — `--all` mode (summed checksum across all 8
  configs).
- `baseline_checksum_fast.txt` — `--fast` mode (canonical config only).

Both serve as correctness references for future optimized variants.
The driver:

- prints `MATCH` if the parameter fingerprint matches the relevant
  baseline and the checksum agrees within `1e-10 * max(|baseline|, 1)`;
- prints `MISMATCH` (with the diff and tol) if the checksum has drifted
  — treat this as a correctness regression;
- skips the comparison when the fingerprint doesn't match (e.g.
  different sizes or `niters`).

To **regenerate** a baseline (e.g. after deliberately changing the
algorithm or input fill pattern):

```bash
# Regenerate the --all baseline
make clean && make && ./driver
cp last_run.txt baseline_checksum.txt
git add baseline_checksum.txt
git commit -m "Regenerate SoilBiogeochemCompetition baseline_checksum.txt"

# Regenerate the --fast baseline
./driver --fast
cp last_run.txt baseline_checksum_fast.txt
git add baseline_checksum_fast.txt
git commit -m "Regenerate SoilBiogeochemCompetition baseline_checksum_fast.txt"
```

## Notes for future optimization stages

- The canonical-path science loops have been factored into per-element
  helper procedures (see the helper layout under
  `SoilBiogeochemCompetition.F90` above). Each helper takes scalar
  args (no `(c,j)` indices); the surrounding `do j; do fc;` loops live
  in the main routine. This shape is OpenACC-friendly: each helper can
  be marked `!$acc routine seq` and the surrounding loop can be
  `!$acc parallel loop` without further restructuring.
- The first attempt at extracting Loop NH4's residual-uptake body wrapped
  the whole branchy loop body in a single `pure subroutine` and caused
  a +27% per-call regression at `-O2` because nvfortran wouldn't inline
  it. The current shape (extract just the pure-math expressions as
  `pure function`s, leave branches and accumulators at the call site)
  is what works — keep helpers small enough that the compiler inlines
  them.
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
