# Choosing a wallclock

**Round up, and scale by resolution — not by whether the test is "global".**

The cost is asymmetric. A generous wallclock costs some queue priority; an underestimate
kills the test after it has already burned the time, and you pay for the whole run twice.

The driver is the land-gridcell count, and the grids in `testlist_clm.xml` span orders of
magnitude: single point, then `f10_f10_mg37` at 10°x15° (the coarsest global grid),
`f45_f45_mg37` at 4°x5°, `f19_g17` at ~1.9°x2.5°, `f09_g17` at ~0.9°x1.25°, and
`ne30pg3_t232` at ~1°. "Global" is not the axis — `f10_f10` is cheap and `f09_g17` is not.

## Derive the value, do not guess it

Copy the wallclock from an existing entry **at the same grid and a comparable duration**,
and round up when there is no close match. Say in the change where the figure came from, so
a reviewer can check the comparison rather than squinting at the number.

To find comparable entries:

```
./cime/scripts/query_testlists --show-options
```

`--show-options` is what prints each test's wallclock. Narrow the listing only if it is
unwieldy, and narrow on something relevant to the comparison you are making —
`--xml-machine`, `--xml-compiler`, `--xml-category`. Reading `testlist_clm.xml` directly
works too, and is easier when what you want is every entry at one grid, which is usually
exactly what you want.

## Anchors

- **`Ld5` is `00:20:00` at nearly every grid** — single point, `f10_f10_mg37`, `f19_g17`,
  `f09_g17` and `ne30pg3_t232` all use it today, because a short test is dominated by
  startup rather than by gridcells. **It is not uniform, though, and the exceptions are not
  the large grids:** a dozen-odd entries run `00:30:00` to `01:00:00`, most of them `ERP_D`
  or `ERS_D` at `f10_f10_mg37`, which is the *cheapest* global grid. Match the test type
  before taking the `00:20:00`.
- **Long durations are where resolution starts to bite, and a single-point figure must not
  be carried across.** At `f10_f10_mg37`, roughly-one-year runs (`Ld396`) spread from
  `00:40:00` to `02:00:00` with the test type, and two- to three-year runs (`Ld765`,
  `Ld1096`, `Ly3`) mostly sit at `01:30:00`–`01:40:00`, with one `ERS_Ld765` at `00:60:00`
  below that and the `Ly2` feature tests at `02:00:00` above it. At `f19_g17`, thirty-odd
  times the gridcells, the one existing `Ly2` carries `02:00:00`.

  The spread within one grid and duration is the point: test type matters nearly as much as
  size, so match both before rounding up.
- **Single-point figures do not interpolate by duration at all.** The list holds a
  single-point `Ly5` at `00:20:00` and another at `01:00:00`, while a single-point run of
  about two years carries `0:50:00` and one of twenty carries `04:00:00`. What the compset
  does dominates, so at a single point find an entry with the same test type and compset
  and take its value; do not scale one duration to another.
- **Exact-restart runs the model about 1.5 times over**, so scale the matched value
  accordingly — a single-point `ERS_Ld731` gets `01:00:00`.

These anchors are the state of the test list, not a specification. Re-read the file rather
than trusting this page if a number matters.

## Format

Write `HH:MM:SS`. Existing entries include oddities — `0:50:00` without the leading zero,
and `00:60:00` where `01:00:00` was meant — so copy the *value* from a comparable entry,
not its spelling.
