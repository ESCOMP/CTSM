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

- **`Ld5` gets `00:20:00`, and that holds everywhere** — single point, `f10_f10_mg37`,
  `f19_g17`, `f09_g17` and `ne30pg3_t232` all use it today. Short tests are dominated by
  startup, not by gridcells.
- **Long durations are where resolution starts to bite, and a single-point figure must not
  be carried across.** At `f10_f10_mg37`, roughly-one-year runs (`Ld396`) spread from
  `00:40:00` to `02:00:00` with the test type, and two- to three-year runs (`Ld765`,
  `Ld1096`, `Ly3`) cluster at `01:30:00`–`01:40:00`. At `f19_g17`, thirty-odd times the
  gridcells, the one existing `Ly2` carries `02:00:00`. Single-point multi-year runs range
  from `0:50:00` for about two years up to `04:00:00` for twenty.

  The spread within one grid and duration is the point: test type matters nearly as much as
  size, so match both before rounding up.
- **Exact-restart runs the model about 1.5 times over**, so scale the matched value
  accordingly — a single-point `ERS_Ld731` gets `01:00:00`.

These anchors are the state of the test list, not a specification. Re-read the file rather
than trusting this page if a number matters.

## Format

Write `HH:MM:SS`. Existing entries include oddities — `0:50:00` without the leading zero,
and `00:60:00` where `01:00:00` was meant — so copy the *value* from a comparable entry,
not its spelling.
