# The `water_type` factory

`unittestWaterTypeFactory` in `src/unit_test_shr/` builds a complete `water_type` so you
do not have to hand-roll one. Its calls have a required order, and one of them does more
than its name suggests. Worked example in the tree: `test_irrigation.pf`.

## The order

```fortran
call this%water_factory%init()
call this%water_factory%setup_before_subgrid( &
     my_nlevsoi = 3, nlevgrnd_additional = 1, my_nlevsno = 3)
! ... build the subgrid here ...
call this%water_factory%setup_after_subgrid(dz = ...)     ! snl = ... is optional
call this%water_factory%create_water_type(this%water_inst)
! ... and in tearDown:
call this%water_factory%teardown(this%water_inst)
```

Those five bindings are the whole public interface.

## `setup_before_subgrid` sets the level counts itself

It sets `nlevsoi`, `nlevgrnd`, `nlevmaxurbgrnd` and `nlevsno` — so do not also set
`nlevsno` directly and expect your value to survive. If you have already set them, pass
the same values back in to make the call a no-op, which is what its own comment tells you
to do.

## `snl` is one value for every column

`setup_after_subgrid`'s optional `snl` sets every column to the same snow-layer count. If
you need columns that differ, set `col%snl` per column instead.
