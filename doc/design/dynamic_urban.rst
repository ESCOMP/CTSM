.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This documents some of the high-level design decisions made during implementation of
dynamic urban landunits.

============================================================================
 The use of dzsoi_decomp for urban landunits to calculate totcolch4 in ch4Mod.F90
============================================================================
During the first test simulation for dynamic urban, we encountered a methane conservation
error the first time PCT_URBAN changed. The dynamic adjustments for conc_ch4_sat_col and
conc_ch4_unsat_col (the column_state_updater in subroutine DynamicColumnAdjustments within
ch4Mod.F90) were distributing non-zero values for roof and walls for layers 1,nlevsoi.
When the total column ch4 is summed over the soil layers (or in this case, urban layers), the
summation is done over nlevsoi, not nlevurb, using dz. dz is 1.e36 for roof/wall layers
that are greater than nlevurb, thus creating an imbalance.

Rather than trying to keep the BGC variables physically meaningful in urban landunits,
we will just pack these variables in a way that should conserve these variables, even if
the values in each of the urban columns is somewhat nonsensical. Specifically: we'll take
col%wtgcell at face value in urban columns in dynColumnStateUpdaterMod - i.e., for the sake
of storing / conserving these BGC variables, we'll act as if that gives the true column
weight on the grid cell. This way we'll end up storing all of the C & N from the vegetated
column in the urban columns, and there shouldn't be any that is lost from the system. If that
urban landunit later shrinks, the stored C & N should be restored symmetrically. It shouldn't
really matter that it was stored in a non-physical way (e.g., with some C & N stored in urban
walls), since the BGC variables are irrelevant over the urban areas and we just want to be able
to restore the amount that was originally stored if an urban landunit grows and then later shrinks.
But for this to work right, we need to treat the relevant BGC variables as having the same dz over
all urban columns as over the soil column. Note that there already seems to be an implicit assumption
that dz is the same for all columns in the dynamic column state updates, in that dz doesn't enter
into the conservation equations. In terms of what needs to change, we think that the only relevant
code is the code that sums up total C / N / CH4 for the sake of balance checks: these balance checks
need to be consistent with the assumptions made in the conservation code. The C and N summations
already use dzsoi_decomp, which is the same for all columns, so this is already what we want.
The only thing that needs to change is the use of dz in totcolch4 in ch4Mod.F90: we've changed that to now use
dzsoi_decomp over urban columns. (This begs the question of why this isn't already using
dzsoi_decomp for consistency with the C & N code; we're not sure about this.)

See issue #1445 for the original discussion on this topic.
