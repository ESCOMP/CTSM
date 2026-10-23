.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This documents some of the high-level design decisions made during implementation of
dynamic landunits.

===================================================================
 Switching between time interpolated or non-interpolated variables
===================================================================
The original intent was to make it easy to switch between using time interpolated or
non-interpolated variables for the variables read from the dynamic landcover file. This
change can be accomplished simply by changing a variable of type
``dyn_var_time_uninterp_type`` to ``dyn_var_time_interp_type``. However, there are
currently some barriers to making this change, as discussed below.

Sometime after implementing that simple switching mechanism, we added annual flux
dribblers. These dribblers turned out to be important when using non-time-interpolated
area changes, to interpolate some results of these changes even though the area changes
themselves were not interpolated. These dribblers currently enforce that changes happen
only on the year boundary (as an error check to ensure that we aren't getting changes at
unexpected times), so the current implementation now forces us to use
non-time-interpolated variables.

Therefore, changing to time interpolated variables would require either removing the
dribblers or removing the check that changes happen only on the year boundary. (Removing
the check can be done by changing the allows_non_annual_delta argument to the dribbler
constructor.) If all area changes were time interpolated then we might not need the
dribblers, but in practice some area changes - in particular glacier areas - will probably
continue to be non-time-interpolated, making it problematic to simply remove the
dribblers. The dribblers are set up to handle fluxes outside the year boundary by simply
applying the flux in the current time step. This means that it could probably work to keep
using the dribblers and just remove the check that prevents changes outside the year
boundary, though some analysis of the use of the dribblers would be warranted to ensure
that this will be okay for all uses of the dribblers. (There is then a small oddity in
that, for an interpolated area change, the change in the first time step of the year ends
up being further spread through the year by the dribbler. So the change in the first time
step of the year will be tiny, and there will instead be a slightly-larger-than-intended
change in every other time step of the year. But the error from this oddity should be very
small, so this is probably okay.)

Because of the use of dribblers with this check preventing non-year-boundary changes, we
also needed to add an explicit ``is_beg_curr_year`` check around calls to various interp
routines in ``dynSubgrid_driver``. This was needed to avoid this dribbler check from
getting triggered for some startup or hybrid runs that start mid-year. (See issue #4165
for more details.) If we used time interpolated variables, we would need to remove this
``is_beg_curr_year`` conditional; this would be possible in conjunction with having the
dribblers stop checking for non-year-boundary changes, as discussed above.

==================================================================================
 The use of dzsoi_decomp for urban landunits to calculate totcolch4 in ch4Mod.F90
==================================================================================
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
