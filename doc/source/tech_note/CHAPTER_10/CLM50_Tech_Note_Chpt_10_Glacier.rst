Glaciers
============

This chapter describes features of CLM that are specific to coupling to
an ice sheet model (in the CESM context, this is the Glimmer-CISM model;
Lipscomb and Sacks (2012) provide documentation and user’s guide for
Glimmer-CISM). General information about glacier land units can be found
elsewhere in this document (see Chapter 2 for an overview).

Overview
-------------

CLM is responsible for computing three quantities that are passed to the
ice sheet model:

#. Surface mass balance (SMB) – the net annual accumulation/ablation of
   mass at the upper surface (section 10.3)

#. Ground surface temperature, which serves as an upper boundary
   condition for Glimmer-CISM’s temperature calculation

#. Surface topography, which currently is fixed in time, and is provided
   on CLM’s surface dataset

The ice sheet model is typically run at much higher resolution than CLM
(e.g., :math:`\sim`\ 5 km rather than :math:`\sim`\ 100 km). To improve
the downscaling from CLM’s grid to the ice sheet grid, the glaciated
portion of each grid cell is divided into multiple elevation classes
(section 10.2). The above quantities are computed separately in each
elevation class. Glimmer-CISM then computes high-resolution quantities
via horizontal and vertical interpolation.

There are several reasons for computing the SMB in CLM rather than in
Glimmer-CISM:

#. It is much cheaper to compute the SMB in CLM for :math:`\sim`\ 10
   elevation classes than in Glimmer-CISM. For example, suppose we are
   running CLM at a resolution of :math:`\sim`\ 50 km and Glimmer at
   :math:`\sim`\ 5 km. Greenland has dimensions of about 1000 x 2000 km.
   For CLM we would have 20 x 40 x 10 = 8,000 columns, whereas for
   Glimmer we would have 200 x 400 = 80,000 columns.

#. We can use the sophisticated snow physics parameterization already in
   CLM instead of implementing a separate scheme for Glimmer-CISM. Any
   improvements to the CLM are applied to ice sheets automatically.

#. The atmosphere model can respond during runtime to ice-sheet surface
   changes. As shown by Pritchard et al. (2008), runtime albedo feedback
   from the ice sheet is critical for simulating ice-sheet retreat on
   paleoclimate time scales. Without this feedback the atmosphere warms
   much less, and the retreat is delayed.

#. Mass is more nearly conserved, given that the rate of surface ice
   growth or melting computed in CLM is equal to the rate seen by the
   dynamic ice sheet model. (Mass conservation is not exact, however,
   because of approximations made in interpolating from the CLM grid to
   the ice-sheet grid.)

#. The improved SMB is available in CLM for all glaciated grid cells
   (e.g., in the Alps, Rockies, Andes, and Himalayas), not just those
   which are part of ice sheets.

The current coupling between CLM and Glimmer-CISM is one-way only. That
is, CLM sends the SMB and surface temperature to Glimmer-CISM but does
not do anything with the fields that are returned. The CLM glacier
fraction and surface topography are therefore fixed in time. One-way
coupling is reasonable for runs of :math:`\sim`\ 100 years or less, in
which ice-sheet elevation changes are modest. For longer runs with
larger elevation changes, two-way coupling is highly desirable. A
two-way coupling scheme is under development.

Multiple elevation class scheme
------------------------------------

In the typical operation of CLM, the glacier land unit contains a single
column (section 2.1.1). However, when running CESM with an active ice
sheet model, the glacier land unit is replaced by a glacier\_mec land
unit, where “mec” denotes “multiple elevation classes”. In most ways,
glacier\_mec land units behave the same as standard glacier land units.
However, each glacier\_mec land unit is divided into a user-defined set
of columns based on surface elevation. The default is 10 elevation
classes whose lower limits are 0, 200, 400, 700, 1000, 1300, 1600, 2000,
2500, and 3000 m. Each column is characterized by a fractional area and
surface elevation that are read in during model initialization. Each
glacier\_mec column within a grid cell has distinct ice and snow
temperatures, snow water content, surface fluxes, and SMB.

Glacier\_mec columns, like glacier columns, are initialized with a
temperature of 250 K. While glacier columns are initialized with a snow
liquid water equivalent (LWE) equal to the maximum allowed value of 1 m,
glacier\_mec columns begin with a snow LWE of 0.5 m so that they will
reach their equilibrium mean snow depth sooner. Glacier\_mec columns
typically require several decades of spin-up to equilibrate with a given
climate.

The atmospheric surface temperature, potential temperature, specific
humidity, density, and pressure are downscaled from the mean grid cell
elevation to the glacier\_mec column elevation using a specified lapse
rate (typically 6.0 deg/km) and an assumption of uniform relative
humidity. At a given time, lower-elevation columns can undergo surface
melting while columns at higher elevations remain frozen. This gives a
more accurate simulation of summer melting, which is a highly nonlinear
function of air temperature. The precipitation rate and radiative fluxes
are not currently downscaled, but could be in the future if care were
taken to preserve the cell-integrated values.

In contrast to most CLM subgrid units, glacier\_mec columns can be
active (i.e., have model calculations run there) even if their area is
zero. This is done because the ice sheet model may require a SMB even
for some grid cells where CLM does not have glacier land units. To allow
for this, grid overlap files have been pre-computed. For given
resolutions of CLM and Glimmer-CISM, these files identify all
land-covered grid cells that overlap any part of the ice sheet grid. In
these overlapping cells, glacier\_mec columns are defined in all
elevation classes. Some columns may have zero area and are called
“virtual” columns. These columns do not affect energy exchange between
the land and the atmosphere, but are included for potential forcing of
Glimmer-CISM.

Computation of the surface mass balance
--------------------------------------------

The SMB of a glacier or ice sheet is the net annual
accumulation/ablation of mass at the upper surface. Ablation is defined
as the mass of water that runs off to the ocean. Not all the surface
meltwater runs off; some of the melt percolates into the snow and
refreezes. Accumulation is primarily by snowfall and deposition, and
ablation is primarily by melting and evaporation/sublimation. CLM uses a
surface-energy-balance (SEB) scheme to compute the SMB. In this scheme,
the melting depends on the sum of the radiative, turbulent, and
conductive fluxes reaching the surface, as described elsewhere in this
document.

CLM has a somewhat unrealistic treatment of accumulation and melting for
standard glacier land units. The snow depth is limited to a prescribed
depth of 1 m liquid water equivalent, with any additional snow assumed
to run off to the ocean. (This amounts to a crude parameterization of
iceberg calving.) Snow melting is treated in a realistic fashion, with
meltwater percolating downward through snow layers as long as the snow
is unsaturated. Once the underlying snow is saturated, any additional
meltwater runs off. When glacier ice melts, however, the meltwater is
assumed to remain in place until it refreezes. In warm parts of the ice
sheet, the meltwater does not refreeze, but stays in place indefinitely.

In the modified glacier\_mec columns, the treatment of melting and
freezing depends on the logical variable *glc\_dyntopo*. This variable
controls whether CLM surface topography changes dynamically as the ice
sheet evolves (i.e., whether the coupling is one-way or two-way). If
*glc\_dyntopo* is true, then CLM receives updated topographic
information from the ice sheet model. In this case, snow in excess of
the prescribed maximum depth is assumed to turn into ice, contributing a
positive SMB to the ice sheet model. Melting ice is assumed to run off
to the ocean, giving a negative SMB. The net SMB associated with ice
formation (by conversion from snow) and melting/runoff is computed for
each column, averaged over the coupling interval, and sent to the
coupler (*qice, mm/s*). If *glc\_dyntopo* is false, then surface runoff
for glacier\_mec land units is computed as for glacier land units: Any
snow in excess of 1 m LWE runs off to the ocean, and Melted ice remains
in place until it refreezes. Excess snow and melted ice still contribute
to positive and negative values, respectively, of *qice*, but only for
the purpose of forcing Glimmer-CISM. Currently, *glc\_dyntopo* = false
is the default, and the only supported option.

Note that the SMB typically is defined as the total accumulation of ice
and snow, minus the total ablation. The *qice* flux passed to
Glimmer-CISM is the mass balance for ice alone, not snow. We can think
of CLM as owning the snow, whereas Glimmer-CISM owns the underlying ice.
Fluctuations in snow depth between 0 and 1 m LWE are not reflected in
the SMB passed to Glimmer-CISM.
