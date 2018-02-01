.. _rst_Glaciers:

Glaciers
============

This chapter describes features of CLM that are specific to coupling to
an ice sheet model (in the CESM context, this is the CISM model;
:ref:`Lipscomb and Sacks (2012)<LipscombSacks2012>` provide
documentation and user’s guide for CISM). General information
about glacier land units can be found elsewhere in this document (see
Chapter :numref:`rst_Surface Characterization, Vertical Discretization,
and Model Input Requirements` for an overview).

.. _Overview Glaciers:

Overview
--------

CLM is responsible for computing two quantities that are passed to the
ice sheet model:

#. Surface mass balance (SMB) - the net annual accumulation/ablation of
   mass at the upper surface (section 
   :numref:`Computation of the surface mass balance`)

#. Ground surface temperature, which serves as an upper boundary
   condition for CISM's temperature calculation

The ice sheet model is typically run at much higher resolution than CLM
(e.g., :math:`\sim`\ 5 km rather than :math:`\sim`\ 100 km). To improve
the downscaling from CLM’s grid to the ice sheet grid, the glaciated
portion of each grid cell is divided into multiple elevation classes
(section :numref:`Multiple elevation class scheme`). The above
quantities are computed separately in each elevation class. The CESM
coupler then computes high-resolution quantities via horizontal and
vertical interpolation, and passes these high-resolution quantities to
CISM.

There are several reasons for computing the SMB in CLM rather than in
CISM:

#. It is much cheaper to compute the SMB in CLM for :math:`\sim`\ 10
   elevation classes than in CISM. For example, suppose we are
   running CLM at a resolution of :math:`\sim`\ 50 km and CISM at
   :math:`\sim`\ 5 km. Greenland has dimensions of about 1000 x 2000 km.
   For CLM we would have 20 x 40 x 10 = 8,000 columns, whereas for
   CISM we would have 200 x 400 = 80,000 columns.

#. We can use the sophisticated snow physics parameterization already in
   CLM instead of implementing a separate scheme for CISM. Any
   improvements to the CLM are applied to ice sheets automatically.

#. The atmosphere model can respond during runtime to ice-sheet surface
   changes (even in the absence of two-way feedbacks with CISM). As
   shown by :ref:`Pritchard et al. (2008)<Pritchardetal2008>`, runtime
   albedo feedback from the ice sheet is critical for simulating
   ice-sheet retreat on paleoclimate time scales. Without this feedback
   the atmosphere warms much less, and the retreat is delayed.

#. The improved SMB is available in CLM for all glaciated grid cells
   (e.g., in the Alps, Rockies, Andes, and Himalayas), not just those
   which are part of ice sheets.

In typical runs, CLM computes the SMB and sends it to CISM, but CISM's
ice sheet geometry remains fixed over the course of the run. In these
runs, CISM serves two roles in the system:

#. Over the CISM domain (typically Greenland in CESM2), CISM dictates
   glacier areas and topographic elevations, overriding the values on
   CLM's surface dataset. CISM also dictates the elevation of
   non-glacier land units in its domain, and atmospheric fields are
   downscaled to non-glacier land units only in this domain. (So if you
   run with a stub glacier model - SGLC - then glacier areas and
   elevations will be taken entirely from CLM's surface dataset, and no
   downscaling will be done over non-glacier land units.)

#. CISM provides the grid onto which SMB is downscaled. (If you run with
   SGLC then SMB will still be computed in CLM, but it won't be
   downscaled to a high-resolution ice sheet grid.)

However, it is also possible to run CESM with an evolving ice sheet. In
this case, CLM responds to CISM's evolution by adjusting the areas of
the glacier landunit and each elevation class within this landunit, as
well as the mean topographic heights of each elevation class. Thus,
CLM's glacier areas and elevations remain in sync with
CISM's. Conservation of mass and energy is done as for other landcover
change (see Chapter :numref:`rst_Transient Landcover Change`).

.. _Glacier regions:

Glacier regions and their behaviors
-----------------------------------

The world's glaciers and ice sheets are broken down into a number of
different regions (four by default) that differ in three respects:

#. Whether the gridcell's glacier landunit contains:

   a. Multiple elevation classes (section :numref:`Multiple elevation
      class scheme`)

   b. Multiple elevation classes plus virtual elevation classes

   c. Just a single elevation class whose elevation matches the
      atmosphere's topographic height (so there is no adjustment in
      atmospheric forcings due to downscaling).

#. Treatment of glacial melt water:

   a. Glacial melt water runs off and is replaced by ice, thus keeping
      the column always frozen. This behavior is discussed in more
      detail in section :numref:`Computation of the surface mass
      balance`.

   b. Glacial melt water remains in place until it refreezes - possibly
      remaining in place indefinitely if the glacier column is in a warm
      climate. With this behavior, ice melt does not result in any
      runoff. Regions with this behavior cannot compute SMB, because
      negative SMB would be meaningless (due to the liquid water on top
      of the ice column). This behavior produces less realistic glacier
      physics. However, it avoids the negative ice runoff that is needed
      for the "replaced by ice" behavior to conserve mass and energy (as
      described in section :numref:`Computation of the surface mass
      balance`). Thus, in regions where CLM has glaciers but the
      atmospheric forcings are too warm to sustain those glaciers, this
      behavior avoids persistent negative ice runoff. This situation can
      often occur for mountain glaciers, where topographic smoothing in
      the atmosphere results in a too-warm climate. There, avoiding
      persistent negative ice runoff can be more important than getting
      the right glacier ice physics.

#. Treatment of ice runoff from snow capping (as described in section
   :numref:`Runoff from glaciers and snow-capped surfaces`). Note that this
   is irrelevant in regions with an evolving, two-way-coupled ice sheet
   (where the snow capping term is sent to CISM rather than running off):

   a. Ice runoff from snow capping remains ice. This is a crude
      parameterization of iceberg calving, and so is appropriate in
      regions where there is substantial iceberg calving in reality.

   b. Ice runoff from snow capping is melted (generating a negative
      sensible heat flux) and runs off as liquid. This matches the
      behavior for non-glacier columns. This is appropriate in regions
      that have little iceberg calving in reality. This can be important
      to avoid unrealistic cooling of the ocean and consequent runaway
      sea ice growth.

The default behaviors for the world's glacier and ice sheet regions are
described in :numref:`Table Glacier region behaviors`. Note that the
standard CISM grid covers Greenland plus enough surrounding area to
allow for ice sheet growth and to have a regular rectangular grid. We
need to have the "replaced by ice" melt behavior within the CISM domain
in order to compute SMB there, and we need virtual elevation classes in
that domain in order to compute SMB for all elevation classes and to
facilitate glacial advance and retreat in the two-way-coupled
case. However, this domain is split into Greenland itself and areas
outside Greenland so that ice runoff in the Canadian archipelago (which
is inside the CISM domain) is melted before reaching the ocean, to avoid
runaway sea ice growth in that region.

.. _Table Glacier region behaviors:

.. table:: Glacier region behaviors

 +---------------+---------------+---------------+---------------+
 | Region        | Elevation     | Glacial melt  | Ice runoff    |
 |               | classes       |               |               |
 +===============+===============+===============+===============+
 | Greenland     | Virtual       | Replaced by   | Remains ice   |
 |               |               | ice           |               |
 +---------------+---------------+---------------+---------------+
 | Inside        | Virtual       | Replaced by   | Melted        |
 | standard CISM |               | ice           |               |
 | grid but      |               |               |               |
 | outside       |               |               |               |
 | Greenland     |               |               |               |
 | itself        |               |               |               |
 +---------------+---------------+---------------+---------------+
 | Antarctica    | Multiple      | Replaced by   | Remains ice   |
 |               |               | ice           |               |
 +---------------+---------------+---------------+---------------+
 | All others    | Single        | Remains in    | Melted        |
 |               |               | place         |               |
 +---------------+---------------+---------------+---------------+


.. _Multiple elevation class scheme:

Multiple elevation class scheme
-------------------------------

The glacier landunit contains multiple columns based on surface
elevation. These are known as elevation classes, and the land unit is
referred to as glacier\_mec. (As described in section :numref:`Glacier
regions`, some regions have only a single elevation class, but they are
still referred to as glacier\_mec landunits.) The default is to have 10
elevation classes whose lower limits are 0, 200, 400, 700, 1000, 1300,
1600, 2000, 2500, and 3000 m. Each column is characterized by a
fractional area and surface elevation that are read in during model
initialization, and then possibly overridden by CISM as the run
progresses. Each glacier\_mec column within a grid cell has distinct ice
and snow temperatures, snow water content, surface fluxes, and SMB.

The atmospheric surface temperature, potential temperature, specific
humidity, density, and pressure are downscaled from the atmosphere's
mean grid cell elevation to the glacier\_mec column elevation using a
specified lapse rate (typically 6.0 deg/km) and an assumption of uniform
relative humidity. Longwave radiation is downscaled by assuming a linear
decrease in downwelling longwave radiation with increasing elevation
(0.032 W m :sup:`-2` m :sup:`-1`, limited to 0.5 - 1.5 times the
gridcell mean value, then normalized to conserve gridcell total energy)
:ref:`(Van Tricht et al., 2016)<VanTrichtetal2016>`. Total precipitation
is partitioned into rain vs. snow as described in Chapter
:numref:`rst_Surface Characterization, Vertical Discretization, and
Model Input Requirements`. The partitioning of precipitation is based on
the downscaled temperature, allowing rain to fall at lower elevations
while snow falls at higher elevations.

This downscaling allows lower-elevation columns to undergo surface
melting while columns at higher elevations remain frozen. This gives a
more accurate simulation of summer melting, which is a highly nonlinear
function of air temperature.

Within the CISM domain, this same downscaling procedure is also applied
to all non-urban land units. The elevation of non-glacier land units is
taken from the mean elevation of ice-free grid cells in CISM. This is
done in order to keep the glaciated and non-glaciated portions of the
CISM domain as consistent as possible.

In contrast to most CLM subgrid units, glacier\_mec columns can be
active (i.e., have model calculations run there) even if their area is
zero. These are known as "virtual" columns. This is done because the ice
sheet model may require a SMB even for some grid cells where CLM does
not have glacier land units. Virtual columns do not affect energy
exchange between the land and the atmosphere, but are included for
potential forcing of CISM.

.. _Computation of the surface mass balance:

Computation of the surface mass balance
---------------------------------------

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
the purpose of forcing CISM. Currently, *glc\_dyntopo* = false
is the default, and the only supported option.

Note that the SMB typically is defined as the total accumulation of ice
and snow, minus the total ablation. The *qice* flux passed to CISM is
the mass balance for ice alone, not snow. We can think of CLM as owning
the snow, whereas CISM owns the underlying ice.  Fluctuations in snow
depth between 0 and 10 m water equivalent are not reflected in the SMB
passed to CISM. In transient runs, this can lead to delays of a few
decades in the onset of accumulation or ablation in a given glacier
column.

In regions where SMB is computed for glaciers, SMB is also computed for
the natural vegetated land unit. Because there is no ice to melt in this
land unit, it can only generate a zero or positive SMB. A positive SMB
is generated once the snow pack reaches its maximum depth. When running
with an evolving ice sheet, this condition triggers glacial inception.

FIXME: Make sure I talk about the positive liquid runoff and negative
ice runoff that result from melted ice.

FIXME: With interactive CISM, snow capping is applied to the surface
mass balance, not ice runoff.
