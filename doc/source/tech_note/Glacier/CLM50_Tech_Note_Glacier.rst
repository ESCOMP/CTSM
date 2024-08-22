.. _rst_Glaciers:

Glaciers
========

This chapter describes features of CLM that are specific to coupling to an ice sheet model (in the CESM context, this is the CISM model; :ref:`Lipscomb and Sacks (2012)<LipscombSacks2012>` provide documentation and user's guide for CISM). General information about glacier land units can be found elsewhere in this document (see Chapter :numref:`rst_Surface Characterization, Vertical Discretization, and Model Input Requirements` for an overview).

.. _Glaciers summary of CLM5.0 updates relative to CLM4.5:

Summary of CLM5.0 updates relative to CLM4.5
--------------------------------------------

Compared with CLM4.5 (:ref:`Oleson et al. 2013 <Olesonetal2013>`), CLM5.0 contains substantial improvements in its capabilities for land-ice science. This section summarizes these improvements, and the following sections provide more details.

- All runs include multiple glacier elevation classes over Greenland and Antarctica and compute ice sheet surface mass balance in those regions.

- A number of namelist parameters offer fine-grained control over glacier behavior in different regions of the world (section :numref:`Glacier regions`). (The options used outside of Greenland and Antarctica reproduce the standard CLM4.5 glacier behavior.)

- CLM can now keep its glacier areas and elevations in sync with CISM when running with an evolving ice sheet. (However, in typical configurations, the ice sheet geometry still remains fixed throughout the run.)

- The downscaling to elevation classes now includes downwelling longwave radiation and partitioning of precipitation into rain vs. snow (section :numref:`Multiple elevation class scheme`).

- Other land units within the CISM domain undergo the same downscaling as the glacier land unit, and surface mass balance is computed for the natural vegetated land unit. This allows CLM to produce glacial inception when running with an evolving ice sheet model.

- There have also been substantial improvements to CLM's snow physics, as described in other chapters of this document.

.. _Overview Glaciers:

Overview
--------

CLM is responsible for computing two quantities that are passed to the ice sheet model:

#. Surface mass balance (SMB) - the net annual accumulation/ablation of mass at the upper surface (section :numref:`Computation of the surface mass balance`)

#. Ground surface temperature, which serves as an upper boundary condition for CISM's temperature calculation The ice sheet model is typically run at much higher resolution than CLM (e.g., :math:`\sim`\ 5 km rather than :math:`\sim`\ 100 km). To improve the downscaling from CLM's grid to the ice sheet grid, the glaciated portion of each grid cell is divided into multiple elevation classes (section :numref:`Multiple elevation class scheme`). The above quantities are computed separately in each elevation class. The CESM coupler then computes high-resolution quantities via horizontal and vertical interpolation, and passes these high-resolution quantities to CISM.

There are several reasons for computing the SMB in CLM rather than in CISM:

#. It is much cheaper to compute the SMB in CLM for :math:`\sim`\ 10 elevation classes than in CISM. For example, suppose we are running CLM at a resolution of :math:`\sim`\ 50 km and CISM at :math:`\sim`\ 5 km. Greenland has dimensions of about 1000 x 2000 km. For CLM we would have 20 x 40 x 10 = 8,000 columns, whereas for CISM we would have 200 x 400 = 80,000 columns.

#. We can use the sophisticated snow physics parameterization already in CLM instead of implementing a separate scheme for CISM. Any improvements to CLM are applied to ice sheets automatically.

#. The atmosphere model can respond during runtime to ice-sheet surface changes (even in the absence of two-way feedbacks with CISM). As shown by :ref:`Pritchard et al. (2008)<Pritchardetal2008>`, runtime albedo feedback from the ice sheet is critical for simulating ice-sheet retreat on paleoclimate time scales. Without this feedback the atmosphere warms much less, and the retreat is delayed.

#. The improved SMB is potentially available in CLM for all glaciated grid cells (e.g., in the Alps, Rockies, Andes, and Himalayas), not just those which are part of ice sheets.

In typical runs, CISM is not evolving; CLM computes the SMB and sends it to CISM, but CISM's ice sheet geometry remains fixed over the course of the run. In these runs, CISM serves two roles in the system:

#. Over the CISM domain (typically Greenland in CESM2), CISM dictates glacier areas and topographic elevations, overriding the values on CLM's surface dataset. CISM also dictates the elevation of non-glacier land units in its domain, and only in this domain are atmospheric fields downscaled to non-glacier land units. (So if you run with a stub glacier model - SGLC - then glacier areas and elevations will be taken entirely from CLM's surface dataset, and no downscaling will be done over non-glacier land units.)

#. CISM provides the grid onto which SMB is downscaled. (If you run with SGLC then SMB will still be computed in CLM, but it won't be downscaled to a high-resolution ice sheet grid.)

It is also possible to run CESM with an evolving ice sheet. In this case, CLM responds to CISM's evolution by adjusting the areas of the glacier land unit and each elevation class within this land unit, as well as the mean topographic heights of each elevation class. Thus, CLM's glacier areas and elevations remain in sync with CISM's. Conservation of mass and energy is done as for other landcover change (see Chapter :numref:`rst_Transient Landcover Change`).

.. _Glacier regions:

Glacier regions and their behaviors
-----------------------------------

The world's glaciers and ice sheets are broken down into a number of different regions (three by default) that differ in three respects:

#. Whether the gridcell's glacier land unit contains:

   a. Multiple elevation classes (section :numref:`Multiple elevation class scheme`)

   b. Multiple elevation classes plus virtual elevation classes 

   c. Just a single elevation class whose elevation matches the atmosphere's topographic height (so there is no adjustment in atmospheric forcings due to downscaling).

#. Treatment of glacial melt water:

   a. Glacial melt water runs off and is replaced by ice, thus keeping the column always frozen. In the absence of a dynamic ice sheet model, this behavior implicitly assumes an infinite store of glacial ice that can be melted (with appropriate adjustments made to ensure mass and energy conservation). This behavior is discussed in more detail in section :numref:`Computation of the surface mass balance`.

   b. Glacial melt water remains in place until it refreezes - possibly remaining in place indefinitely if the glacier column is in a warm climate. With this behavior, ice melt does not result in any runoff. Regions with this behavior cannot compute SMB, because negative SMB would be meaningless (due to the liquid water on top of the ice column). This behavior produces less realistic glacier physics. However, it avoids the negative ice runoff that is needed for the "replaced by ice" behavior to conserve mass and energy (as described in section :numref:`Computation of the surface mass balance`). Thus, in regions where CLM has glaciers but the atmospheric forcings are too warm to sustain those glaciers, this behavior avoids persistent negative ice runoff. This situation can often occur for mountain glaciers, where topographic smoothing in the atmosphere results in a too-warm climate. There, avoiding persistent negative ice runoff can be more important than getting the right glacier ice physics.

#. Treatment of ice runoff from snow capping (as described in section :numref:`Runoff from glaciers and snow-capped surfaces`). Note that this is irrelevant in regions with an evolving, two-way-coupled ice sheet (where the snow capping term is sent to CISM rather than running off):

   a. Ice runoff from snow capping remains ice. This is a crude parameterization of iceberg calving, and so is appropriate in regions where there is substantial iceberg calving in reality.

   b. Ice runoff from snow capping is melted (generating a negative sensible heat flux) and runs off as liquid. This matches the behavior for non-glacier columns. This is appropriate in regions that have little iceberg calving in reality. This can be important to avoid unrealistic cooling of the ocean and consequent runaway sea ice growth.

The default behaviors for the world's glacier and ice sheet regions are described in :numref:`Table Glacier region behaviors`. Note that the Greenland region stops at the edge of Greenland as defined by CISM. This means that, by default, SMB is not computed for grid cells outside Greenland but within the CISM domain. (This treatment of the non-Greenland portion of the CISM domain as being the same as the world's mountain glaciers rather than like Greenland itself is mainly for the sake of avoiding unrealistic fluxes from the Canadian archipelago that can potentially result in runaway sea ice growth in that region.)

.. _Table Glacier region behaviors:

.. table:: Glacier region behaviors

 +---------------+---------------+---------------+---------------+
 | Region        | Elevation     | Glacial melt  | Ice runoff    |
 |               | classes       |               |               |
 +===============+===============+===============+===============+
 | Greenland     | Virtual       | Replaced by   | Remains ice   |
 |               |               | ice           |               |
 +---------------+---------------+---------------+---------------+
 | Antarctica    | Multiple      | Replaced by   | Remains ice   |
 |               |               | ice           |               |
 +---------------+---------------+---------------+---------------+
 | All others    | Single        | Remains in    | Melted        |
 |               |               | place         |               |
 +---------------+---------------+---------------+---------------+

.. note::

It is possible to have non-virtual, non-SMB-computing areas within the CISM domain (as is the case for the portion of CISM's Greenland domain outside of Greenland itself). However, these areas will send 0 SMB and will not be able to adjust to CISM-dictated changes in glacier area. Therefore, it is best to set up the glacier regions and their behaviors so that as much of the CISM domain as possible is covered by virtual, SMB-computing areas.

.. note::

 The combination of the ``Glacial melt = Replaced by ice`` and the ``Ice runoff = Melted`` behaviors results in particularly non-physical behavior: During periods of glacial melt, a negative ice runoff is generated (due to the ``Glacial melt = Replaced by ice`` behavior); this negative ice runoff is converted to a negative liquid runoff plus a positive sensible heat flux (due to the ``Ice runoff = Melted`` behavior). The net result is zero runoff but a positive sensible heat flux generated from glacial melt. Because of how physically unrealistic this is, CLM does not allow this combination of behaviors.

.. _Multiple elevation class scheme:

Multiple elevation class scheme
-------------------------------

The glacier land unit contains multiple columns based on surface elevation. These are known as elevation classes, and the land unit is referred to as *glacier\_mec*. (As described in section :numref:`Glacier regions`, some regions have only a single elevation class, but they are still referred to as *glacier\_mec* land units.) The default is to have 10 elevation classes whose lower limits are 0, 200, 400, 700, 1000, 1300, 1600, 2000, 2500, and 3000 m. Each column is characterized by a fractional area and surface elevation that are read in during model initialization, and then possibly overridden by CISM as the run progresses. Each *glacier\_mec* column within a grid cell has distinct ice and snow temperatures, snow water content, surface fluxes, and SMB.

The atmospheric surface temperature, potential temperature, specific humidity, density, and pressure are downscaled from the atmosphere's mean grid cell elevation to the *glacier\_mec* column elevation using a specified lapse rate (typically 6.0 deg/km) and an assumption of uniform relative humidity. Longwave radiation is downscaled by assuming a linear decrease in downwelling longwave radiation with increasing elevation (0.032 W m\ :sup:`-2` m\ :sup:`-1`, limited to 0.5 - 1.5 times the gridcell mean value, then normalized to conserve gridcell total energy) :ref:`(Van Tricht et al., 2016)<VanTrichtetal2016>`. Total precipitation is partitioned into rain vs. snow as described in Chapter :numref:`rst_Surface Characterization, Vertical Discretization, and Model Input Requirements`. The partitioning of precipitation is based on the downscaled temperature, allowing rain to fall at lower elevations while snow falls at higher elevations.

This downscaling allows lower-elevation columns to undergo surface melting while columns at higher elevations remain frozen. This gives a more accurate simulation of summer melting, which is a highly nonlinear function of air temperature.

Within the CISM domain, this same downscaling procedure is also applied to all non-urban land units. The elevation of non-glacier land units is taken from the mean elevation of ice-free grid cells in CISM. This is done in order to keep the glaciated and non-glaciated portions of the CISM domain as consistent as possible.

In contrast to most CLM subgrid units, glacier\_mec columns can be active (i.e., have model calculations run there) even if their area is zero. These are known as "virtual" columns. This is done because the ice sheet model may require a SMB for some grid cells where CLM has zero glacier area in that elevation range. Virtual columns also facilitate glacial advance and retreat in the two-way coupled case. Virtual columns do not affect energy exchange between the land and the atmosphere.

.. _Computation of the surface mass balance:

Computation of the surface mass balance
---------------------------------------

This section describes the computation of surface mass balance and associated runoff terms. The description here only applies to regions where glacial melt runs off and is replaced by ice, not to regions where glacial melt remains in place. Thus, by default, this only applies to Greenland and Antarctica, not to mountain glaciers elsewhere in the world. (See also section :numref:`Glacier regions`.)

The SMB of a glacier or ice sheet is the net annual accumulation/ablation of mass at the upper surface. Ablation is defined as the mass of water that runs off to the ocean. Not all the surface meltwater runs off; some of the melt percolates into the snow and refreezes. Accumulation is primarily by snowfall and deposition, and ablation is primarily by melting and evaporation/sublimation. CLM uses a surface-energy-balance (SEB) scheme to compute the SMB. In this scheme, the melting depends on the sum of the radiative, turbulent, and conductive fluxes reaching the surface, as described elsewhere in this document.

Note that the SMB typically is defined as the total accumulation of ice and snow, minus the total ablation. The SMB flux passed to CISM is the mass balance for ice alone, not snow. We can think of CLM as owning the snow, whereas CISM owns the underlying ice. Fluctuations in snow depth between 0 and 10 m water equivalent are not reflected in the SMB passed to CISM. In transient runs, this can lead to delays of a few decades in the onset of accumulation or ablation in a given glacier column.

SMB is computed and sent to the CESM coupler regardless of whether and where CISM is operating. However, the effect of SMB terms on runoff fluxes differs depending on whether and where CISM is evolving in two-way-coupled mode. This is described by the variable *glc\_dyn\_runoff\_routing*. (This is real-valued in the code to handle the edge case where a CLM grid cell partially overlaps with the CISM grid, but we describe it as a logical variable here for simplicity.) In typical cases where CISM is not evolving, *glc\_dyn\_runoff\_routing* will be false everywhere; in these cases, CISM's mass is not considered to be part of the coupled system. In cases where CISM is evolving and sending its own calving flux to the coupler, *glc\_dyn\_runoff\_routing* will be true over the CISM domain and false elsewhere.

Any snow capping (section :numref:`Runoff from glaciers and snow-capped surfaces`) is added to :math:`q_{ice,frz}`. Any liquid water (i.e., melted ice) below the snow pack in the glacier column is added to :math:`q_{ice,melt}`, then is converted back to ice to maintain a pure-ice column. Then the total SMB is given by :math:`q_{ice,tot}`:

.. math::
   :label: 13.1

   q_{ice,tot} = q_{ice,frz} - q_{ice,melt}

CLM is responsible for generating glacial surface melt, even when running with an evolving ice sheet. Thus, :math:`q_{ice,melt}` is always added to liquid runoff (:math:`q_{rgwl}`), regardless of *glc\_dyn\_runoff\_routing*. However, the ice runoff flux depends on *glc\_dyn\_runoff\_routing*. If *glc\_dyn\_runoff\_routing* is true, then CISM controls the fate of the snow capping mass in :math:`q_{ice,frz}` (e.g., eventually transporting it to lower elevations where it can be melted or calved). Since CISM will now own this mass, the snow capping flux does *not* contribute to any runoff fluxes generated by CLM in this case.

If *glc\_dyn\_runoff\_routing* is false, then CLM sends the snow capping flux as runoff, as a crude representation of ice calving (see also sections :numref:`Runoff from glaciers and snow-capped surfaces` and :numref:`Glacier regions`). However, this ice runoff flux is reduced by :math:`q_{ice,melt}`. This reduction is needed for conservation; its need is subtle, but can be understood with either of these explanations:

- When ice melts, we let the liquid run off and replace it with new ice. That new ice needs to come from somewhere to keep the coupled system in water balance. We "request" the new ice from the ocean by generating a negative ice runoff equivalent to the amount we have melted.

- Ice melt removes mass from the system, as it should. But the snow capping flux also removes mass from the system. The latter is a crude parameterization of calving, assuming steady state - i.e., all ice gain is balanced by ice loss. This removal of mass due to both accumulation and melt represents a double-counting. Each unit of melt indicates that one unit of accumulation should not have made it to the ocean as ice, but instead melted before it got there. So we need to correct for this double-counting by removing one unit of ice runoff for each unit of melt.

For a given point in space or time, this reduction can result in negative ice runoff. However, when integrated over space and time, for an ice sheet that is near equilibrium, this just serves to decrease the too-high positive ice runoff from snow capping. (The treatment of snow capping with *glc\_dyn\_runoff\_routing* false is based on this near-equilibrium assumption - i.e., that ice accumulation is roughly balanced by :math:`calving + melt`, integrated across space and time. For glaciers and ice sheets that violate this assumption, either because they are far out of equilibrium with the climate or because the model is being run for hundreds of years, there are two ways to avoid the unrealistic ice runoff from snow capping: by running with an evolving, two-way-coupled ice sheet or by changing a glacier region's ice runoff behavior as described in section :numref:`Glacier regions`.)

In regions where SMB is computed for glaciers, SMB is also computed for the natural vegetated land unit. Because there is no ice to melt in this land unit, it can only generate a zero or positive SMB. A positive SMB is generated once the snow pack reaches its maximum depth. When running with an evolving ice sheet, this condition triggers glacial inception.

