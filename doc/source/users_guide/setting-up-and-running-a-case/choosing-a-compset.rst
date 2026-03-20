.. include:: ../substitutions.rst

.. _choosing-a-compset:

====================
Choosing a compset
====================

When setting up a new case one of the first choices to make is which "component set" (or compset) to use. The compset refers to which component models are used as well as specific settings for them. We label the different types of compsets with a different letter of the alphabet from "A" (for all data model) to "X" (for all dead model). The compsets of interest when working with CLM are the "I" compsets (which contain CLM with a data atmosphere model and a stub ocean, and stub sea-ice models), "E" and "F" compsets (which contain CLM with the active atmosphere model (CAM), prescribed sea-ice model, and a data ocean model), and "B" compsets which have all active components. Below we go into details on the "I" compsets which emphasize CLM as the only active model, and just mention the two other categories.

To run CLM coupled to CAM ("E" or "F" compsets) or fully coupled ("B" compsets) you need to be running CLM from a CESM checkout rather than a CTSM checkout (see :ref:`ctsm_vs_cesm_checkout`).

When working with CLM you usually want to start with a relevant "I" compset before moving to the more complex cases that involve other active model components. The "I" compsets can exercise CLM in a way that is similar to the coupled modes, but with much lower computational cost and faster turnaround times.

Compsets coupled to data atmosphere and stub ocean/sea-ice ("I" compsets)
-------------------------------------------------------------------------

The entire list of compsets available out-of-the-box can be browsed at the `CESM Component Set Definitions <https://docs.cesm.ucar.edu/models/cesm2/config/compsets.html>`_ page. Note that using the compset longnames, even more combinations are possible than those listed. That webpage also includes information on whether each compset has been tested and/or scientifically validated.

To get a list of the compsets use the ``query_config`` command as follows:
::

    $CTSMROOT/cime/scripts/query_config --compsets clm

Compsets with different choices for the River "Runoff OutFlow" (ROF) Model
--------------------------------------------------------------------------

CTSM can be run with four different options for the ROF model: stub, MOSART, RTM, or the newly available component mizuRoute. The default for compsets is MOSART and as such it isn't mentioned in the compset aliases. Compsets with the stub ROF model usually have a "Rs" in the name to designate that a stub ROF is being used (for example the single point tower site compset I1PtClm60SpRs).
Compset aliases for MOSART as it's the default don't have it in the name. Compsets with clm4_5 physics are with RTM as RTM was the default ROF model when CLM4.5 was created. Also since Paleo climate work uses RTM, the "NoAnthro" compset aliases use RTM. Compset aliases with "Miz" in the name use mizuRoute.

Both MOSART and RTM run on a default half degree grid, that is selected as part of the standard grid aliases. mizuRoute also can use the standard grid aliases and will be run on it's half degree grid if so. See the next section for on the other options for mizuRoute grids.

Compsets and grids with the mizuRoute ROF model
-----------------------------------------------

Compset aliases with "Miz" in the name use mizuRoute as the ROF model. For example, the compset alias I2000Clm60SpMizGs which is for present day with clm6_0 physics using Satellite Phenology with a stub glacier model.
As mizuRoute is a ROF model grid alias for special grids for it include a "_r*" in the middle of the compset name (between the atmosphere/land grid and the ocean grid/mask).

To get a list of the grid alises use the ``query_config`` command as follows:
::

    $CTSMROOT/cime/scripts/query_config --grids

There are five mizuRoute ROF grids available (r05, rUSGS, rHDMA, rHDMAlk, rMERIT)

- r05 is the default half degree grid
- rUSGS is the lowest resolution mizuRoute HRU grid and only covering Continental US
- rHDMA is the medium resolution mizuRoute HRU grid
- rHDMAlk is the medium resolution mizuRoute HRU grid that includes lakes
- rMERIT is the highest resolution mizuRoute HRU grid

- r05 is the same as the default grid for MOSART and RTM, but just over land.
- rUSGS is the USGS Geospatial Fabric for National Hydrologic Modeling. It has over 110k HRU's over the continental US.
- rHDMA is the HydroSHEDS Derived Medium Resolution Global Hydrological Model. It has about 300k HRU's.
- rHDMAlk is the HDMA grid with lakes added in. It has about 300k HRU's as well.
- rMERIT is the Multi-Error-Removed Improved-Terrain global hydrological model. It has about 3 million HRU's.

Unlike the other ROF models, the mizuRoute grids require mapping files in order to do the mapping in the coupler. This is because the grids have small overlaps and holes that the ESMF on the fly regridding currently has trouble with.
As such we supply a limited number of mapping files for common atmosphere/land resolutions (NLDAS i.e. continental US, half-degree, 1-degree, and 2-degree using the CRU half degree grid hcru, and the finite-volume grids f09 and f19).
The mizuRoute r05 and HDMA grids over the very limited 5x5 Amazon region are also available.

Compsets coupled to active atmosphere with data ocean
-----------------------------------------------------
CAM compsets are compsets that start with "E" or "F" in the name. They are described more fully in the scripts documentation or the CAM documentation. "E" compsets have a slab ocean model while "F" compsets have a data ocean model.

Fully coupled compsets with fully active ocean, sea-ice, and atmosphere
-----------------------------------------------------------------------
Fully coupled compsets are compsets that start with "B" in the name. They are described more fully in the scripts documentation.

Conclusion to choosing a compset
--------------------------------
We've introduced the basic type of compsets that use CLM and given some further details for the "standalone CLM" (or "I" compsets). `$CTSMROOT/cime_config/config_compsets.xml <CLM-https://github.com/ESCOMP/CTSM/blob/master/cime_config/config_compsets.xml>`_ lists all of the compsets and gives a full description of each of them. In the next section we look into customizing the setup time options for compsets using CLM.
