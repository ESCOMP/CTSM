.. _scientific-validiation:

========================
 Scientific Validation
========================

In this section we go over what has been extensively tested and scientifically validated with |version|, and maybe more importantly what has NOT been tested and may NOT be scientifically validated. You can use all features of CLM, but need to realize that some things haven't been tested extensively or validated scientifically. When you use these features you may run into trouble doing so, and will need to do your own work to make sure the science is reasonable.

--------------------------------------------------------------
Standard Configuration and Namelist Options that are Validated
--------------------------------------------------------------

See 
`http://www.cesm.ucar.edu/models/cesm1.2/clm/CLM_configurations_CESM1.2.pdf <http://www.cesm.ucar.edu/models/cesm1.2/clm/CLM_configurations_CESM1.2.pdf>`_ for an explanation of what configurations are scientifically validated for |version|. For CLM4.0 changes to the science of the model are minimal since CESM1.1.1 so we expect answers to be very similar to using it.

In the sections below we go through configuration and/or namelist options or modes that the user should be especially wary of using. You are of course free to use these options, and you may find that they work functionally. Although in some cases you will find issues even with functionality of using them. If so you will need to test, debug and find solutions for these issues on your own. But in every case you will need to go through more extensive work to validate these options from a scientific standpoint. Some of these options are only for |version| while others are for both CLM4.0 AND |version| we explicitly say which they apply to.

---------------------------------------------------------------------------------------------------------------
Configure Modes NOT scientifically validated, documented, supported or, in some cases, even advised to be used:
---------------------------------------------------------------------------------------------------------------

These are options that you would add to ``CLM_CONFIG_OPTS``.

1. exlaklayers on[|version| only] This mode is NOT tested and may NOT be even functional.

2. snicar_frc on[CLM4.0 AND |version|] This mode is tested and functional, but is NOT constantly scientifically validated, and should be considered experimental.

3. vichydro on[|version| only] This mode is tested and functional, but does NOT have long scientific validation simulations run with it so, should be considered experimental.

4. vsoilc_centbgc[|version| only] This option is extensively tested for both "on" and "off". The "no-vert" option has limited testing performed on it, but isn't scientifically validated (and it currently has a bug -- see 1746 and 1672 in `$CTSMROOT/doc/KnownBugs <CLM-URL>`_). The "no-cent" and "no-nitrif" options are NOT tested and as such may NOT ben even functional.

----------------------------------------------
Namelist options that should NOT be exercised:
----------------------------------------------

----------------------------------------------------
Build-Namelist options that should NOT be exercised:
----------------------------------------------------
1. -irrig with -bgc cn and -phys clm4_0 We have only run the irrigation model with CLMSP (i.e. without the CN model). We recommend that if you want to run the irrigation model with CN, that you do a spinup. But, more than that you may need to make adjustments to irrig_factor in $CTSMROOT/src/biogeophys/CanopyFluxesMod.F90. See the notes on this in the description of the irrigation model in the 
`Technical Descriptions of the Interactive Crop Management and Interactive Irrigation Models <CLM-URL>`_.

2. -irrig with -crop on and -phys clm4_0 Irrigation doesn't work with the prognostic crop model. Irrigation is only applied to generic crop currently, which negates it's practical usage. We also have a known problem when both are on (see bug 1326 in the `$CTSMROOT/doc/KnownBugs <CLM-URL>`_ file). If you try to run in this mode, the CLM build-namelist will return with an error.

--------------------------------------------
Namelist items that should NOT be exercised:
--------------------------------------------

suplnitro='ALL' The suplnitro namelist option to the CN Biogeochemistry model supplies unlimited nitrogen and therefore vegetation is over-productive in this mode.

urban_traffic:Not currently functional

allowlakeprod:Considered experimental.

anoxia_wtsat:Considered experimental (deprecated will be removed).

atm_c14_filename:Considered experimental (dataset not provided).

exponential_rooting_profile:Considered experimental.

fin_use_fsat:Considered experimental.

glc_dyntopo:Not currently functional.

lake_decomp_fact:Considered experimental.

more_vertlayers:Considered experimental.

no_frozen_nitrif_denitrif:Considered experimental.

perchroot:Considered experimental.

perchroot_alt:Considered experimental.

replenishlakec:Considered experimental.

use_c14_bombspike:Considered experimental (dataset not provided).

usefrootc:Considered experimental.
