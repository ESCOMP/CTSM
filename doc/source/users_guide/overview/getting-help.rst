.. include:: ../substitutions.rst

.. _getting-help:

==============
 Getting Help
==============
In addition to this users-guide there are several other resources that are available to help you use |version|. The first one is the |cesmrelease| User's-Guide, which documents the entire process of creating cases with |cesmrelease|. And next is the CIME User's Guide which goes over the scripts and infrastructure used for running |version| in |cesmrelease|. The CESM bulletin board which is a web-site for exchanging information between users of CESM. There are also CLM web-pages specific for CLM, and finally there is an email address to report bugs that you find in |cesmrelease|.

---------------------------
The CESM User's-Guide
---------------------------
|release|
|release| in |cesmrelease| is always run from within the standard |cesmrelease| build and run scripts. Therefore, the user of |version| should familiarize themselves with the |cesmrelease| scripts and understand how to work with them. User's-Guide documentation on the |cesmrelease| scripts are available from the following web-page. The purpose of this |version| in |cesmrelease| User's Guide is to give the |version| user more complete details on how to work with CLM and the set of tools that support CLM, as well as to give examples that are unique to the use of CLM. However, the |cesmrelease| Scripts User's-Guide remains the primary source to get detailed information on how to build and run the CESM system.

`|cesmrelease| Quickstart Guide <https://escomp.github.io/cesm/release-cesm2/>`_

---------------------------
The CIME User's-Guide
---------------------------

The CIME Users'-Guide goes into the how to use the scripts and infrastructure of the CESM. `CIME Users Guide <http://esmci.github.io/cime/>`_

-----------------------
The CESM Bulletin Board
-----------------------

There is a rich and diverse set of people that use the CESM, and often it is useful to be in contact with others to get help in solving problems or trying something new. To facilitate this we have an online Bulletin Board for questions on the CESM. There are also different sections in the Bulletin Board for the different component models or for different topics.

`CTSM Forum <https://bb.cgd.ucar.edu/cesm/forums/ctsm-clm-mosart-rtm.134/>`_

`All CESM Forums (including forums for infrastructure/porting questions, etc.) <http://bb.cgd.ucar.edu/cesm/>`_

-----------------
The CLM web pages
-----------------

The main `CLM web page <http://www.cgd.ucar.edu/tss/clm/>`_ contains information on the CLM, its history, developers, as well as downloads for previous model versions. Some other links are available at the `CESM2 land component webpage <http://www.cesm.ucar.edu/models/cesm2/land>`. There are also documentation text files in the `$CTSMROOT/doc directory <https://github.com/ESCOMP/CTSM/tree/master/doc>`_ that give some quick information on using CLM.

Also note that several of the XML database files can help with namelist options, namelist defaults, or compsets. For the most recent release:

- `$CTSMROOT/bld/namelist_files/namelist_definition_ctsm.xml <https://github.com/ESCOMP/CTSM/blob/master/bld/namelist_files/namelist_definition_ctsm.xml>`_ -- definition of latest CTSM namelist items.
- `$CTSMROOT/bld/namelist_files/namelist_defaults_ctsm.xml <https://github.com/ESCOMP/CTSM/blob/master/bld/namelist_files/namelist_defaults_ctsm.xml>`_ -- default values for latest CTSM namelist items.
- `$CTSMROOT/cime_config/config_component.xml <https://github.com/ESCOMP/CTSM/blob/master/cime_config/config_component.xml>`_ -- definition of all the CLM specific XML variables.
- `$CTSMROOT/cime_config/config_compsets.xml <https://github.com/ESCOMP/CTSM/blob/master/cime_config/config_compsets.xml>`_ -- definition of all the CLM compsets.

Some archives are available for previous versions:

- `Archive of namelist_definition_clm4_0.xml <https://github.com/ESCOMP/CTSM/blob/clm5.0.000/bld/namelist_files/namelist_definition_clm4_0.xml>`_ -- definition of CLM4.0 namelist items.
- `Archive of namelist_definition_clm4_5.xml <https://github.com/ESCOMP/CTSM/blob/clm5.0.000/bld/namelist_files/namelist_definition_clm4_5.xml>`_ -- definition of CLM4.5/CLM5.0 namelist items.
- `Archive of namelist_defaults_clm4_0.xml <https://github.com/ESCOMP/CTSM/blob/clm5.0.000/bld/namelist_files/namelist_defaults_clm4_0.xml>`_ -- default values for CLM4.0 namelist items.
- `Archive of namelist_defaults_clm4_5.xml <https://github.com/ESCOMP/CTSM/blob/clm5.0.000/bld/namelist_files/namelist_defaults_clm4_5.xml>`_ -- default values for CLM4.5/CLM5.0 namelist items.


----------------------------
Reporting bugs in |version|
----------------------------

If you have any problems, additional questions, bug reports, or any other feedback, please report it as an issue on GitHub https://github.com/ESCOMP/ctsm/issues or for CIME scripts and infrastructure to https://github.com/ESMCI/CIME/issues. Or send an email to <`cesmhelp@cgd.ucar.edu <cesmhelp@cgd.ucar.edu>`_> or <`ctsm-software@ucar.edu <ctsm-software@ucar.edu>`_>. If you find bad, wrong, or misleading information in this users guide report it as an issue on CTSM.

.. _acronyms-and-terms:

---------------------------------------
Some Acronym's and Terms We'll be Using
---------------------------------------

CAM
  Community Atmosphere Model (CAM). The prognostically active atmosphere model component of CESM.

CESM
  Community Earth System Model (CESM). The coupled earth system model that CLM is a component of.

CIME
  The Common Infrastructure for Modeling the Earth (CIME - pronounced "SEAM") provides a Case Control System for configuring, compiling and executing Earth system models, data and stub model components, a driver and associated tools and libraries.

CLM
  Community Land Model (CLM). The prognostically active land model component of CESM.

CLMBGC
  Community Land Model (|version|) with BGC Biogeochemistry. Uses CN Biogeochemistry with vertically resolved soil Carbon, CENTURY model like pools, and Nitrification/De-Nitrification. The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="phys clm5_0 -bgc bgc``

CLMBGC-Crop
  Community Land Model (|version|) with BGC Biogeochemistry and prognotic crop.  The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="phys clm5_0 -bgc bgc -crop``

CLMCN
  Community Land Model (CLM) with Carbon Nitrogen (CN) Biogeochemistry (either CLM4.0, CLM4.5 or |version|) The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-bgc cn" -append``

CLMSP
  Community Land Model (CLM) with Satellite Phenology (SP) (either CLM4.0, CLM4.5 or |version|) The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-bgc sp" -append``

CLMU
  Community Land Model (CLM) Urban Model (either CLM4.0, CLM4.5 or |version|). The urban model component of CLM is ALWAYS active (unless you create special surface datasets that have zero urban percent, or for regional/single-point simulations for a non-urban area).

CRUNCEP
  The Climate Research Unit (CRU) analysis of the NCEP atmosphere reanalysis atmosphere forcing data. This can be used to drive CLM with atmosphere forcing from 1901 to 2016. This data is updated every year, the version we are currently using is Version-7. The las CESM1.2.2 release used Version-4 data.

CTSM
  The Community Terrestrial Systems Model, of which |version| and CLM4.5 are namelist option sets of. CTSM is a wider community
  that includes using CTSM for Numerical Weather Prediction (NWP) as well as climate.

DATM
  Data Atmosphere Model (DATM) the prescribed data atmosphere component for CESM. Forcing data that we provide are either the CRUNCEP, Qian, or GSWP3 forcing datasets (see below).

DV
  Dynamic global vegetation, where fractional PFT (see PFT below) changes in time prognostically. Can NOT be used with prescribed transient PFT (requires either CLMBGC or CLMCN for either CLM4.0, CLM4.5 or |version|). The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-bgc cndv" -append``

  This option is being phased out for the different methodology of FATES (see below). DV is not currently scientifically validated
  and as such should be considered experimental.

ESMF
  Earth System Modeling Framework (ESMF). They are a software project that provides a software library to support Earth System modeling. We provide interfaces for ESMF as well as use their regridding capabilities for offline CLM tools.

FATES
  Functionally Assembled Terrestrial Ecosystem Simulator. This is being developed by the Next Generation Ecosystem Experiment Tropics' (NGEE-T)
  project and uses both |version| and the land model component of E3SM (Energy Exascale Earth System Model).

FUN
  Fixation and Uptake of Nitrogen model, a parameter option of |version|.

GSWP3
  Global Soil Wetness Project (GSPW3) atmospheric forcing data. It is a 3-hourly 0.5° global forcing product (1901-2014) that is based on the NCEP 20th Century Reanalysis, with additional bias corrections added by GSWP3.

LUNA
  Leaf Utilization of Nitrogen for Assimilation parameterization option as part of |version|.

NCAR
  National Center for Atmospheric Research (NCAR). This is the research facility that maintains CLM with contributions from other national labs and Universities.

NCEP
  The National Center for Environmental Prediction (NCEP). In this document this normally refers to the reanalysis atmosphere data produced by NCEP.

MOSART
  Model for Scale Adaptive River Transport, ROF model component option added as part of |version|. It is the standard
  ROF model used in |version| compsets.

PFT
  Plant Function Type (PFT). A type of vegetation that CLM parameterizes.

ROF
  River runOff Model to route flow of surface water over land out to the ocean. |cesmrelease| has two components options for this
  the new model MOSART and previous model RTM.

RTM
  River Transport Model, ROF model component option that has been a part of all versions of CESM. It is the standard ROF
  model used in CLM4.5 and CLM4.0 compsets.

SCRIP
  Spherical Coordinate Remapping and Interpolation Package (SCRIP). We use it's file format for specifying both grid coordinates as well as mapping between different grids.

VIC
  Variable Infiltration Capacity (VIC) model for hydrology. This is an option to |version| in place of the standard |version| hydrology. The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-vichydro on" -append``
