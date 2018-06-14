.. _getting-help:

==============
 Getting Help
==============
In addition to this users-guide there are several other resources that are available to help you use CLM5.0. The first one is the CESM1.2.0 User's-Guide, which documents the entire process of creating cases with CESM1.2.0. The next is the CESM bulletin board which is a web-site for exchanging information between users of CESM. There are also CLM web-pages specific for CLM, and finally there is an email address to report bugs that you find in CESM1.2.0.

---------------------------
The CESM User's-Guide
---------------------------
+|release|
+|release| in +|cesmrelease| is always run from within the standard +|cesmrelease| build and run scripts. Therefore, the user of CLM4.5 should familiarize themselves with the +|cesmrelease| scripts and understand how to work with them. User's-Guide documentation on the +|cesmrelease| scripts are available from the following web-page. The purpose of this +|version| in +|cesmrelease| User's Guide is to give the +|version| user more complete details on how to work with CLM and the set of tools that support CLM, as well as to give examples that are unique to the use of CLM. However, the +|cesmrelease| Scripts User's-Guide remains the primary source to get detailed information on how to build and run the CESM system.

`+|cesmrelease| Quickstart Guide <https://escomp.github.io/cesm/release-cesm2/>`_

-----------------------
The CESM Bulletin Board
-----------------------

There is a rich and diverse set of people that use the CESM, and often it is useful to be in contact with others to get help in solving problems or trying something new. To facilitate this we have an online Bulletin Board for questions on the CESM. There are also different sections in the Bulletin Board for the different component models or for different topics.

`CESM Online Bulletin Board <http://bb.cgd.ucar.edu/>`_

-----------------
The CLM web pages
-----------------

The main CLM web page contains information on the CLM, it's history, developers, as well as downloads for previous model versions. There are also documentation text files in the models/lnd/clm/doc directory that give some quick information on using CLM.

`CLM web page <http://www.cgd.ucar.edu/tss/clm/>`_
`CLM Documentation Text Files <CLM-URL>`_

Also note that several of the XML database files can be viewed in a web browser to get a nice table of namelist options, namelist defaults, or compsets. Simply view them as a local file and bring up one of the following files:

- `models/lnd/clm/bld/namelist_files/namelist_definition_clm4_0.xml <CLM-URL>`_ -- definition of CLM4.0 namelist items.
- `models/lnd/clm/bld/namelist_files/namelist_definition_clm4_5.xml <CLM-URL>`_ -- definition of CLM4.0 namelist items.
- `models/lnd/clm/bld/namelist_files/namelist_defaults_clm4_0.xml <CLM-URL>`_ -- default values for CLM4.0 namelist items.
- `models/lnd/clm/bld/namelist_files/namelist_defaults_clm4_5.xml <CLM-URL>`_ -- default values for CLM4.5 namelist items.
- `scripts/ccsm_utils/Case.template/config_definition.xml <CLM-URL>`_ -- definition of all env_*.xml items.
- `scripts/ccsm_utils/Case.template/config_compsets.xml <CLM-URL>`_ -- definition of all the compsets.
- `models/lnd/clm/bld/namelist_files/history_fields_clm4_0.xml <CLM-URL>`_ -- definition of CLM4.0 history fields.
- `models/lnd/clm/bld/namelist_files/history_fields_clm4_5.xml <CLM-URL>`_ -- definition of CLM4.5 history fields.

------------------------
Reporting bugs in CLM4.5
------------------------

If you have any problems, additional questions, bug reports, or any other feedback, please send an email to <`cesmhelp@cgd.ucar.edu <cesmhelp@cgd.ucar.edu>`_>. If you find bad, wrong, or misleading information in this users guide send an email to <`erik@ucar.edu <mailto:erik@ucar.edu>`_>. The current list of known issues for CLM4.5 in CESM1.2.0 is in the models/lnd/clm/doc/KnownBugs file, and the list of issues for CESM1.2.0 is at... 
`http://www.cesm.ucar.edu/models/cesm1.2//tags/cesm1_2_0/#PROBLEMS <http://www.cesm.ucar.edu/models/cesm1.2//tags/cesm1_2_0/#PROBLEMS>`_.

---------------------------------------
Some Acronym's and Terms We'll be Using
---------------------------------------

CAM
  Community Atmosphere Model (CAM). The prognostically active atmosphere model component of CESM.

CESM
  Community Earth System Model (CESM). The coupled earth system model that CLM is a component of.

CLM
  Community Land Model (CLM). The prognostically active land model component of CESM.

CLMBGC
  Community Land Model (CLM4.5) with BGC Biogeochemistry. Uses CN Biogeochemistry with vertically resolved soil Carbon, CENTURY model like pools, and Nitrification/De-Nitrification. The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="phys clm4_5 -bgc cn -vsoilc_centbgc on -clm4me on"``

CLMCN
  Community Land Model (CLM) with Carbon Nitrogen (CN) Biogeochemistry (either CLM4.0 or CLM4.5) The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-bgc cn" -append``

CLMSP
  Community Land Model (CLM) with Satellite Phenology (SP) (either CLM4.0 or CLM4.5) The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-bgc none" -append``

CLMU
  Community Land Model (CLM) Urban Model (either CLM4.0 or CLM4.5). The urban model component of CLM is ALWAYS active (unless you create special surface datasets that have zero urban percent, or for regional/single-point simulations for a non-urban area).

CRUNCEP
  The Climate Research Unit (CRU) analysis of the NCEP atmosphere reanalysis atmosphere forcing data. This can be used to drive CLM with atmosphere forcing from 1901 to 2010. We also DO expect to be able to update this dataset beyond 2010 as newer data becomes available.

DATM
  Data Atmosphere Model (DATM) the prescribed data atmosphere component for CESM. Forcing data that we provide are either the Qian or CRUNCEP forcing datasets (see below).

DV
  Dynamic global vegetation, where fractional PFT (see PFT below) changes in time prognostically. Can NOT be used with prescribed transient PFT (requires either CLM4.5-BGC or CLMCN for either CLM4.0 or CLM4.5). The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-bgc cndv" -append``

ESMF
  Earth System Modeling Framework (ESMF). They are a software project that provides a software library to support Earth System modeling. We provide interfaces for ESMF as well as use their regridding capabilities for offline CLM tools.

NCAR
  National Center for Atmospheric Research (NCAR). This is the research facility that maintains CLM with contributions from other national labs and Universities.

NCEP
  The National Center for Environmental Prediction (NCEP). In this document this normally refers to the reanalysis atmosphere data produced by NCEP.

PFT
  Plant Function Type (PFT). A type of vegetation that CLM parameterizes.

PTCLM
  PoinT CLM (PTCLM) a python script that operates on top of CLM for CLM4.5 to run single point simulations for CLM.

Qian
  The Qian et. al. analysis of the NCEP forcing data. This can be used to drive CLM with atmosphere forcing from 1948 to 2004. We do NOT expect to be able to update this dataset beyond 2004.

SCRIP
  Spherical Coordinate Remapping and Interpolation Package (SCRIP). We use it's file format for specifying both grid coordinates as well as mapping between different grids.

VIC
  Variable Infiltration Capacity (VIC) model for hydrology. This is an option to CLM4.5 in place of the standard CLM4.5 hydrology. The CLM_CONFIG_OPTS option for this is

  ``./xmlchange CLM_CONFIG_OPTS="-vichydro on" -append``
