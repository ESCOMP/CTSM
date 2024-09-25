Purpose and description of changes since ctsm5.2.005
----------------------------------------------------

Bring in updates needed for the CESM3.0 science capability/functionality "chill". Most importantly bringing
in: CN Matrix to speed up spinup for the BGC model, updated surface datasets, updated Leung 2023 dust emissions,
explicit Air Conditioning for the Urban model, updates to crop calendars. For clm6_0 physics these options are now
default turned on in addition to Sturm snow, and excess ice.

Changes to CTSM Infrastructure:
===============================

 - manage_externals removed and replaced by git-fleximod
 - Ability to handle CAM7 in LND_TUNING_MODE

Changes to CTSM Answers:
========================

 Changes to defaults for clm6_0 physics:
  - Urban explicit A/C turned on
  - Snow thermal conductivity is now Sturm_1997
  - New IC file for f09 1850
  - New crop calendars
  - Dust emissions is now Leung_2023
  - Excess ice is turned on
  - Updates to MEGAN for BVOC's
  - Updates to BGC fire method

 Changes for all physics versions:

  - Parameter files updated
  - FATES parameter file updated
  - Glacier region 1 is now undefined
  - Update in FATES transient Land use
  - Pass active glacier (CISM) runoff directly to river model (MOSART)
  - Add the option for using matrix for Carbon/Nitrogen BGC spinup

New surface datasets:
=====================

- With new surface datasets the following GLC fields have region "1" set to UNSET:
     glacier_region_behavior, glacier_region_melt_behavior, glacier_region_ice_runoff_behavior
- Updates to allow creating transient landuse timeseries files going back to 1700.
- Fix an important bug on soil fields that was there since ctsm5.2.0. This results in mksurfdata_esmf now giving identical answers with a change in number of processors, as it should.
- Add in creation of ne0np4.POLARCAP.ne30x4 surface datasets.
- Add version to the surface datasets.
- Remove the --hires_pft option from mksurfdata_esmf as we don't have the datasets for it.
- Remove VIC fields from surface datasets.

New input datasets to mksurfdata_esmf:
======================================

- Updates in PFT/LAI/soil-color raw datasets (now from the TRENDY2024 timeseries that ends in 2023), as well as two fire datasets (AG fire, peatland), and the glacier behavior dataset.

