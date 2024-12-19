# What's new in CTSM 5.3 (tag `ctsm5.3.0`)

## Purpose and description of changes since CTSM 5.2 (tag `ctsm5.2.005`)
- Adds CN Matrix method to speed up spinup for the BGC model.
- Updates surface datasets.
- Brings in new Leung 2023 dust emissions.
- Adds explicit air conditioning for the urban model.
- Updates crop calendars.
- Updates fire model with various improvements, including default parameterization against climate reanalysis from CRU-JRA instead of GSWP.
- FATES compsets can now be run with transient land use.

These changes were needed ahead of the CESM3 capability/functionality "chill". For `clm6_0` physics these options are now turned on by default, in addition to Sturm snow and excess ice.

## Changes to CTSM infrastructure
 - `manage_externals` removed and replaced by `git-fleximod`
 - Ability to handle CAM7 in `LND_TUNING_MODE`

## Changes to CTSM answers

 Changes to defaults for `clm6_0` physics:
  - Urban explicit A/C turned on
  - Snow thermal conductivity method is now `Sturm_1997`
  - New initial conditions file for f09 ("1-degree") 1850, with more in progress
  - New crop calendars
  - Dust emissions method is now `Leung_2023`
  - Excess ice is turned on
  - Updates to MEGAN for BVOCs
  - Updates to BGC fire method

 Changes for all physics versions:
  - Parameter files updated
  - FATES parameter file updated
  - Glacier region 1 is now undefined
  - Update in FATES transient land use
  - Pass active glacier (CISM) runoff directly to river model (MOSART)
  - Add the option for using Matrix CN method for Carbon/Nitrogen BGC spinup

## New surface datasets

- With new surface datasets the following GLC fields have region "1" set to UNSET: glacier_region_behavior, glacier_region_melt_behavior, glacier_region_ice_runoff_behavior
- Updates to allow creating transient landuse timeseries files going back to 1700.
- Fix an important bug on soil fields that was there since `ctsm5.2.0`. This has the side effect of `mksurfdata_esmf` now giving identical answers with a change in number of processors, as it should.
- Surface datasets now provided for the `ne0np4.POLARCAP.ne30x4` grid.
- Surface datasets now have their version number embedded to prevent mismatch of surface dataset and CTSM version.
- Remove the `--hires_pft` option from `mksurfdata_esmf` as we don't have the datasets for it.
- Remove `VIC` fields from surface datasets.
- Updates to input datasets in PFT/LAI/soil-color raw datasets (now from the TRENDY2024 timeseries that ends in 2023), as well as two fire datasets (crop fire peak month, peatland fraction), and the glacier behavior dataset.

