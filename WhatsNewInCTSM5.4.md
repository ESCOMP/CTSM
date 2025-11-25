# What's new in CTSM 5.4 (tag `ctsm5.4.0xx`)

# Purpose and description of changes since CTSM 5.3 (tag `ctsm5.3.021`)

## New features

* New surface datasets from CMIP7 data including PFT and urban distributions, land use transitions, population density, and atmospheric C isotopes.  These data are only available through the historical record (1850-2023), and  
  * are not available for future periods (presently known as SSP),  
  * for future periods and N deposition we continue to use CMIP6 data from CESM2.  
* Option to use CRUJRA2024 atmospheric driver data with clm6 and clm5 physics options ([PR #2956](https://github.com/ESCOMP/ctsm/pull/2956)), this is the default data-atmosphere (DATM) for clm6. This CRUJRA dataset covers 1901-2023, whereas previous GSWP3 only covers 1901-2014.  
* Capability to run single-point PLUMBER tower sites, similar to the NEON tower capability ([issue #1487](https://github.com/ESCOMP/CTSM/issues/1487)). Initial conditions are not provided for PLUMBER sites.  
* New CLM\_CMIP\_ERA flag in env\_run.xml. Valid options are cmip7 and cmip6. Defaults to cmip7 except in compsets containing SSP for which it defaults to cmip6 because there are no future-period datasets yet available for CMIP7.  
* Automatic, more flexible use of anomaly forcings for CMIP6 ISSP cases, which also use the cmip6 CLM\_CMIP\_ERA flag: [Documentation](https://escomp.github.io/CTSM/users_guide/running-special-cases/Running-with-anomaly-forcing.html)

* Unsupported script that checks for spinup equilibrium in `tools/contrib/` for spectral element grids ([PR #2991](https://github.com/ESCOMP/ctsm/pull/2991)).  
* New paramfile tools that allow users to query and modify CLM parameter files ([documentation](https://escomp.github.io/CTSM/users_guide/using-clm-tools/paramfile-tools.html))  
* Optional time-evolving \`leafcn\_target\`. More under “Additional detail” below.  
* New vertical movement scheme for soil nitrate, which is off by default (PR [#2992](https://github.com/ESCOMP/CTSM/pull/2992)).  
* Documentation improvements and new URL: https://escomp.github.io/CTSM/index.html.  
* FATES:  
  * Grazing ([sci.1.81.0\_api.37.1.0](https://github.com/NGEET/fates/releases/tag/sci.1.81.0_api.37.1.0)).  
  * Johnson and Berry 2021 electron transport model ([sci.1.85.0\_api.40.0.0](https://github.com/NGEET/fates/releases/tag/sci.1.85.0_api.40.0.0)).  
  * Managed Fire ([sci.1.87.0\_api.41.0.0](https://github.com/NGEET/fates/releases/tag/sci.1.87.0_api.41.0.0)).

## Answer changes

Changes to defaults for \`clm6\` physics:

* New CMIP7 surface and landuse timeseries datasets (see in Additional Details below).  
* New namelist variables \`snow\_thermal\_cond\_glc\_method\` and \`snow\_thermal\_cond\_lake\_method\` ([PR #3072](https://github.com/ESCOMP/CTSM/pull/3072)). Snow thermal conductivity uses Jordan1991 over glaciers to reduce Greenland melt rates by default and Sturm over land and lake land units.  
* Bytnerowicz is now the default nfix\_method for clm6 (https://github.com/ESCOMP/ctsm/pull/2972) which revises the temperature function for nitrogen fixation, replacing the Houlton *et al.* function.  
* Updates to MEGAN for BVOCs (https://github.com/ESCOMP/CTSM/pull/3065 https://github.com/ESCOMP/CTSM/pull/3309). Removes dependence on soil moisture from clm6 physics.  
* New model parameter values that were calibrated to improve carbon cycle representation with CRUJRA.  
* New model parameter values that were calibrated to improve the fire model. Now using li2024 fire code.  
* New initial conditions files for f09 ("1-degree" 1850, 2000), f19 (“2-degree” 1850), and ne30 (1850, 1979, 2000\) resolutions.  
* Change default for glcmec\_downscale\_longwave to FALSE for clm6 physics as turning off the LW downscaling improves the melt and runoff biases.  
* See “Changes to FATES and the FATES parameter file” below.  
* Namelist defaults change so that  
  * use\_c13/use\_c14 are on only for HistClm60Bgc compsets with CRUJRA2024 or CAM7 forcing; examples of when use\_c13/use\_c14 are now off include SSP and single-point compsets, as well as cases using older forcings, such as CAM6, GSWP3v1, Qian, and CRUv7  
  * when use\_c13 or use\_c14 is on, turn on the corresponding time series file  
  * irrigation is on for transient cases (1850-2000, 1850-2100, but not for clm4\_5).

Changes for all physics versions:

* Parameters updated: Added MIMICS parameter \`mimics\_fi\` (fraction of litter inputs that bypass litter pools, directly contributing to SOM) and updated other MIMICS parameters (https://github.com/ESCOMP/CTSM/pull/2365) to remove NPP control on turnover, fix density dependent control on turnover, add litterfall fluxes that bypass litter pools and contribute directly to soil organic matter.  
* FATES parameter file updated: ([PR \#2965](https://github.com/ESCOMP/CTSM/pull/2965), [PR \#2904](https://github.com/ESCOMP/CTSM/pull/2904), [PR \#1344](https://github.com/NGEET/fates/pull/1344), [PR \#3087](https://github.com/ESCOMP/CTSM/pull/3087)). See “FATES parameter file” section below for details.  
* New surface datasets and landuse timeseries files (see “surface datasets” section below).

## Heads up

* History tapes now split into two files from hX to hXi and hXa, where X is the tape number (e.g. h0i/h0a) and where "i" stands for history file containing instantaneous fields, while "a" stands for history file containing non-instantaneous fields. Details in the “history files” section below and in the PRs \[\#2445\](https://github.com/ESCOMP/ctsm/pull/2445) \[\#117\](https://github.com/ESCOMP/MOSART/pull/117) \[\#61\](https://github.com/ESCOMP/RTM/pull/61) and the corresponding issues.  
* Adding time to 1d weighting fields in transient simulations PR \\[\#3328](https://github.com/ESCOMP/CTSM/pull/3328)  
* Regarding CMIP7 vs. CMIP6 inputs:  
  * We supply only CMIP7 C13/C14 isotope datasets, so these get used regardless of CLM\_CMIP\_ERA setting.  
  * We supply only CMIP7 population density with clm6 physics in non-SSP cases, because the fire model is calibrated to that; conversely, we supply only CMIP6 population density for pre-clm6 physics and for SSP cases.  
  * We supply only CESM2 nitrogen deposition (ndep), so this gets used regardless of CLM\_CMIP\_ERA setting.  
  * For DATM we supply only CMIP6 aerosols.  
  * For DATM we supply only CMIP6 CO2.  
* Issue with DOUT\_S\_SAVE\_INTERIM\_REST [https://github.com/ESCOMP/CTSM/issues/3351](https://github.com/ESCOMP/CTSM/issues/3351) was fixed.  
* As of ctsm5.3.040, the new ctsm\_pylib conda environment is incompatible with our tools from before ctsm5.3.040 and vice versa. More under “Additional detail” below.

# Additional detail

## Changes related to history files

(Note 1: The same information in this section applies to MOSART and RTM.  
Note 2: The gist of the information in this section also appears in the [CTSM User’s Guide](https://escomp.github.io/CTSM/users_guide/setting-up-and-running-a-case/customizing-the-clm-namelist.html#various-ways-to-change-history-output-averaging-flags)).

Following ctsm5.3.018 "Change history time to be the middle of the time bounds" and keeping CLM history consistent with CAM history, the CTSM5.4 change intends to prevent confusion associated with the time corresponding to instantaneous history fields by putting them on separate files than non-instantaneous fields.

The now separate instantaneous history files represent the exact time step when they were written and do not include a time\_bounds variable. Conversely, non-instantaneous history files represent the period of their time\_bounds variable. As a result, time data on non-instantaneous history files are now read correctly during post processing (e.g. by xarray). Special handling may still be needed for instantaneous history files, whose timestamps represent the date and time at the END of the history timestep. So, e.g., an instantaneous variable saved at the end of year 2023 will get the timestamp 2024-01-01 00:00:00.

Users will now see:

1\) Two history files per clm, mosart, and rtm history tape:  
 tape h0 becomes h0a and h0i  
 tape h1 becomes h1a and h1i  
 ...  
 tape hX becomes hXa and hXi

2\) Two history-restart files per history restart tape:  
 rh0 becomes rh0a and rh0i  
 rh1 becomes rh1a and rh1i  
 ...  
 rhX becomes rhXa and rhXi

The CLM handles empty history (and corresponding history-restart) files by not generating them, while rtm and mosart give an error. Instead of refactoring rtm and mosart to behave like the clm (considered out of scope), we have introduced one active instantaneous field in mosart and one in rtm to bypass the "empty file" error.

## New surface datasets and landuse timeseries files (\[PR \#3482\](https://github.com/ESCOMP/CTSM/pull/3482))

* Transient landuse timeseries files going back to 1700 made for f09 and 360x720 grids.  
* New resolutions now supported: ne3np4.pg3, mpasa30, ne0np4.NATL.ne30x8 (\[PR \#3482\](https://github.com/ESCOMP/CTSM/pull/3482))  
* Updates to input datasets (also referred to as raw datasets):  
  * PFT/LAI/soil-color raw datasets; now from the CMIP7 timeseries that ends in 2023 (Issues [\#2851](https://github.com/ESCOMP/CTSM/issues/2851)).  
  * Two fire datasets: crop fire peak month and population density (Issue   
  * \[\#2701\]([https://github.com/ESCOMP/CTSM/issues/2701](https://github.com/ESCOMP/CTSM/issues/2701)) \[\#3302\](https://github.com/ESCOMP/CTSM/issues/3302)).  
  * Transient (historical) urban datasets are now based on CMIP7 urban data, partitioned into TBD, HD, and MD classes in proportion to GaoOneill present day classification.

## Changes to FATES and the FATES parameter file

* See [HLM-FATES compatibility table](https://fates-users-guide.readthedocs.io/en/latest/user/release-tags-compat-table.html) in the FATES user’s guide for all FATES tags associated with CTSM tag updates  
* FATES answer changing updates  
  * The default hydro solver is updated to 2D Picard from 1D Taylor ([ctsm5.3.027](https://github.com/ESCOMP/CTSM/releases/tag/ctsm5.3.027))  
  * Simplified leaf sun-shade fraction for two-stream radiation ([sci.1.83.0\_api.39.0.0](https://github.com/NGEET/fates/releases/tag/sci.1.83.0_api.39.0.0))  
  * Default maximum canopy layer updated from 2 to 3 ([sci.1.87.1\_api.41.0.0](https://github.com/NGEET/fates/releases/tag/sci.1.87.1_api.41.0.0))  
  * Various bug fixes (see compatibility table)  
* FATES Parameter File Updates   
  * ctsm5.3.025 (API 37\)  
    * Adds pft-dependent btran model switches  
    * Adds parameters for land use grazing  
    * Updates the FATES z0mr turbulence parameters for consistency with CLM  
  * ctsm5.3.027 (API 38\)  
    * Migrates a number of global parameter file variables to the namelist  
    * Adds \`fates\_leaf\_fnps\` parameter for the electron transport model  
    * \`fates\_leaf\_theta\_cj\_c3\` and \`fates\_leaf\_theta\_cj\_c4\` depricated  
  * ctsm5.3.045 (API 40\)  
    * Changes to the default competitive exclusion parameter from probabilistic to rank-ordered sorting of cohorts by default  
    * Sets the logging default to clear cut  
    * Refactors the pft-specific phenology habit selection into a single parameter  
  * ctsm5.3.070 (API 41\)  
    * Add parameters for the managed fire feature addition  
    * Corrects the fates landuse crop pft to c3 cool grass

## New ctsm\_pylib conda environment

If you have a ctsm\_pylib conda environment installed from before ctsm5.3.040, you may want to keep that under a different name. We suggest the following command for doing this in a local copy of ctsm5.3.040 or later:

```shell
./py_env_create -r ctsm_pylib_old
```

This first renames your existing ctsm\_pylib to ctsm\_pylib\_old and then installs the Python 3.13.2 version as ctsm\_pylib. If you are unsure whether you already have ctsm\_pylib installed, use the same command regardless, as it will skip the renaming step if necessary.

Information about additional py\_env\_create options — including how to install a fresh copy of the old conda environment — is available as follows:

```shell
./py_env_create --help
```

## Potentially time-evolving \`leafcn\_target\` replaces time-constant \`leafcn\`

The former is calculated as a function of the latter and can be time-evolving depending on new paramfile parameter \`leafcn\_co2\_slope\` https://github.com/ESCOMP/ctsm/pull/1654. The time-evolving effect defaults to off with \`leafcn\_co2\_slope\` \= 0 on the parameter file.

# Simulations supporting this release by providing initial conditions

* f19 \`Clm60BgcCruJra\` 16pft: https://github.com/NCAR/LMWG_dev/issues/125  
* f09 with \`Clm60BgcCropCruJra\`: https://github.com/NCAR/LMWG_dev/issues/124  
* ne30 with \`Clm60BgcCropCruJra\`: https://github.com/NCAR/LMWG_dev/issues/123 (123\_HIST\_popDens)  
* ne30 SP https://github.com/NCAR/LMWG_dev/issues/126  
* f09 SP https://github.com/NCAR/LMWG_dev/issues/127

