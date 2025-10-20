# What's new in CTSM 5.4 (tag `ctsm5.4.0xx`)

## Purpose and description of changes since CTSM 5.3 (tag `ctsm5.3.021`)

REMOVE THESE NOTES WHEN DONE AND READY TO RELEASE
- As of 2025/7/29 slevis has gone through the ChangeLog from cctsm5.3.022 to ctsm5.3.065
- Ask for reviewers to browse/read this doc for accuracy, omissions, and redundancies
- For omissions, request contributions from the relevant developers

### New features

* Can now choose the CRUJRA2024 atmospheric driver data with clm6 and clm5 [PR \#2956](https://github.com/ESCOMP/ctsm/pull/2956))
* Can now run PLUMBER towers, similar to the NEON tower capability ([issue \#1487](https://github.com/ESCOMP/CTSM/issues/1487))
* New namelist variables `snow_thermal_cond_glc_method` and `snow_thermal_cond_lake_method` ([PR \#3072](https://https://github.com/ESCOMP/CTSM/pull/3072))
* Potentially time-evolving `leafcn_target` replaces time-constant `leafcn`: the former is calculated as a function of the latter and can be time-evolving depending on new paramfile parameter `leafcn_co2_slope` ([PR \#1654](https://github.com/ESCOMP/ctsm/pull/1654)
* Easier, more flexible use of anomaly forcings for ISSP cases ([PR \#2686](https://github.com/ESCOMP/CTSM/pull/2686) [PR \#3212](https://github.com/ESCOMP/CTSM/pull/3212))
* New equilibrium script in /tools/contrib for spectral element grids ([PR \#2991](https://github.com/ESCOMP/ctsm/pull/2991))

### Answer changes

Changes to defaults for `clm6_0` physics:

* Bytnerowicz is now the default nfix_method for clm6 ([PR \#2972](https://github.com/ESCOMP/ctsm/pull/2972))
* New initial conditions files for f09 ("1-degree" 1850, 2000), f19 (“2-degree” 1850), and ne30 (1850, 1979, 2000) resolutions?
* Updates to MEGAN for BVOCs ([PR \#3065](https://github.com/ESCOMP/CTSM/pull/3065) [PR \#3309](https://github.com/ESCOMP/CTSM/pull/3309))

Changes for all physics versions:

* Parameters updated: Added MIMICS parameter `mimics_fi`. ([PR \#2365](https://github.com/ESCOMP/CTSM/pull/2365))
* FATES parameter file updated: ([PR \#2965](https://github.com/ESCOMP/CTSM/pull/2965) [PR \#2904](https://github.com/ESCOMP/CTSM/pull/2904) [PR \#1344](https://https://github.com/NGEET/fates/pull/1344) [PR \#3087](https://github.com/ESCOMP/CTSM/pull/3087))
* New surface datasets and landuse timeseries files (see section below).

### Heads up

* History tapes now split into two files from hX to hXi and hXa, where X is the tape number (e.g. h0i/h0a) and where "i" stands for history file containing instantaneous fields, while "a" stands for history file containing non-instantaneous fields. Details below and in the PRs [\#2445](https://github.com/ESCOMP/ctsm/pull/2445) [\#117](https://github.com/ESCOMP/MOSART/pull/117) [\#61](https://github.com/ESCOMP/RTM/pull/61) and the correspondng issues.
* As of ctsm5.3.040, the new ctsm_pylib conda environment is incompatible with our tools from before ctsm5.3.040 and vice versa. If you have a ctsm_pylib conda environment installed from before ctsm5.3.040, keep that under a different name. We suggest the following command for doing this in a local copy of ctsm5.3.040 or later:

./py_env_create -r ctsm_pylib_old

This first renames your existing ctsm_pylib to ctsm_pylib_old and then installs the Python 3.13.2 version as ctsm_pylib. If you are unsure whether you already have ctsm_pylib installed, use the same command regardless, as it will skip the rename step if necessary.

Information about additional py_env_create options — including how to install a fresh copy of the old conda environment — is available as follows:

./py_env_create --help

##

## Additional detail

### Changes related to history files

(Note that the same information in this section applies to MOSART and RTM.)

Following ctsm5.3.018 "Change history time to be the middle of the time bounds" and to keep CLM history consistent with CAM history, the CTSM5.4 change intends to prevent confusion associated with the time corresponding to instantaneous history fields by putting them on separate files than non-instantaneous fields. The result is

1) two history files per clm, mosart, and rtm history tape:
 tape h0 becomes h0a and h0i
 tape h1 becomes h1a and h1i
 ...
 tape hX becomes hXa and hXi

2) two history restart files per history restart tape:
 rh0 becomes rh0a and rh0i
 rh1 becomes rh1a and rh1i
 ...
 rhX becomes rhXa and rhXi

The CLM handles empty history (and corresponding history restart) files by not generating them, while rtm and mosart give an error. Instead of refactoring rtm and mosart to behave like the clm (considered out of scope), we have introduced one active instantaneous field in mosart and one in rtm to bypass the "empty file" error.

### New surface datasets and landuse timeseries files ([PR \#xxxx](https://github.com/ESCOMP/CTSM/pull/xxxx))

* Transient landuse timeseries files going back to 1700 made for f09.
* Surface datasets now provided for the `xxxx` grid. ([PR \#xxxx](https://github.com/ESCOMP/CTSM/pull/xxxx), [issue \#xxxx](https://github.com/ESCOMP/CTSM/issues/xxxx))
* Updates to input datasets:
  * PFT/LAI/soil-color raw datasets; now from the xxxx timeseries that ends in 2023. (Issues [\#xxxx](https://github.com/ESCOMP/CTSM/issues/xxxx)
  * Two fire datasets: crop fire peak month and population density. (Issue [\#xxxx](https://github.com/ESCOMP/CTSM/issues/xxxx))

### Changes to FATES parameter file

  * 

### Changes to rpointer files


## Simulations supporting this release

- f19 `Clm60Bgc` 16pft: [https://github.com/NCAR/LMWG\_dev/issues/xxx](https://github.com/NCAR/LMWG_dev/issues/xxx)
- f09 with `Clm60BgcCrop`: [https://github.com/NCAR/LMWG\_dev/issues/xxx](https://github.com/NCAR/LMWG_dev/issues/xxx)
- ne30 with `Clm60BgcCrop`: [https://github.com/NCAR/LMWG\_dev/issues/xxx](https://github.com/NCAR/LMWG_dev/issues/xxx)

