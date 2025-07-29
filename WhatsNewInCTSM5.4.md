# What's new in CTSM 5.4 (tag `ctsm5.4.0xx`)

## Purpose and description of changes since CTSM 5.3 (tag `ctsm5.3.021`)

### New features

* 

### Answer changes

Changes to defaults for `clm6_0` physics:

* New initial conditions files for f09 ("1-degree" 1850, 2000), f19 (“2-degree” 1850), and ne30 (1850, 1979, 2000) resolutions.
* Updates to MEGAN for BVOCs. ([PR \#xxxx](https://github.com/ESCOMP/CTSM/pull/xxxx))

Changes for all physics versions:

* Parameters updated: PPE-based modifications were made to the parameters `xxxx`. ([PR \#xxxx](https://github.com/ESCOMP/CTSM/pull/xxxx))
* FATES parameter file updated?
* New surface datasets and landuse timeseries files (see section below).

### Heads up

* History tapes now split into two files from hX to hXi and hXa, where X is the tape number (e.g. h0i/h0a) and where "i" stands for history file containing instantaneous fields, while "a" stands for history file containing non-instantaneous fields. Details below and in the PRs [\#2445](https://github.com/ESCOMP/ctsm/pull/2445) [\#117](https://github.com/ESCOMP/MOSART/pull/117) [\#61](https://github.com/ESCOMP/RTM/pull/61) and the correspondng issues.

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

