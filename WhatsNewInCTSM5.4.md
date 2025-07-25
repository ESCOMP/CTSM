# What's new in CTSM 5.3 (tag `ctsm5.3.021`)

## Purpose and description of changes since CTSM 5.2 (tag `ctsm5.2.005`)

### New features

* `manage_externals` replaced by [`git-fleximod`](https://github.com/ESMCI/git-fleximod/blob/main/README.md). ([PR \#2559](https://github.com/ESCOMP/CTSM/pull/2559))
* No longer runs the 0th time step in first segment of startup and hybrid runs; branch and continue runs never had this 0th time step. ([PR \#2084](https://github.com/ESCOMP/CTSM/pull/2084))
* New CN Matrix method speeds up spinup for the BGC model. ([PR \#640](https://github.com/ESCOMP/CTSM/pull/640); [Liao et al. 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023MS003625)).
* New `Leung_2023` dust emissions. ([PR \#1897](https://github.com/ESCOMP/CTSM/pull/1897); [Leung et al 2023](https://doi.org/10.5194/acp-23-6487-2023), [Leung et al. 2024](https://doi.org/10.5194/acp-24-2287-2024))
* Explicit air conditioning for the urban model. ([PR \#2275](https://github.com/ESCOMP/CTSM/pull/2275); [Li et al. 2024](https://agupubs.onlinelibrary.wiley.com/share/NY4AYPREB8Y8BUDP7DXD?target=10.1029/2023MS004107))
* FATES compsets can now be run with transient land use; off by default. ([PR \#2507](https://github.com/ESCOMP/CTSM/pull/2507))
* Ability to handle CAM7 in `LND_TUNING_MODE`. ([PR \#2632](https://github.com/ESCOMP/CTSM/pull/2632))

### Answer changes

Changes to defaults for `clm6_0` physics:

* Urban explicit A/C turned on (links above).
* Snow thermal conductivity method is now `Sturm1997`. ([PR \#2348](https://github.com/ESCOMP/CTSM/pull/2348); see also [discussion \#1960](https://github.com/ESCOMP/CTSM/discussions/1960))
* New initial conditions files for f09 ("1-degree" 1850, 2000), f19 (“2-degree” 1850), and ne30 (1850, 1979, 2000) resolutions.
* New crop calendars. ([PR \#2664](https://github.com/ESCOMP/CTSM/pull/2664); informed by [Rabin et al., 2023](https://gmd.copernicus.org/articles/16/7253/2023/gmd-16-7253-2023.html))
* Dust emissions method is now `Leung_2023` (links above).
* Excess ice is turned on. ([PR \#1787](https://github.com/ESCOMP/CTSM/pull/1787))
* Updates to MEGAN for BVOCs. ([PR \#2588](https://github.com/ESCOMP/CTSM/pull/2588))
* New BGC fire method `li2024crujra`: Avoid crop fires during growing season; allow lightning ignitions in tropical closed forests; add effect of landscape fragmentation on ignitions and duration; recalibrate against GFED5 burned area and with CRU-JRA climate. ([PR \#2684](https://github.com/ESCOMP/CTSM/pull/2684), [PR \#2711](https://github.com/ESCOMP/CTSM/pull/2711), [PR \#2715](https://github.com/ESCOMP/CTSM/issues/2715))

Changes for all physics versions:

* Parameters updated for CRU-JRA forcing: PPE-based modifications were made to the parameters `leafcn`, `slatop`, `froot_leaf`, `medlynslope`, and `kmax`. Lowers LAI and biomass in boreal and tropical forests, without reducing latent heat in the tropics. Affected PFTs: NET temperate, NET boreal, BET tropical, BDS boreal, C3 arctic grass. ([PR \#2500](https://github.com/ESCOMP/CTSM/pull/2500))
* FATES parameter file updated (see section below).
* Pass active glacier (CISM) runoff directly to river model (MOSART) ([MOSART PR \#94](https://github.com/ESCOMP/MOSART/pull/94))
* New surface datasets and landuse timeseries files (see section below).
* CNMatrix is new default spinup method (links above).

### Heads up

* Small glacier changes mean that you can’t use a 5.3 surface dataset with pre-5.3 code and vice versa anymore. (Merged with [PR \#2500](https://github.com/ESCOMP/CTSM/pull/2500))
* Updates the definition of history variable “time” from *end* of `time_bounds` to *middle* of `time_bounds`. ([PR \#2838](https://github.com/ESCOMP/CTSM/pull/2838); see section below)
* Standardizes history variable attributes and a history dimension name. ([PR \#2052](https://github.com/ESCOMP/CTSM/pull/2052); see section below)

##

## Additional detail

### Changes related to time and history files

(Note that the same information in this section applies to MOSART and RTM.)

Startup and hybrid runs no longer run the 0th time step, consistent with the same change in CAM. (Branch and continue runs never had this 0th time step.) This means you will not get an extraneous initial history file anymore. In some circumstances this may also affect the names of history files.

In most cases, the history `time` variable is now defined as the middle of a history file’s `time_bounds` instead of the end, for consistency with the same change in CAM. The exception is if you specify `hist_avgflag_pertape = 'I'` for that file, in which case it will be treated as an “instantaneous” file. Instantaneous history files (a) have their `time` coordinate set to the end of the last timestep (as did all history files before this tag) and (b) do not include `time_bounds`.

The history dimension name `hist_interval` (of output variable `time_bounds`) is standardized to be `nbnd`. History variables `time_bounds`, `mcdate`, `mcsec`, `mdcur`, and `mscur` are standardized to include the calendar attribute.

### New surface datasets and landuse timeseries files ([PR \#2500](https://github.com/ESCOMP/CTSM/pull/2500))

* Transient landuse timeseries files going back to 1700 now possible (and made for f09).
* Fix an important bug on soil fields that was there since `ctsm5.2.0`. This has the side effect of `mksurfdata_esmf` now giving identical answers with a change in number of processors, as it should. ([Issue \#2744](https://github.com/ESCOMP/CTSM/issues/2744))
* Surface datasets now provided for the `ne0np4.POLARCAP.ne30x4` grid. ([PR \#2716](https://github.com/ESCOMP/CTSM/pull/2716), [issue \#2720](https://github.com/ESCOMP/CTSM/issues/2720))
* Surface datasets now have their version number embedded to prevent mismatch of surface dataset and CTSM version. ([Issue \#2723](https://github.com/ESCOMP/CTSM/issues/2723))
* Remove outdated hydrology `VIC` (Variable Infiltration Capacity Hydrology model) fields from surface datasets.
* Updates to input datasets:
  * PFT/LAI/soil-color raw datasets; now from the TRENDY2024 timeseries that ends in 2023. (Issues [\#2570](https://github.com/ESCOMP/CTSM/issues/2570) and [\#2452](https://github.com/ESCOMP/CTSM/issues/2452))
  * Two fire datasets: crop fire peak month and peatland fraction. (Issue [\#2618](https://github.com/ESCOMP/CTSM/issues/2618))
  * Glacier behavior dataset (related to how non-Greenland glaciers are handled). (Issue [\#423](https://github.com/ESCOMP/CTSM/issues/423))

### Changes to FATES parameter file

* [PR \#2507](https://github.com/ESCOMP/CTSM/pull/2507) ([ctsm5.2.013](https://github.com/ESCOMP/CTSM/releases/tag/ctsm5.2.013)) / [FATES PR \#1116](https://github.com/NGEET/fates/pull/1116) ([sci.1.77.1\_api.36.0.0](https://github.com/NGEET/fates/releases/tag/sci.1.77.0_api.36.0.0))
  * Adds new parameters for new land use harvest mode and land use by PFT capabilities
* [PR \#2700](https://github.com/ESCOMP/CTSM/pull/2700) ([ctsm5.3.003](https://github.com/ESCOMP/CTSM/releases/tag/ctsm5.3.003)) / [FATES PR \#1255](https://github.com/NGEET/fates/pull/1255) ([sci.1.78.3\_api.36.1.0](https://github.com/NGEET/fates/releases/tag/sci.1.78.3_api.36.1.0))
  * Adds two arctic shrub PFTs, increasing the number of default PFTs to 14 and update arctic grass parameters ([FATES PR\#1236](https://github.com/NGEET/fates/pull/1236))
  * Splits `fates_turnover_leaf` parameter into canopy and understory specific turnover rates and provides new values for both ([FATES PR\#1136](https://github.com/NGEET/fates/pull/1136))
  * Updates grass allometry parameters ([FATES PR\#1206](https://github.com/NGEET/fates/pull/1206))
  * Changes the prescribe nutrient uptake defaults from 1 to 0 for all PFTs (only relevant to ELM-FATES)

### Changes to rpointer files

The rpointer files are simple text files that CESM uses to keep track of how far simulations have progressed, pointing to the filename of the latest restart file for that component. There is one such file for each component, so for CTSM `I` cases that's `lnd`, `cpl`, and `atm` (and `rof` if it's active). Normally, when the user is just extending the length  of simulations, there’s no need to worry about these files.

However, if there was a problem when a simulation shut down, it's possible that different components will have mismatched restarts and rpointer files. In the past, this meant figuring out what restart file should be pointed to in each component rpointer file and correcting it by hand in an editor. There was only the final set of rpointer files that was kept for a case.

Now, with this update, the `lnd`, `cpl`, and `atm` rpointer files have the simulation date in the filenames, so it's easy to spot if the restarts are mismatched for one of the components. Also, since there are matching rpointer files for each time restarts are created, it's now easier to (a) make sure restarts and rpointer files are all correctly matched and (b) for a user to take a set of restarts and matching rpointer files to start up from for any part of an existing simulation. This means you don't have to hand-edit the rpointer files, making sure you don't make a mistake when you do.

Old rpointer filenames:

* `rpointer.atm`
* `rpointer.cpl`
* `rpointer.lnd`

New names:

* `rpointer.atm.YYYY-MM-DD-SSSSS`
* `rpointer.cpl.YYYY-MM-DD-SSSSS`
* `rpointer.lnd.YYYY-MM-DD-SSSSS`

Where `YYYY-MM-DD-SSSSS` is the year, month, day, and seconds into the day for the simulation timestamp. For example, `rpointer.lnd.2000-01-01-00000` for a rpointer file for starting at midnight (beginning of the day) on January 1st, 2000\.

Note that this is backwards-compatible, so for all the components you can use either the new or old format for the rpointer filenames. Thus, if you are restarting from an existing case before `ctsm5.3.016`, you CAN use the rpointer filenames that don't have the timestamps in the name.

## Simulations supporting this release

- f19 `Clm60Bgc` 16pft: [https://github.com/NCAR/LMWG\_dev/issues/70](https://github.com/NCAR/LMWG_dev/issues/70)
- f09 with `Clm60BgcCrop`: [https://github.com/NCAR/LMWG\_dev/issues/69](https://github.com/NCAR/LMWG_dev/issues/69)
- ne30 with `Clm60BgcCrop`: [https://github.com/NCAR/LMWG\_dev/issues/68](https://github.com/NCAR/LMWG_dev/issues/68)

Note: Dust emissions in CTSM 5.3 will be different from the above simulations because of a tuning update that came in after those.