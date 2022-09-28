# Software Design Documentation

<span style="font-size:larger;"><b>Project Name</b></span>
<!-- Note you can add this to an issue or pull request in github. That is the normal usage of this template. -->

**Date:**: 09/27/22

**Written By**: Adrianna Foster (@adrifoster)

## Introduction
---------------------------------------

As laid out in CTSM Issue [#270](https://github.com/ESCOMP/CTSM/issues/270) we want to downscale input ozone partial pressure (mol/mol) temporally if we receive a multi-day average input ozone (usually in DATM mode).

We have the infrastructure laid out for CAM and DATM to inform CTSM which type of input we are receiving (i.e., 'multiday_average' or 'subdaily'). We should only apply this downscaling when receiving multiday average ozone.

We have a gridded file provided by Louisa Emmons for a diurnal ozone variation factor. The data is additionally dimensioned by 'secs' (seconds of day), and depending on the model time of day, can be used as a multiplicative factor for converting multi-day average ozone to sub-daily ozone partial pressure.

## Solutions
---------------------------------------

Read in the diurnal anomaly file as a streams file, being careful to not assume too much so that other files may be provided in the future. Create a diurnal ozone type that has an `Interp` method to downscale and interpolate multi-day average ozone data to sub-daily based on a factor attribute. Inside the existing `CalcOzoneUptake` routine, if we are using multi-day average ozone, interpolate/downscale the `forc_o3` array using the `Interp` subroutine.

**Some alternate ideas**: 
1. have the existing `ozone_type` type have this diurnal ozone anomaly type as an attribute that gets initiated only if we are reading in multi-day average ozone. 
2. have the diurnal ozone type be fed into relevant `ozone_type` methods as arguments. We would only want to do this if other modules are going to use the diurnal ozone type. This has a possibility, but no clear plans are on the horizon for it.

**Consider going with the first solution for now**

## Design Considerations
---------------------------------------

### Assumptions and Dependencies:

Right now we are assuming that the units of the input diurnal file are going to be in seconds, and that the file covers the whole day. 

There is no need to assume that the dimension will be 1-24 (i.e. on the hour) or that the dimensions be equal intervals.


## Design and Architecture
---------------------------------------

### System diagram or flowchart

#### Diurnal Factor Streams File

We need to read in the diurnal ozone factor file provided by Louisa Emmons as a streams file. Some modifications are required on the file so that it can be accurately read and parsed by ESMF and to facilitate our chosen implementation strategy:

1. Add a new 'time' variable (even though the data is not dimensioned by time, except for seconds of day) so that ESMF can parse it. The date is arbitrarily set to 2000-01-01, and the dimensions are set to 'UNLIMITED'.

2. Shift the 'secs' array up 1800 seconds so that it provides the midpoint of the averaged time interval, rather than the beginning. This conversion is done to facilitate interpolation between the time-of-day dimensions.

These file modifications can be seen in the Jupyter notebook /glade/u/home/afoster/Diurnal_ozone.ipynb.

The file is read in using a new `src/share_esmf/diurnalOzoneStreamMod` module, with associated updates to `bld/CLMBuildNamelist.pm`, `bld/namelist_files/namelist_defaults_ctsm.xml` and `bld/namelist_files/namelist_definition_csvm.xml`

Right now the new variable `use_do3_streams` is set in the `user_nl_clm` file. And the other namelist variables (`stream_fldfilename_do3`, `stream_meshfile_do3`, and `do3_mapalgo`) are inside the `lnd_in` file under a `do3_streams` namelist. 

**I'm not sure setting up the `use_do3_streams` namelist variable is the best way to go about this, since we essentially want this feature on whenever the ozone frequency is 'multiday_average'**

Right now we are hard-coding the `stream_var_name`, and `stream_lev_dimname`. I think this is okay. Otherwise, we would want them to be namelist variables that get read in.

Because the diurnal anomaly file is not dimensioned by time other than seconds of day, we just need to read it in and advance one time. This is all currently done in the `diurnalOzoneStreamMod`'s `read_O3_stream` subroutine. We set the date to an arbitrary date (the same as on the file, I'm not sure if this is necessary?).

We also need to grab the `secs` array to do interpolation/downscaling with. 

Both are read from the file and then an instance of `diurnal_ozone_anom_type` is initialized and relevant arrays are set to the ozone factor and seconds data.

### Diurnal Ozone Anomaly Type

The module `DiurnalOzoneType` sets up a `diurnal_ozone_anom_type` which has only a few attributes and methods:

**Attributes**:
1. `ntimes` - size of time/seconds-of-day dimension (private)
2. `o3_anomaly_grc` - o3 anomaly data, 2d, gridcells x ntimes
3. `time_arr` - time/seconds of day array, 1d, ntimes

**Methods**:
1. `Init` - Initializes the anomaly data structures, calls `InitAllocate`
2. `InitAllocate` - allocates arrays and sets them to nan
3. `Interp` - Interpolates/downscales an input multi-day average ozone (`forc_o3`) using the o3 anomaly data and outputs a downscaled `forc_o3_down` array. See below for `Interp` algorithm/pseudo code

### Implementation within OzoneMod

We add an instance of the `diurnal_ozone_anom_type` as an attribute of the `ozone_base_type`. 

Within the `Init` method of the `ozone_base_type`, we read the `drv_flds_in` file to get the `atom_ozone_frequency_val`. If this value is `atm_ozone_frequency_multiday_average` then we initialize and read in the o3 anomaly data by calling the `read_O3_stream` described above. *Do we need an else decision here?*

We also add an integer `atm_ozone_freq` as a new attribute, and set it in `Init`.

Within the `CalcOzoneUptake` method, we check for the `atm_ozone_freq` flag and if it is `atm_ozone_frequency_multiday_average` we call the `Interp` method.

### Algorithm and Pseudo code for Interp



## Rollout Plan
---------------------------------------
Define the roll-out phases and tests you plan to do


## Review Sign-off
---------------------------------------
* Reviewer(s):

*Sign-off Completed on YYYY-MM-DD*
