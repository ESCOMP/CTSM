# New single-point documentation

## 1. Subset the data
For single-point cases, you need to subset a surface dataset and (optionally) DATM data. The Python script to subset this data can be found in the CTSM repository at `tools/site_and_regional/subset_data`.

Note that you will need to have a python environment set up that includes the packages `scipy`, `xarray`, and `numpy`. If you have conda or miniconda installed, you can create a conda environment for this and other CTSM python tools using the script `py_env_create` at the top level of your CTSM checkout.

To subset surface data and climate forcings (DATM) for a single point, use the command:

```shell
tools/site_and_regional/subset_data point --lat $my_lat --lon $my_lon --site $my_site_name --create-surface --create-datm --datm-syr $my_start_year --datm-eyr $my_end_year --create-user-mods --outdir $my_output_dir
```

- `$my_lat`: latitude of point, *must be between -90 and 90 degrees*.
- `$my_lon`: longitude of point, *must be between 0 and 360 degrees*
- `$my_site_name`: name of site, *used for file naming*
- `$my_start_year`: start year for DATM data to subset, *default between 1901 and 2014*
- `$my_end_year`: end year for DATM data to subset, *default between 1901 and 2014*
- `$my_output_dir`: output directory to place the subset data and user_mods directory. This should be something specific to _just_ your data for `$my_site_name`.

```shell
my_lat=40
my_lon=55
my_site_name=boulder_test_202503
my_start_year=2010
my_end_year=2011
my_output_dir=/glade/work/samrabin/subset_data_outputs
tools/site_and_regional/subset_data point --lat $my_lat --lon $my_lon --site $my_site_name --create-surface --create-datm --datm-syr $my_start_year --datm-eyr $my_end_year --create-user-mods --outdir $my_output_dir --overwrite
```

You can also have the script subset landuse data. See the help (`tools/site_and_regional/subset_data --help`) for all argument options.

**_Note that this script defaults to subsetting specific surface, domain, and landuse files and GSWP3 DATM data, and can currently only be run as-is on Derecho. To update the files and file locations, you will need to modify the files and/or directories in the [`default_data_2000.cfg`](https://github.com/ESCOMP/CTSM/blob/master/tools/site_and_regional/default_data_2000.cfg) file._**

The `--create-user-mods` command tells the script to set up a user mods directory in your specified `$my_output_dir` and to specify the required `PTS_LAT` and `PTS_LON` parameters. You can then use this user mods directory to set up your CTSM case, as described below.

## 2. Create the case

You can use the user mods directory set up in the previous subset data step to tell CIME/CTSM where your subset files are located.

```shell
cime/scripts/create_newcase --case $my_case_name --res CLM_USRDAT --compset $compset --run-unsupported --user-mods-dirs $my_output_dir/user_mods
```

```shell
cime/scripts/create_newcase --case boulder_test_202503 --res CLM_USRDAT --compset I1850Clm60Bgc --run-unsupported --user-mods-dirs /glade/work/samrabin/subset_data_outputs/boulder_test_202503/user_mods
```

- `$my_case_name`: the path of the case directory you want to create
- `$compset`: the compset you would like to use (see the [above-mentioned tutorial](https://github.com/NCAR/CTSM-Tutorial-2022/blob/main/notebooks/Day2a_GenericSinglePoint.ipynb) for an example)
Note the use of `$my_output_dir/user_mods` which is the `user_mods/` directory that the subset data script set up within your specified `$my_output_dir`.

Following this, you should be able to update any other case-specific parameters you want to change (e.g., `STOP_N`, `STOP_OPTION`, etc.). You should set the values `DATM_YR_ALIGN`, `DATM_YR_START`, and `DATM_YR_END` to match your `$my_end_year` and `$my_start_year` values from above: **Shouldn't this just be in `shell_commands`?**

```shell
./xmlchange DATM_YR_ALIGN=$my_start_year
./xmlchange DATM_YR_START=$my_start_year
./xmlchange DATM_YR_END=$my_end_year
```

```shell
./xmlchange DATM_YR_ALIGN=2010
./xmlchange DATM_YR_START=2010
./xmlchange DATM_YR_END=2011
```

Note that `./case.setup` on Derecho will automatically set queue to `develop` and walltime to one hour. You might need a longer walltime, but the maximum walltime for `develop` is one hour. To change it to two hours on Derecho:
```shell
./xmlchange --subgroup case.run JOB_QUEUE=main,JOB_WALLCLOCK_TIME=2:00:00
```
