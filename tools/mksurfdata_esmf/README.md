# Instructions for Using mksurfdata_esmf to Create Surface Datasets
#### $CTSMROOT/tools/mksurfdata_esmf/README.md

## Table of contents
1. Purpose NOW IN THE USER'S GUIDE https://escomp.github.io/CTSM/users_guide/using-clm-tools/creating-surface-datasets.html#mksurfdata-esmf-purpose
2. Build Requirements NOW IN THE USER'S GUIDE https://escomp.github.io/CTSM/users_guide/using-clm-tools/creating-surface-datasets.html#build-requirements
3. [Building the executable](#building-the-executable)
4. [Running a Single Submission](#running-for-a-single-submission)
5. [Running for Multiple Datasets](#running-for-the-generation-of-multiple-datasets)
6. [Notes](#notes)

<!-- ================== -->
### Building the executable
<!-- ================== -->

 Before starting, be sure that you have run

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 ./bin/git-fleximod update  # Assuming at the top level of the CTSM/CESM checkout
```

This will bring in CIME and ccs_config which are required for building.

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
setenv DEBUG TRUE  # only if debugging and your shell is tcsh (in bash use: export DEBUG=TRUE)
 ./gen_mksurfdata_build         # For machines with a cime build
```

 Note: The pio_iotype value gets set and written to a simple .txt file
 by this build script. The value depends on your machine. If not running
 on derecho, casper, or izumi, you may need to update this, though
 a default value does get set for other machines.

<!-- ========================= -->
## Running for a single submission
<!-- ========================= -->

### Setup ctsm_pylib
 Work in the ctsm_pylib environment, which requires the following steps when
 running on Derecho. On other machines it will be similar but might be different
 in order to get conda in your path and activate the ctsm_pylib environment.

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 module load conda
 cd ../..  # or ../../../.. for a CESM checkout)
 ./py_env_create    # Assuming at the top level of the CTSM/CESM checkout
 conda activate ctsm_pylib
```

to generate your target namelist:

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 ./gen_mksurfdata_namelist --help
```

for example try --res 1.9x2.5 --start-year 1850 --end-year 1850:

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 ./gen_mksurfdata_namelist --res <resolution> --start-year <year1> --end-year <year2>
```

> [!TIP]
> **IF FILES ARE MISSING FROM** /inputdata, a target namelist will be generated
> but with a generic name and with warning to run `./download_input_data` next.
> **IF A SMALLER SET OF FILES IS STILL MISSING AFTER RUNNING** `./download_input_data`
> and rerunning `./gen_mksurfdata_namelist`, then rerun
> `./gen_mksurfdata_namelist with your options needed.
> and rerun `./download_input_data` until
> `./gen_mksurfdata_namelist` finds all files.

 Example, to generate your target jobscript (again use --help for instructions):

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 ./gen_mksurfdata_jobscript_single --number-of-nodes 2 --tasks-per-node 128 --namelist-file target.namelist
 qsub mksurfdata_jobscript_single.sh
```

 Read note about regional grids at the end.

<!-- ========================================= -->
## Running for the generation of multiple datasets
<!-- ========================================= -->
 Work in the ctsm_pylib environment, as explained in earlier section.
 gen_mksurfdata_jobscript_multi runs `./gen_mksurfdata_namelist` for you

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 ./gen_mksurfdata_jobscript_multi --number-of-nodes 2 --scenario global-present
 qsub mksurfdata_jobscript_multi.sh
```

 If you are looking to generate all (or a large number of) the datasets or the
 single-point (1x1) datasets, you are best off using the Makefile. For example

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
 make all  # ...or
 make all-subset
```

 As of 2024/9/12 one needs to generate NEON and PLUMBER2 fsurdat files by
 running ./neon_surf_wrapper and ./plumber2_surf_wrapper manually in the
 /tools/site_and_regional directory.

<!-- = -->
## NOTES
<!-- = -->

# Guidelines for input datasets to mksurfdata_esmf

> [!TIP]
> ALL raw datasets \*.nc **FILES MUST NOT BE NetCDF4**.

Example to convert to CDF5

``` shell
nccopy -k cdf5 oldfile newfile
```

> [!TIP]
> The LAI raw dataset \*.nc **FILE MUST HAVE** an "unlimited" time dimension

Example to change time to unlimted dimension using the NCO operator ncks.

``` shell
ncks --mk_rec_dmn time file_with_time_equals_12.nc -o file_with_time_unlimited.nc
```

### IMPORTANT THERE HAVE BEEN PROBLEMS with REGIONAL grids!!

> [!CAUTION]
> See 
>
> https://github.com/ESCOMP/CTSM/issues/2430

In general we recommend using subset_data and/or fsurdat_modifier
for regional grids.
