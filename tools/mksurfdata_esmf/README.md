# Instructions for Using mksurfdata_esmf to Create Surface Datasets

## Table of contents
1. [Purpose](#purpose)
1. [Building](#building)
1. [Running a Single Submission](#running-for-a-single-submission)
1. [Running for Multiple Datasets](#running-for-the-generation-of-multiple-datasets)
1. [Notes](#notes)

<!-- ======= -->
## Purpose
<!-- ======= -->
This tool is intended to generate fsurdat files (surface datasets) for the
CTSM. It can generate global, regional, and single-point fsurdat files, as long
as a mesh file is available for the grid.

The subset_data tool allows users to make fsurdat files from existing fsurdat
files when a mesh file is unavailable. Generally, users should consider the
subset_data tool for generating regional and single-point fsurdat files.

<!-- ======= -->
## Building
<!-- ======= -->

<!-- ============= -->
### Build Requirements
<!-- ============= -->

mksurfdata_esmf is a distributed memory parallel program (using Message Passing
Interface -- MPI) that utilizes both ESMF (Earth System Modelling Framework)
for regridding as well as PIO (Parallel I/O) and NetCDF output. As
such, libraries must be built for the following:

1. MPI
2. NetCDF
3. PIO
4. ESMF

In addition for the build: python, bash-shell, CMake and GNU-Make are required

These libraries need to be built such that they can all work together in the
same executable. Hence, the above order may be required in building them.

CTSM submodules that are required are: cime and ccs_config. See [Building](#building-the-executable) on getting
those. A python environment that includes particular packages is also required
we demonstrate how to use the ctsm_pylib environment that we support in CTSM.

Note, PNETCDF is an optional library that can be used, but is NOT required.

#### Use cime to manage the build requirements

See [IMPORTANT NOTE](important note-only-working-on-derecho-currently)

For users working on cime machines you can use the build script to build the
tool. On other machines you'll need to do a port to cime and tell how to build
for that machine. That's talked about in the cime documentation.
And you'll have to make some modifications to the build script.

https://github.com/ESMCI/cime/wiki/Porting-Overview

Machines that already run CTSM or CESM have been ported to cime. So if you can
run the model on your machine, you will be able to build the tool there.

To get a list of the machines that have been ported to cime: 

``` shell
# Assuming pwd is the tools/mksurfdata_esmf directory
cd ../../cime/scripts  # or ../../../../cime/scripts for a CESM checkout
./query_config --machines
```

#### NOTE:
In addition to having a port to cime, the machine also needs to have PIO built
and able to be referenced with the env variable PIO which will need to be in
the porting instructions for the machine. An independent PIO library
is available on supported CESM machines.

<!-- ============================================== -->
#### IMPORTANT NOTE: ONLY WORKING ON DERECHO CURRENTLY
<!-- ============================================== -->


> [!IMPORTANT]
> Currently we have run and tested mksurfdata_esmf on Derecho. Please see this github issue about mksurfdata_esmf on other CESM machines:

https://github.com/ESCOMP/CTSM/issues/2341

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

