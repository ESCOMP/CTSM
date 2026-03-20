.. include:: ../substitutions.rst

.. _creating-surface-datasets:

===========================
 Creating Surface Datasets
===========================

<!-- ======= -->
## mksurfdata_esmf Purpose
<!-- ======= -->
This tool is intended to generate fsurdat files (surface datasets) and landuse files for the
CTSM. It can generate global, regional, and single-point fsurdat files, as long
as a mesh file is available for the grid.

The subset_data tool allows users to make fsurdat files from existing fsurdat
files when a mesh file is unavailable. Generally, users are encouraged to use the
subset_data tool for generating regional and single-point fsurdat files (see note about regional grids at the end).

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

The ``$CTSMROOT/tools/README.md`` goes through the complete process for creating input files needed to run CLM. We repeat that file here:

.. include:: ../../../../tools/README.md
   :code: markdown

