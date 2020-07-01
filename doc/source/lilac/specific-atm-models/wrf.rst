.. _wrf:

.. highlight:: shell

=====================
 Using CTSM with WRF
=====================

This section describes the procedure for building and running the CTSM
library and its dependencies, and linking to these libraries in the WRF
model's build via LILAC. As such this section repeats some information
from earlier sections but in recipe form and with minimal explanation.

.. important::

  This section assumes use of a machine that has been ported to CIME.
  In this example we assume NCAR’s cheyenne computer in particular.

Clone CTSM Repository and Build CTSM
------------------------------------

Clone the CTSM repository::

    git clone https://github.com/ESCOMP/CTSM.git
    cd CTSM
    git checkout lilac_cap
    ./manage_externals/checkout_externals

.. todo::

    Remove "git checkout lilac_cap" from the above when ready

Build CTSM and its dependencies based on the instructions from previous sections,
for example on Cheyenne::

    ./lilac/build_ctsm /glade/scratch/$USER/ctsm_build_dir --compiler intel --machine cheyenne

Source ctsm_build_environment.sh (bash environment)::

    source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.sh

or ctsm_build_environment.csh (Cshell environment):

.. code-block:: Tcsh

    source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.csh

.. note::

  For further detail on preparing the CTSM, including how to
  recompile when making code changes to the CTSM, read Section 3.2:
  https:../obtaining-building-and-running/index.html
  By the way, do not let Section 3.2.2 confuse you. We address that step
  right after compiling the WRF model (next).

Building WRF with CTSM
----------------------

.. todo::

    update the git address to WRF feature branch...

Clone the WRF CTSM feature branch::

    git clone https://github.com/billsacks/WRF.git
    cd WRF
    git checkout lilac_dev


Set makefile variables from CTSM needed for the WRF build by setting the following environment.
For example for Bash::

    export WRF_CTSM_MKFILE=/glade/scratch/$USER/ctsm_build_dir/bld/ctsm.mk


.. code-block:: Tcsh

    setenv WRF_CTSM_MKFILE /glade/scratch/$USER/ctsm_build_dir/bld/ctsm.mk

The next two environment settings for building WRF may help if you
encounter compilation errors, but should be unnecessary for completing
the current example on cheyenne.

Explicitly define which model core to build by::

    export WRF_EM_CORE=1

Explicilty turn off data assimilation by::

    export WRF_DA_CORE=0

Now configure and build WRF for your machine and intended compiler.
The ./clean command is necessary after any modification of WRF code::

    ./clean -a
    ./configure

At the prompt choose one of the options, similar to the compiler used 
for building CTSM. The specific example has been tested successfully by
choosing 15 here.

.. todo::

    Negin, by "similar to" do you mean "same as" in the above?

The next prompt requests an option for nesting. Currently nesting is not
available for WRF-CTSM so enter 1.

Now compile em_real and save the log::

    ./compile em_real >& compile.log


.. note::

    Optional: One may use tmux or nohup for configuring and compiling.
    Try "man nohup" for more information.

.. note::

    Check the bottom of your log file for a successful compilation message
    or search the file for the string "Error" with a capital E.

.. note::

    The ./compile step may take more than 30 minutes to complete.
    While you wait, follow the instructions in Section 3.2.2 (next)

Now follow the instructions in this Section::

 https:../obtaining-building-and-running/setting-ctsm-runtime-options.html

In step 3 of that Section we used for this example::

 lnd_domain_file = /glade/work/slevis/barlage_wrf_ctsm/conus/gen_domain_files/domain.lnd.wrf2ctsm_lnd_wrf2ctsm_ocn.191211.nc
 fsurdat = /glade/work/slevis/git_wrf/ctsm_surf/surfdata_conus_hist_16pfts_Irrig_CMIP6_simyr2000_c191212.nc
 finidat = /glade/work/slevis/git_wrf/ctsm_init/finidat_interp_dest_wrfinit_snow_ERAI_12month.nc

In step 4 of that Section we used for this example::

 atm_mesh_filename = '/glade/work/slevis/barlage_wrf_ctsm/conus/mesh/wrf2ctsm_land_conus_ESMFMesh_c20191216.nc'
 lnd_mesh_filename = '/glade/work/slevis/barlage_wrf_ctsm/conus/mesh/wrf2ctsm_land_conus_ESMFMesh_c20191216.nc' 

In step 6 of that Section you will copy some files to your WRF/run
directory. Then you will be ready to continue.

.. note::

 If you wish to merge your WRF initial conditions from a wrfinput file
 into the existing CTSM initial condition file, complete the following step.

Type::

 module load ncl
 ncl transfer_wrfinput_to_ctsm_with_snow.ncl 'finidat="the_existing_finidat_file.nc"' 'wrfinput="your_wrfinput_file"' 'merged="the_merged_finidat_file.nc"'

.. todo::

 Make the above ncl script available.

Compile WRF Preprocessing System (WPS)
--------------------------------------

The WRF Preprocessing System (WPS) is a set of programs to prepare
inputs to the real program executable (real.exe) for WRF real-data simulations.
If you wish to complete the offered example with preexisting inputs, then
skip to the next section, which is titled "Run WRF."

.. note::

    Building WPS requires that WRF be already built successfully.


Get WPS from:

https://www2.mmm.ucar.edu/wrf/users/download/wrf-regist_or_download.php

New users must complete a registration form in this step.

Then compile WPS similar to the way WRF was built. In summary::

    cd WPS
    ./configure

At the prompt choose your intended compiler, similar to your WRF build.
After configuring, check configure.wps to make sure all the libs and paths 
are set correctly.

.. todo::

    Negin, by "similar to" do you mean "same as" in the above?

Then, compile WPS::

    ./compile >& compile.log

.. note::

    If wps builds succesfully you should see geogrid.exe, ungrib.exe, and metgrid.exe.
    Alternatively, you can check the log for successful build message.


Run WPS
-------

Edit namelist.wps for your domain of interest, which should be the same
domain as used in your WRF namelist.

First, use geogrid.exe to define the domain and interpolate static geographical data
to the grids::

    ./geogrid.exe >& log.geogrid

If the geogrid step finishes successfully, you should see the following message in the log file::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Successful completion of geogrid.  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Next, run ungrib to get gribbed data into usable format to be ingested by WRF.

To run ungrib.exe, first link the GRIB data files that are going to be used::

    ./link_grib.csh $your_GRIB_data_path

Based on your GRIB data type, link or copy the appropriate VTable to your WPS directory.
WRF has some prepared VTable under /ungrib/Variable_tables/ folder.

Extract meteorological fields from GRIB-formatted files::

    ./ungrib.exe >& log.ungrib

Check ungrib log for the following message showing successful completion of ungrib step::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Successful completion of ungrib.   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


At this point, you should see ungrib output (intermediate files) in your WPS directory.

Horizontally interpolate the meteorological fields extracted by ungrib to
the model grids defined in geogrid::

    ./metgrid.exe >& log.metgrid


Check the metgrid log for the following message showing successful completion of
metgrid step::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Successful completion of metgrid.  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



Run real.exe
------------

Run real.exe to generate initial and boundary conditions.

Follow WRF instructions for creating initial and boundary conditions. 
In summary, complete the following steps:

Move or link WPS output files (met_em.d01* files) to your WRF/run directory. 

Edit namelist.input for your WRF domain and desirable configurations.
This should be the same domain as in the namelist used in WPS.


.. todo::

    update the option number of wrf namelist.


To run WRF-CTSM, change land-surface option to 51::

  sf_surface_physics = 51

.. note::

  sf_surface_physics values for running WRF-Noah and WRF-NoahMP are
  2 and 4, respectively.

.. todo::

    add the link and adding some note that nested run is not possible....

Run real.exe (if compiled parallel submit a batch job) to generate
wrfinput and wrfbdy files.


Check the last line of the real log file for the following message::

    SUCCESS COMPLETE REAL_EM INIT


Run WRF
-------

If real.exe completed successfully, we should have wrfinput and wrfbdy files
in our directory. If you plan to use this example's preexisting files, copy
the following files to your WRF/run directory::

 /glade/work/slevis/git_wrf/WRF/test/em_real/namelist.input.ctsm.2013.d01.12month
 /glade/work/slevis/git_wrf/WRF/test/em_real/wrfinput_d01.ERAI.12month
 /glade/work/slevis/git_wrf/WRF/test/em_real/wrfbdy_d01.ERAI.12month

Now run WRF-CTSM. On Cheyenne this means submitting a batch job to PBS (Pro workload management system).
For detailed instructions on running a batch job on Cheyenne, please check:
https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs

A simple PBS script to run WRF-CTSM on Cheyenne looks like this:

.. code-block:: Tcsh

    #!/bin/tcsh
    #PBS -N your_job_name
    #PBS -A your_project_code
    #PBS -l walltime=01:00:00
    #PBS -q queue_name
    #PBS -j oe
    #PBS -k eod
    #PBS -m abe
    #PBS -M your_email_address
    #PBS -l select=2:ncpus=36:mpiprocs=36

    ### Set TMPDIR as recommended
    setenv TMPDIR /glade/scratch/$USER/temp
    mkdir -p $TMPDIR

    ### Run the executable
    mpiexec_mpt ./wrf.exe

If you named this script run_wrf_ctsm.csh, submit the job like this::

    qsub run_wrf_ctsm.csh

If your terminal windows have logged off, repeat
source ctsm_build_environment.sh (bash environment) before submitting
the job::

    source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.sh

or ctsm_build_environment.csh (Cshell environment):

.. code-block:: Tcsh

    source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.csh

