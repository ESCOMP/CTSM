.. _wrf:

.. highlight:: shell

=====================
 Using CTSM with WRF
=====================

This section describes the procedure for building and running the CTSM
library and its dependencies, and linking to these libraries in the WRF
model's build via LILAC. As such this section repeats some information
from earlier sections but in recipe form and with minimal detail.

.. important::

  This section assumes use of a machine that has been ported to CIME.
  In this example we assume NCAR’s cheyenne computer in particular.

Preparing the CTSM
==================

.. todo::

  I think we don't need some of the following since too much instructions
  are usually more confusing...
  I tried removing incorrect or redundant information.

  For example, people who do not use cheyenne might get confused by the 
  following commands on scratch.

Decide where you will work, for example::

  cd /glade/scratch/$USER
  mkdir git_wrf_ctsm
  cd git_wrf_ctsm

.. note::

  Discs other than /glade/scratch may provide insufficient space for
  output from simulations longer than one or two months.


Obtain CTSM by running::

  git clone https://github.com/ESCOMP/ctsm.git
  cd ctsm
  git checkout lilac_cap
  ./manage_externals/checkout_externals -v

Build CTSM and its dependencies, for example for cheyenne::

  ./lilac/build_ctsm /glade/scratch/$USER/ctsm_build_dir --compiler intel --machine cheyenne


Set environment similar to environments used for your CTSM build for bash::

  source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.sh

or Cshell::

  source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.csh


.. note::

  For additional details on preparing the CTSM, including how to
  recompile when making code changes to the CTSM, read section
  _obtaining-and-building-ctsm. <-- CREATED LINK TO THE CORRECT SECTION?

Building the WRF model with CTSM
=======================
.. todo::

  update the git address to WRF feature branch...

Clone WRF CTSM branch into your directory::

  git clone git@github.com:billsacks/WRF.git
  cd WRF
  git checkout lilac_dev


For building WRF using CTSM, we should set makefile variables from CTSM needed for 
WRF build by (BASH)::

  export WRF_CTSM_MKFILE=/glade/scratch/$USER/ctsm_build_dir/bld/ctsm.mk

or::

  setenv WRF_CTSM_MKFILE /glade/scratch/$USER/ctsm_build_dir/bld/ctsm.mk

.. todo::
    Bill and Sam do we need the following still:?

The following is needed in order to undo an undesired setting in that env_mach_specific file::

  export MPI_USE_ARRAY=None

or::

  setenv MPI_USE_ARRAY None

There are also few other environmental setting that should be set for building WRF.
Some of these are not required, but might help if you face any compilation errors.

Explicitly define which model core to build by::

    export WRF_EM_CORE=1

or::

    setenv WRF_EM_CORE 1

Explicilty turn off data assimilation by::

    export WRF_DA_CORE=0

or::

    setenv WRF_DA_CORE 0


Make sure you set NETCDF environment variable by::

    setenv NETCDF /usr/local/netcdf/ (or wherever you have netcdf compiled.)

Then configure and build WRF for your machine and intended compiler by::

    ./clean -a
    ./configure

Choose one of the options, similar to the compiler used for building CTSM.

Next, choose one of the options for nesting. Currently nesting is not available for WRF-CTSM,
therefore we should use 1.

Then compile em_real and save the log::

    ./compile em_real >& compile.log


.. note::

    The ./compile step might take more than 30 minutes to complete.


.. note::

    Check the bottom of your log file for a successful compilation message
    or search the file for the string "Error" with a capital E.

.. note::

    Optional: One may use tmux or nohup for configuring and compiling.
    Try "man nohup" for more information.



Compile WRF Preprocessing System (WPS)
==================================================

The WRF Preprocessing System (WPS) is a set of programs to prepare
input to the real program for WRF real-data simulations.

.. note::
    Building WPS requires that WRF be already built successfully.


Get WPS zipped tar file from: 
http://www2.mmm.ucar.edu/wrf/users/download/get_source.html

Untar WPS tar file::

    gzip -cd WPSV4.0.TAR.gz | tar -xf -


Then we should compile WPS similar to the way we build WRF. In summary::

    cd WPS
    ./configure

Here choose one option, for your intended compiler, similar to your WRF build.
After configuring, you can check configure.wps for making sure all the libs and paths 
are set correctly.

Then, compile WPS::
    ./compile >& compile.log

.. note::
    If wps build is succsfully you should see geogrid.exe, ungrib.exe, and metgrid.exe.
    Alternatively, you can check the log for successful build message.


Run WRF Preprocessing System (WPS) Steps
==================================================

Edit namelist.wps for your domain of interest, which should be the same
domain as used in your WRF namelist.

Define the domain and interpolate static geographical data to the grids::

  ./geogrid.exe >& log.geogrid

Link in the input GFS data files::

  ./link_grib.csh $path_where_you_placed_GFS_files

Extract meteorological fields from GRIB-formatted files::

  ./ungrib.exe

Horizontally interpolate the metrological fields extracted by ungrib to
the model grids defined in geogrid::

  ./metgrid.exe >& log.metgrid

You should now have met_em.d01* files.


Run Real program
==================================================
Run real.exe to generate initial and boundary conditions. 

Follow WRF instructions for creating initial and boundary
conditions. In summary, complete the following steps: 

Move or link WPS output files (met_em.d01* files) to your WRF/run directory. 

Edit namelist.input for your WRF domain and desirable configurations.
This should be the same domain as in the namelist used in WPS. 
To run WRF-CTSM, change land-surface option to 51::

  sf_surface_physics = 51

.. note::

  sf_surface_physics values for running WRF-Noah and WRF-NoahMP are
  2 and 4, respectively.

Run real.exe (if compiled parallel submit a batch job) to generate
wrfinput and wrfbdy files.


Create input namelists for CTSM and LILAC
=========================================

Introduce the following diffs to ./git_wrf_ctsm/ctsm/lilac/atm_driver/<file>
by replacing the entries preceded by minus signs with the entries
preceded by plus signs.

diff ./lilac/atm_driver/atm_driver_in ./lilac/atm_driver/atm_driver_in:

.. code-block:: diff

  -  atm_mesh_file = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc'
  -  atm_global_nx = 72
  -  atm_global_ny = 46
  +  atm_mesh_file = '/glade/work/slevis/barlage_wrf_ctsm/conus/mesh/wrf2ctsm_land_conus_ESMFMesh_c20191216.nc'
  +  atm_global_nx = 199
  +  atm_global_ny = 139

diff ./lilac/atm_driver/ctsm.cfg ./lilac/atm_driver/ctsm.cfg:

.. code-block:: diff

  -configuration     = clm
  -structure         = standard
  -clm_bldnml_opts   = -bgc sp
  -gridmask          = gx3v7
  -lnd_grid          = 4x5 
  -lnd_domain_file   = domain.lnd.fv4x5_gx3v7.091218.nc
  -lnd_domain_path   = /glade/p/cesmdata/cseg/inputdata/share/domains
  -clm_namelist_opts = hist_nhtfrq=-24 hist_mfilt=1 hist_ndens=1
  +configuration     = nwp
  +structure         = fast
  +clm_bldnml_opts   = -bgc sp -clm_usr_name wrf2ctsm
  +gridmask          = null
  +lnd_grid          = wrf2ctsm
  +lnd_domain_file   = domain.lnd.wrf2ctsm_lnd_wrf2ctsm_ocn.191211.nc
  +lnd_domain_path   = /glade/work/slevis/barlage_wrf_ctsm/conus/gen_domain_files
  +clm_namelist_opts = hist_nhtfrq=1 hist_mfilt=1 hist_ndens=1 fsurdat="/glade/work/barlage/ctsm/conus/surfdata_conus/surfdata_conus_hist_16pfts_Irrig_CMIP6_simyr2000_c191212.nc" finidat="/glade/scratch/sacks/wrf_code/WRF/test/em_real/nldas_nwp_0109a.clm2.r.2000-04-01-64800.nc" use_init_interp=.true.

diff ./lilac/atm_driver/lilac_in ./lilac/atm_driver/lilac_in:

.. code-block:: diff

  - atm_mesh_filename = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc'
  + atm_mesh_filename = '/glade/work/slevis/barlage_wrf_ctsm/conus/mesh/wrf2ctsm_land_conus_ESMFMesh_c20191216.nc'

  - lnd_mesh_filename = '/glade/p/cesmdata/cseg/inputdata/share/meshes/fv4x5_050615_polemod_ESMFmesh.nc'
  + lnd_mesh_filename = '/glade/work/slevis/barlage_wrf_ctsm/conus/mesh/wrf2ctsm_land_conus_ESMFMesh_c20191216.nc'

Generate the lnd_in file by running the following from
./git_wrf_ctsm/ctsm/lilac/atm_driver::

  ../../lilac_config/buildnml 

Copy lilac_in, lnd_in, and lnd_modelio.nml to the WRF/run directory.




Run WRF
=================

If real program is completed successfully, we should see wrfinput and wrfbdy files
in our directory.

Next, we should run WRF via batch job.
For Cheyenne, we should submit a batch job to PBS (Pro workload management system).
For more instructions on running a batch job on Cheyenne, please check:
https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs


A sample of basic PBS job for Cheyenne::

    #!/bin/tcsh
    #PBS -N job_name
    #PBS -A project_code
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


