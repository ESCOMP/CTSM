.. _wrf:

.. highlight:: shell

=====================
 Using CTSM with WRF
=====================

This section includes instructions on how to run WRF coupled with CTSM via LILAC
framework.

The procedure for running WRF with CTSM is similar to the
workflow for running WRF real cases, except that it requires
additional steps to (1) clone the CTSM repository, (2) build
CTSM and LILAC, and (3) define namelist options reuired for CTSM.

A full description of all steps for a WRF-CTSM run are included here.

Specific new steps that would not be completed in a standard WRF real case
are described in sections :numref:`clone-WRF-CTSM-repositories`,
:numref:`build-CTSM-and-dependencies` , 
and :numref:`wrf-set-ctsm-runtime-options`.

.. important::

  This section assumes use of a machine that has been ported to CIME.
  If CIME is not ported to your machine, please see `instructions on porting CIME
  <https://esmci.github.io/cime/versions/master/html/users_guide/porting-cime.html#porting>`_.

  In this example we assume NCAR’s ``Cheyenne`` HPC system in particular.


.. _clone-WRF-CTSM-repositories:

Clone WRF and CTSM Repositories
-------------------------------

Clone the WRF repository and checkout  ``develop`` branch::

    git clone https://github.com/wrf-model/WRF.git WRF-CTSM
    cd WRF-CTSM
    git checkout develop


Clone the CTSM repository::

    git clone https://github.com/ESCOMP/CTSM.git
    cd CTSM
    ./manage_externals/checkout_externals


.. _build-CTSM-and-dependencies:

Build CTSM and its dependencies
-------------------------------

In your CTSM directory, build CTSM and its dependencies based on the 
instructions from section :numref:`obtaining-and-building-ctsm`::

    ./lilac/build_ctsm /PATH/TO/CTSM/BUILD --machine MACHINE --compiler COMPILER

For example on ``Cheyenne`` and for ``Intel`` compiler::

    ./lilac/build_ctsm ctsm_build_dir --compiler intel --machine cheyenne


.. warning::

    The directory you provided for the build script (``/PATH/TO/CTSM/BUILD``) must *not* exist.
    Alternatively, you can use ``--rebuild`` option.

.. note::

    Run ``./lilac/build_ctsm -h`` to see all options available.
    For example if you would like to run with threading support you can use ``--build-with-openmp``.


Building WRF with CTSM
----------------------

First, load the same modules and set the same environments as used for CTSM build by
sourcing ``ctsm_build_environment.sh`` for Bash::

    source ctsm_build_dir/ctsm_build_environment.sh 

or sourcing ``ctsm_build_environment.csh`` for Cshell:

.. code-block:: Tcsh

    source ctsm_build_dir/ctsm_build_environment.csh 

Set makefile variables from CTSM needed for the WRF build by setting the following environment.
For example for Bash::

    export WRF_CTSM_MKFILE=/glade/scratch/$USER/WRF-CTSM/CTSM/ctsm_build_dir/bld/ctsm.mk

or for Cshell:

.. code-block:: Tcsh

    setenv WRF_CTSM_MKFILE /glade/scratch/$USER/WRF-CTSM/CTSM/ctsm_build_dir/bld/ctsm.mk

.. warning::
    Please note that you should point to the absolute path of the ``ctsm.mk`` file.

There are also few other environmental setting that should be set for building WRF.
Some of these are not required, but might help if you face any compilation errors.
For more information check
`WRF Users' Guide <https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_v4/v4.2/WRFUsersGuide_v42.pdf>`_.


Explicitly define the model core to build by (Bash)::

    export WRF_EM_CORE=1

or (Cshell):	

.. code-block:: Tcsh	

    setenv WRF_EM_CORE 1


Explicilty turn off data assimilation by (Bash)::

    export WRF_DA_CORE=0

or (Cshell):	

.. code-block:: Tcsh	

    setenv WRF_DA_CORE 0

Now in your WRF directory configure and build WRF for your machine 
and intended compiler::

    ./clean -a
    ./configure

At the prompt choose one of the options, based on the compiler used
for building CTSM. Then you should choose if you'd like to build serially or
in parallel. For example, you can choose to build with ``intel`` compiler with 
distributed memory parallelization (``dmpar``).

.. tip::

    ``dmpar`` or distributed memory parallelization is the most highly tested and
    recommended for compiling WRF.

The next prompt requests an option for nesting. Currently nesting is not
available for WRF-CTSM so select option ``1 (basic)``.


Now compile em_real and save the log::

    ./compile em_real >& compile.log


Check the bottom of your log file for a successful compilation message.

.. note::

    The ``./compile`` step may take more than 30 minutes to complete.

.. tip::

    Optional: One may use ``tmux`` or ``nohup`` for configuring and compiling.
    Try ``man nohup`` for more information.

Compile WRF Preprocessing System (WPS)
--------------------------------------

The WRF Preprocessing System (WPS) is a set of programs to prepare
inputs to the real program executable (real.exe) for WRF real-data simulations.
If you wish to complete the offered example with preexisting inputs, then
skip to section :numref:`wrf-set-ctsm-runtime-options`.

.. warning::

    Building WPS requires that WRF be already built successfully.


Get WPS from this website::

    https://www2.mmm.ucar.edu/wrf/users/download/wrf-regist_or_download.php

New users must complete a registration form in this step.

Then compile WPS similar to the way WRF was built. In summary::

    cd WPS
    ./configure

At the prompt choose your intended compiler and parallelization method,
similar to the steps in your WRF build.

Then, compile WPS::

    ./compile >& compile.log

.. note::

    If wps builds succesfully you should see ``geogrid.exe``, ``ungrib.exe``, and ``metgrid.exe``.
    Alternatively, you can check the log for successful build messages.


Run WPS Programs
----------------

Edit ``namelist.wps`` for your domain of interest, which should be the same
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
WRF has some prepared VTable under ``/ungrib/Variable_tables/`` folder.

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

Run ``real.exe`` to generate initial and boundary conditions.

Follow WRF instructions for creating initial and boundary conditions. 
In summary, complete the following steps:

Move or link WPS output files (``met_em.d01*`` files) to your WRF test directory. 

Edit namelist.input for your WRF domain and desirable configurations.
This should be the same domain as WPS namelist.


To run WRF-CTSM, in your namelist change land-surface option to 6::

    sf_surface_physics = 6


Run real.exe (if compiled parallel submit a batch job) to generate
``wrfinput`` and ``wrfbdy`` files.


Check the last line of the real log file for the following message::

    SUCCESS COMPLETE REAL_EM INIT

.. _wrf-set-ctsm-runtime-options:


Set CTSM runtime options
------------------------

.. seealso::

    The instructions for setting CTSM runtime options, are discussed in depth
    in section :numref:`setting-ctsm-runtime-options`.

The goal here is to create files that determine CTSM runtime options which
are defined within these three files:

- ``lnd_in``: This is the main namelist input file for CTSM

- ``lnd_modelio.nml``: This sets CTSM's PIO (parallel I/O library) configuration settings

- ``lilac_in``: This namelist controls the operation of LILAC


The basic process for creating the necessary input files are summarized as
follows:

Go to  ``runtime_inputs`` directory::

    cd CTSM/ctsm_build_dir/runtime_inputs

Next, modify and fill in the ``ctsm.cfg`` file to set high-level options to CTSM.
A few options should be filled in; most can be left at their default values or changed if
desired.

The following is the recommended CTSM options to run WRF::

    configuration     = nwp
    structure         = fast
    bgc_mode          = sp

In ``ctsm.cfg`` you should specify CTSM domain file, surface dataset and finidat file.
For this example (US domain), you can use the following settings::

 lnd_domain_file = /glade/work/slevis/git_wrf/ctsm_domain/domain.lnd.wrf2clm_lnd_noneg_wrf2clm_ocn_noneg.201117.nc
 fsurdat = /glade/work/slevis/git_wrf/ctsm_surf/surfdata_conus_hist_16pfts_Irrig_CMIP6_simyr2000_c210119.nc
 finidat = /glade/work/slevis/git_wrf/ctsm_init/finidat_interp_dest_wrfinit_snow_ERAI_12month.nc

File ``user_nl_ctsm`` allows you to override individual CTSM namelist variables
and set any extra namelist items you would like to appear in your ``lnd_in``.
For this example, we recommend adding the following options in
``user_nl_ctsm``::

    use_init_interp = .true.
    init_interp_fill_missing_with_natveg = .true.

Run the script ``make_runtime_inputs`` to create ``lnd_in`` and
``clm.input_data_list``::

    ./make_runtime_inputs

Modify ``lilac_in`` as needed. For this example, you can use the following options::

 atm_mesh_filename = '/glade/scratch/negins/wrf_ctsm_files/wrf2ctsm_land_conus_ESMFMesh_c20201110.nc'
 lnd_mesh_filename = '/glade/scratch/negins/wrf_ctsm_files/wrf2ctsm_land_conus_ESMFMesh_c20201110.nc' 


Run ``download_input_data`` script to download any of CTSM's standard input
files that are needed based on settings in ``lnd_in`` and ``lilac_in``::

    ./download_input_data

Next, copy or link ``lnd_in``, ``lnd_modelio.nml`` and ``lilac_in`` to the direcotory
from which you will be running the model (e.g. WRF/run) ::

    ln -sf /glade/scratch/$USER/WRF-CTSM/CTSM/ctsm_build_dir/runtime_inputs/lnd_in .
    ln -sf /glade/scratch/$USER/WRF-CTSM/CTSM/ctsm_build_dir/runtime_inputs/lilac_in .
    ln -sf /glade/scratch/$USER/WRF-CTSM/CTSM/ctsm_build_dir/runtime_inputs/lnd_modelio.nml .

Run wrf.exe
-----------

If real.exe completed successfully, we should have ``wrfinput`` and ``wrfbdy`` files
in our directory. 

If you plan to use this example's preexisting files, copy
the following files to your WRF run directory::

    cp /glade/scratch/negins/wrf_ctsm_files/namelist.input . 
    cp /glade/scratch/negins/wrf_ctsm_files/wrfinput_d01 .
    cp /glade/scratch/negins/wrf_ctsm_files/wrfbdy_d01 .

Now run WRF-CTSM. On Cheyenne this means submitting a batch job to PBS (Pro workload management system).
Please check NCAR CISL's `instructions on running a batch job on Cheyenne. 
<https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne/running-jobs/submitting-jobs-pbs>`__

A simple PBS script to run WRF-CTSM on ``Cheyenne`` looks like this:

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

    ### Run the executable
    mpiexec_mpt ./wrf.exe

To submit a batch job to the ``Cheyenne`` queues, use ``qsub`` command followed
by the PBS script name. 
For example, if you named this script ``run_wrf_ctsm.csh``, submit the job like this::

    qsub run_wrf_ctsm.csh


