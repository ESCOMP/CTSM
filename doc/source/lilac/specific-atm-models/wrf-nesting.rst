.. _wrf:

.. highlight:: shell

========================================
 Using CTSM with WRF (Nested Model Runs)
========================================

This section includes instructions on how to run WRF coupled with CTSM for a nested domain.

A nested domain is usually used to have a finer-resolution domain within the coarser model domain. A nested simulation enables running at a higher resolution over a smaller domain.

.. note::
    A nest should cover a portion of the parent domain and is fully contained by
    the parent domain, so that it is driven along its lateral boundaries by the
    parent domain.

.. todo::
    Negin wants to add a flowchart showing the workflow of a nested case.

There are currently two types of nesting available within WRF:

#.  **One-way nesting:**
    In One-way nesting, the boundary conditions are fed to the inner (child) domain from the outer (parent) domain.

#.  **Two-way nesting:**
    In two-way nesting, two things are being done:

        - Exactly similar to 1-way nesting the boundary conditions are fed to the inner domain from the outer (parent) domain.
        - Also, the averaged values from the inner domain are being sent back to the outer domain at the corresponding grid points.

.. important::
    Currently, the WRF-CTSM coupling infrastructure only support one-way nesting.
    This example clarifies the workflow for running a nested WRF-CTSM case using one-way nesting with ``ndown.exe``.

The procedure for running a nested simulation for WRF with CTSM is
similar to the workflow for running WRF real cases, except that it requires additional steps to (1) clone the CTSM repository, (2) build CTSM and LILAC, and (3) define namelist options reuired for CTSM.

A full description of all steps for a WRF-CTSM run are included here.

.. important::

  This section assumes the user has completed all the steps required for
  WRF-CTSM simulation single domain.
  Therefore, we are not repeating the steps necessary for building WRF and
  CTSM.

In this example we use a nested domain over the CONUS as shown below:

.. _Figure ctsm-ndown:

.. todo::
    Replace missing ndown_ctsm_diagram.svg
    
Flowchart for WRF-CTSM one-way nested simulations

Nested Simulations : Pre-processing (geogrid.exe)
-------------------------------------------------
In the WPS/ directory, edit ``namelist.wps`` for a nested simulation over your
desired domains. Make sure to change ``max_dom=2``.

First, use geogrid.exe to define the domain and interpolate static geographical data to the grids::

    ./geogrid.exe >& log.geogrid

This step creates two files, ``geo_em.d01.nc`` and ``geo_em.d02.nc``, which include the domain definition for each domain.

If the geogrid step finishes successfully, you should see the following message in the log file::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Successful completion of geogrid.  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

The basic difference here with a non-nested case is the namelist.wps should have a column for each domain with ``max_dom=2``. For example:

::
    
    &share
     wrf_core = 'ARW',
     max_dom = 2,
     start_date = '2013-04-01_00:00:00','2013-04-01_00:00:00',
     end_date   = '2013-04-30_00:00:00','2013-04-30_00:00:00',
     interval_seconds = 21600
     io_form_geogrid = 2,
    /

    &geogrid
     parent_id         =   1,   1,
     parent_grid_ratio =   1,   3,
     i_parent_start    =   1,  61,
     j_parent_start    =   1,  57,
     e_we              =  200, 103,
     e_sn              =  140, 103,

Therefore ``geogrid.exe`` creates two files corresponding to each domain.

Nested Simulations : Pre-processing (ungrib.exe)
-------------------------------------------------
As mentioned previously, the purpose of the ungrib script is to unpack GRIB meteorological data and pack it into an intermediate file format. This step is exactly identical to a non-nested simulation.

Run ungrib to get gribbed data into usable format to be ingested by WRF.

To run ungrib.exe, first link the GRIB data files that are going to be used::

    ./link_grib.csh $your_GRIB_data_path

Based on your GRIB data type, link or copy the appropriate VTable to your WPS directory. WRF has some prepared VTable under ``/ungrib/Variable_tables/`` folder.

Extract meteorological fields from GRIB-formatted files::

    ./ungrib.exe >& log.ungrib

Check ungrib log for the following message showing successful completion of ungrib step::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Successful completion of ungrib.   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

At this point, you should see ungrib output (intermediate files) in your WPS directory.

Nested Simulations : Pre-processing (metgrid.exe)
-------------------------------------------------
Ensure that the `start_date` and `end_date` for domain two is set correctly for your simulation. Next, run ``metgrid.exe``::

    ./metgrid.exe >& log.metgrid

Check the metgrid log for the following message showing successful completion of metgrid step::

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !  Successful completion of metgrid.  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Running metgrid for two domains will create files like below::

    met_em.d01.*
    met_em.d02.*

Nested Simulations : real.exe
------------------------------

In this step, run ``real.exe`` to generate initial and boundary conditions for both domains.

In summary, complete the following steps:

Move or link WPS output files (``met_em.d01*`` and ``met_em.d02`` files) to your WRF test directory.

Edit namelist.input for your WRF domain and desirable configurations. This should be the same domain as WPS namelist. Make sure you set ``max_dom = 2,`` in the namelist.

To run WRF-CTSM, in your namelist change land-surface option to 6 for both domains::

    sf_surface_physics = 6, 6,

Run real.exe (if compiled parallel submit a batch job) to generate initial and boundary condition files for both domain. Make sure the following three files have been created in your directory::

    wrfinput_d01
    wrfinput_d02
    wrfbdy_d01

The boundary condition file is only created for the outer domain.

Check the last line of the real log file for the following message:

.. todo:: What message?

Rename wrfinput_d02
-------------------
Next, rename the ``wrfinput_d02`` file to ``wrfndi_d02``::

    mv wrfinput_d02 wrfndi_d02

Run ndown.exe
-------------
In this step, we run ndown.exe to create initial and boundary condition for domain 2 based on the domain 1 (outer domain).

Add the following into your namelist.input file under ``&time_control``::

     io_form_auxinput2 = 2

Run ndown.exe to create ``wrfinput_d02`` and ``wrfbdy_d02``.

Run WRF for coarser domain
---------------------------
In this step, run WRF for the outer domain. Make sure that ``max_dom = 1`` to run only for the coarser domain.

This step is exactly identical as the previous example and only creates the ``wrfout*`` files for the coarser domain.

Please make sure to copy ``lnd_in`` , ``lilac_in``, and ``lnd_modelio`` for the coarser domain in this directory.

Create CTSM runtime files for the fine domain
---------------------------------------------
This step is in addition creating CTSM runtime files for coarser domain which was explained here. For succesfully completing the previous step you should have already created these files for the coarser domain.

.. seealso::

    The instructions for setting CTSM runtime options, are discussed in depth
    in section :numref:`setting-ctsm-runtime-options`. For creating the runtime
    files for the finer domain you should follow the steps in section
    :numref:`setting-ctsm-runtime-options`.

Again, the goal here is to create files that determine CTSM runtime options which are defined within these three files:

- ``lnd_in``: This is the main namelist input file for CTSM inner domain

- ``lnd_modelio.nml``: This sets CTSM's PIO (parallel I/O library) configuration settings

- ``lilac_in``: This namelist controls the operation of LILAC

Run WRF for the finer domain
-----------------------------
First, save (rename or move) the data from the coarser domain simulation (``wrfout_d01_*`` files). Next, rename ``wrfinput_d02`` and ``wrfbdy_d02`` to ``wrfinput_d01`` and ``wrfbdy_d01``, respectively.

Edit namelist.input, moving all of the fine-grid domain data from column 2 to column 1 so that this run will be for the fine-grid domain only. Make sure you set ``max_dom=1`` and set your ``time_step`` based on the finer domain.

.. note::
    It may be beneficial to save namelist.input to something else prior to this step in case you need to repeat this
    process in the future. Save the newly-edited namelist as namelist.input .

Now run wrf.exe by submitting a job similar to a not-nested case.

.. important::

    The output for the finer domain is wrfout_d01_* not wrfout_d02_* and although
    in the name it is saying d01 it is technically d02 domain.

