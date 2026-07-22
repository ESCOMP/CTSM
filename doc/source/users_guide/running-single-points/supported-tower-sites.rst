.. include:: ../substitutions.rst

.. _supported-tower-sites:

********************************************
Supported tower sites for single-point runs
********************************************

CTSM has functionality within the ``run_tower`` tool for running single-point cases at particular supported tower sites using forcing data from those sites.

This tool was developed as a collaboration between NCAR's modeling capabilities and NEON's measurement network that could drive scientific discovery at the confluence of geosciences and biological sciences. The tool was then expanded to include PLUMBER2 sites to support a wider variety of ecological research projects.

Broadly, this tool can be used to probe questions such as:

    * What biases in NCAR models can current observations address?
    * How can NCAR models inform observational data streams?
    * What new hypotheses of atmospheric science and macroscale ecology can be tested with observations and NCAR models to increase our understanding of the biosphere-atmosphere system and its response to global environmental change?
    * Can Earth system prediction be extended to ecological forecasts?

====================================================
General Information on Running Supported Tower Sites
====================================================

The ``run_tower`` capability allows users to run Community Land Model (CLM) simulations at NEON and PLUMBER2 tower sites in a streamlined manner by setting up the appropriate model configurations, datasets, and initial conditions. This script can run for one or more (NEON or PLUMBER2) tower sites. It will do the following:

    1) Create a generic base case for cloning.
    2) Make the case for the specific neon or plumber site(s).
    3) Make changes to the case, for
        a. AD spinup
        b. post-AD spinup
        c. transient
        d. SASU or Matrix spinup
    4) Build and submit the case.

The available options, a description of those options, and details on default values can be shown by running ``run_tower --help``.

A `tutorial <https://ncar.github.io/ncar-neon-books/notebooks/NEON_Simulation_Tutorial.html>`_ on running and evaluating data from ``run_tower`` is also available.

.. warning:: Note that the run_tower base case must be of same run type as a requested clone, as described by this `issue ticket <https://github.com/ESCOMP/CTSM/issues/1926>`_.

=========================================
NEON Tower Single Point Simulations
=========================================

With this tool, CLM uses gap-filled meteorology from NEON tower sites, the dominant plant species is mapped to the appropriate model plant functional type (PFT), and soil characteristics used in the simulations are updated to match observations from NEON's soil megapits. Gap-filled NEON tower flux data are also available for model evaluation. Additionally, all the commands to run the model are combined into a script that you can easily call from a single line of code.

Currently supported NEON sites can be found by running ``run_tower --help``.

.. note:: If you choose to run ``all``, single point simulations at all NEON sites will be run. This is a useful feature, but we recommend testing out running just one site first.

Information on the specific sites can be found on the `NEON webpage <https://www.neonscience.org/field-sites>`_.

.. note:: For NEON tower site simulations, the default run type is ``transient``.

To run CTSM at a NEON site, change directories to where the run_tower tool is located, and then run the ``run_tower`` command. You can also add any additional arguments as described by the ``help`` options. These steps will look something like this::

 cd CTSM/tools/site_and_regional
 run_tower --neon-sites ABBY

When a simulation completes, the data are stored in the archive directory under ``CTSM/tools/site_and_regional/archive``. In this directory you will find files that include data for every day of the simulation, as well as files that average model variables monthly. The output file names are automatically generated and are composed of the simulation name, which includes the site name, type of simulation (eg, ``transient``), and the date of simulated data.
The tower simulations generate two types of files:

1) ``h0a`` Variables that are averaged monthly. One file is available for every month of the simulation. These files include hundreds of variables.

2) ``h1a`` Variables that are recorded every 30 minutes. Values are aggregated into one file for each day of the simulation. Each file includes 48 data points for selected variables.

=========================================
PLUMBER2 Tower Single Point Simulations
=========================================

.. note:: A few important notes regarding the PLUMBER2 tower site simulations are that the default run type is ``ad``; additionally, PLUMBER2 cases all start in different years; and the atmospheric forcing data is in local time, consequently, the model is started at GMT time corresponding to local midnight.

Currently, the ``run_tower`` tool supports running CTSM at PLUMBER2 sites using forcing data from the PLUMBER2 project. Detailed site information is provided in `Ukkola et al. 2022 <https://doi.org/10.5194/essd-14-449-2022>`_ , and the description of the experiment and results are provided in `Abramowitz et al. 2024 <https://doi.org/10.5194/bg-21-5517-2024>`_.
Information on the specific sites can be found `here <https://researchdata.edu.au/plumber2-forcing-evaluation-surface-models/1656048>`_.

Currently supported PLUMBER2 sites can be found by running ``run_tower --help``. Keep in mind that the experiment was designed to run 170 sites, however, `Abramowitz et al. 2024 <https://doi.org/10.5194/bg-21-5517-2024>`_ identified issues with some of the sites (e.g., sites with precipitation reported in one unit while the metadata reported a different unit), and most of that article only uses 156 sites.

To run CTSM at a PLUMBER2 site, change directories to where the ``run_tower`` tool is located, and then run the ``run_tower`` command. You can also add any additional arguments as described by the ``--help`` option. In the simplest form, these steps will look like this:
::

    conda activate ctsm_pylib
    cd CTSM/tools/site_and_regional
    run_tower --plumber-sites AR-SLu

The history output for a PLUMBER2 case will be set up and archived similarly to the output for a NEON case, as described above.

--------------------
General Notes
--------------------

A few points regarding the PLUMBER2 simulations using the ``run_tower`` tool:

#. The PLUMBER2 simulations make extensive use of usermods (located at ``cime_config/usermods_dirs/clm/PLUMBER2``) to configure the simulation for each site and type of simulation. This is accomplished by shell commands (``shell_commands``) and user_nl_* files (e.g., ``user_nl_clm`` and ``user_nl_datm_streams``).  The usermods have already been generated by default using the ``tools/site_and_regional/plumber2_usermods`` script which reads a csv file (``tools/site_and_regional/PLUMBER2_sites.csv``) and generates the usermods.

#. The surface datasets have already been generated for all sites and are available in the local inputdata directory. If desired, users can create their own surface datasets by first modifying the csv file (``tools/site_and_regional/PLUMBER2_sites.csv``) and then using the ``python/ctsm/site_and_regional/plumber2_surf_wrapper.py`` tool which makes use of ``subset_data`` to generate the surface dataset.

#. When first becoming familiar with the PLUMBER2 simulations, we suggest using the ``--setup-only`` flag to examine how the site files are being configured before running the simulations. Keep in mind that if you use this option you will then need to build the base case manually (using ``./case.build``) because the cloned case uses the build from the base case, and submit the cloned case (``ad``, ``postad``, or ``transient``) manually (using ``./case.submit``).
   ::

       ./run_tower --plumber-sites AR-SLu --setup-only

#. By default, ``baseflow_scalar=0`` is modified in the parameter file for the wetland sites to prevent them from drying out.

#. Currently, the tool is designed only for cases with active biogeochemistry (BGC). Implementing capability for SP simulations is in progress.

--------------------
Spin-up notes
--------------------

#. By default, all simulations (e.g., ad, postad, transient) are setup to run for 1 hour of wall clock time. For the spin-up simulations, it may be necessary to extend the time to match the desired length of the simulations. For example, `Lombardozzi et al. 2023 <https://doi.org/10.5194/gmd-16-5979-2023>`_ used 300 years for ad and 100 years for postad. Users should keep in mind that colder areas might need longer spin-up times. For example, the following will setup and run the AR-SLu site for three 100-year segments in ``ad`` spinup mode for a total of 300 years:
   ::

       cd  CTSM/tools/site_and_regional
       run_tower --plumber-sites AR-SLu --run-type ad --xmlchange STOP_OPTION=nyears,STOP_N=100,RESUBMIT=2,JOB_WALLCLOCK_TIME=06:00:00

#. The initial and final years of the simulations may differ for the spin-up (e.g., ad and postad) and transient simulations. This is to avoid calendar errors encountered when the start time is local midnight.

#. By default, the NO_LEAP calendar is used for the ad and postad cases to avoid issues with the transition between Feb 28th and 29th on leap years, this is accomplished in the usermods via
   ::

       ./xmlchange CALENDAR=NO_LEAP

#. In ad and postad simulations, ``dtlimit`` is set to ``50`` in the user_nl_datm_streams files to avoid issues encountered when cycling over the datm forcing data.

-----------------------
A practical example
-----------------------

Here is an example of setting up and running a PLUMBER2 site with the ``run_tower`` tool, from ad through transient mode.

:: 

    cd CTSM/tools/site_and_regional
    setenv tower AR-SLu  # (tcsh) or export tower=AR-SLu (bash)
    ./run_tower --plumber-sites ${tower} --run-type ad --xmlchange STOP_OPTION=nyears,STOP_N=100,RESUBMIT=2,JOB_WALLCLOCK_TIME=06:00:00
    # CHECK FOR SUCCESSFUL COMPLETION of ad MODE
    ./run_tower --plumber-sites ${tower} --run-type postad --xmlchange STOP_OPTION=nyears,STOP_N=100,JOB_WALLCLOCK_TIME=06:00:00
    # CHECK FOR SUCCESSFUL COMPLETION of postad MODE
    ./run_tower --plumber-sites ${tower} --run-type transient


