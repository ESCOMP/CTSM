.. include:: ../substitutions.rst

.. _supported-tower-sites:

********************************************
Supported tower sites for single-point runs
********************************************

CTSM has functionality within the ``run_tower`` tool for running single-point cases at particular supported tower sites using forcing data from those sites.

This tool was developed as a collaboration between NCAR's modeling capabilities and NEON's measurement network that could drive scientific discovery at the confluence of geosciences and biological sciences. The tool was then expanded to include PLUMBER sites to support a wider variety of ecological research projects.

Broadly, this tool can be used to probe questions such as:

    * What biases in NCAR models can current observations address?
    * How can NCAR models inform observational data streams?
    * What new hypotheses of atmospheric science and macroscale ecology can be tested with observations and NCAR models to increase our understanding of the biosphere-atmosphere system and its response to global environmental change?
    * Can Earth system prediction be extended to ecological forecasts?

====================================================
General Information on Running Supported Tower Sites
====================================================

The ``run_tower`` capability allows users to run Community Land Model (CLM) simulations at NEON and PLUMBER tower sites in a streamlined manner by setting up the appropriate model configurations, datasets, and initial conditions. This script can run for one or more (NEON or PLUMBER) tower sites. It will do the following:

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

1) ``h0`` Variables that are averaged monthly. One file is available for every month of the simulation. These files include hundreds of variables.

2) ``h1`` Variables that are recorded every 30 minutes. Values are aggregated into one file for each day of the simulation. Each file includes 48 data points for selected variables.

=========================================
PLUMBER Tower Single Point Simulations
=========================================

.. note:: A few important notes regarding the PLUMBER tower site simulations are that the default run type is ``ad``; additionally, PLUMBER cases all start in different years.

Currently, the ``run_tower`` tool supports running CTSM at PLUMBER 2 sites using forcing data from the PLUMBER2 projects. Detailed site information is provided in `Ukkola et al. 2022 <https://doi.org/10.5194/essd-14-449-2022>`_ , and the description of the experiment and its results is provided in `Abramowitz et al. 2024 <https://doi.org/10.5194/bg-21-5517-2024>`_.
Information on the specific sites can be found `here <https://researchdata.edu.au/plumber2-forcing-evaluation-surface-models/1656048>`_.

Currently supported PLUMBER Sites can be found by running ``run_tower --help``. Keep in mind that the experiment was designed to run 170 sites; however, they identified different issues with the sites, and most of the article only uses 156 sites.

To run CTSM at a PLUMBER site, change directories to where the run_tower tool is located, and then run the ``run_tower`` command. You can also add any additional arguments as described by the ``help`` options. These steps will look something like this::

 cd CTSM/tools/site_and_regional
 run_tower --plumber-sites AR-SLu

The output for a PLUMBER case will be set up similarly to the output for a NEON case, as described above.

A few notes regarding the PLUMBER 2 simulations using the ``run_tower`` tool:


1) By default, all simulations (e.g., ad, post-ad, transient) are designed to run for 1 hour of wall clock time. For the spin-up simulations, it is recommended to extend the time to match the length of the simulations. 
2) By default, the initial period of the simulations differs for the spin-up (e.g., ad and post-ad) and transient simulations. It is important to revisit for each site the best period to use for the spin-up simulations. It is suggested to run the spin-up simulations using the no_leap calendar to avoid issues with Feb 28/29 on leap years. It is also recommended to set a dtlimit value of 50 in the spin-up simulation of the user_nl_datm_strems files. 
3) By default, the tools call for surfdata files in a default location. These might not be available for all sites. Users can create their own surfdata files, using tools/site_and_regional/plumber2_usermods. 
4) It is suggested to use the flags-- setup-only to examine how the site files are being configured before running the simulations.
5) It is suggested to add ``echo "baseflow_scalar = 0" >> user_nl_clm`` to the user_nl_clm file to avoid issues with the baseflow in the wetland sites.
6) Currently, the tool is designed only for cases with active biochemistry. For SP simulations, it is recommended to use the user mods directory ``/CTSM/cime_config/usermods_dirs/clm/PLUMBER2/${site}``, which contains the necessary modifications to run the simulations in SP mode. 