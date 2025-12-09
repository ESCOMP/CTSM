.. include:: ../substitutions.rst

.. _single-point-regional-configurations:

*****************************************************
Introduction to Single-Point and Regional Grid Setups
*****************************************************

CTSM is designed to support a wide range of spatial scales, ranging from global simulations to regional runs to highly resolved single-point cases. Setting up and running single-point and regional simulations is useful for a variety of purposes including: running quick cases for testing, evaluating specific vegetation types, or running with observed data from a specific site to generate and test hypotheses.

Single-point cases allow users to run CTSM at a specific location such as a flux tower or ecological field site. Single-point runs are especially useful where high-resolution meteorological forcing data and site-specific observations are available and require minimal computational resources.

Regional configurations support simulations over broader geographic areas defined by a user-specified domain. Regional runs require additional input data such as meteorological forcing for the region. You can either extract regional subsets from global datasets or create custom datasets for your region of interest.

.. _options-for-single-points:

=========================================
 Choosing the right single point options
=========================================

There are several different ways to set up single-point and regional cases.

For supported tower sites:
--------------------------

You can run at a supported tower site if one of the supported single-point/regional datasets is your site of interest (see :ref:`supported-tower-sites`). All the datasets are created for you, and you can easily select one and run it out of the box using a supported resolution from the top level of the CESM scripts. You can also use this method for your own datasets, but you have to create the datasets, and add them to the XML database in scripts, CLM and to the DATM. This is worthwhile if you want to repeat many multiple cases for a given point or region.

Next, using ``subset_data`` is the best way to setup cases quickly where you can use a simple tool to create your own datasets (see :ref:`generic_single_point_runs`). With this method you don't have to change DATM or add files to the XML database. ``subset_data`` will create a usermod directory where you can store your files and the files needed to directly run a case.

For unsupported tower sites:
----------------------------

If you have meteorology data that you want to force your CLM simulations with, you'll need to setup cases as described in :ref:`pre-defined-single-pt-regional-resolutions`. You'll need to create CLM datasets either according to ``CLM_USRDAT_NAME``. You may also need to modify DATM to use your forcing data. And you'll need to change your forcing data to be in a format that DATM can use.

================
Spinning up CTSM
================

We make steady state assumptions about the initial state of ecosystem properties including temperature water, snow, ice, carbon & nitrogen. This is the equilibrium state of the model, given the forcing data. Spinning up the model brings internal state variables into equilibrium with environmental forcing conditions so that the results are not influenced by the initial conditions of state variables (such as soil C). In runs with active biogoechemistry, we need to get the ecosystem carbon and nitrogen pools with long turnover times into steady state.

Specifically, spinning up CTSM consists of 3 parts including:

1. AD, or accelerated decomposition: The turnover and decomposition of the slow pools of C and N that normally have a long residence time in ecosystems is mathematically accelerated, where we make the slow pools spin up more quickly by increasing their turnover time. This includes:
-Accelerating turnover of wood, litter and soil pools
-Accelerating advection and diffusion terms
-Calculating this as a function of latitude so that spinup is more accelerated in high latitude regions.

2. postAD, which occurs after AD spinup: During postAD runs we take away accelerated decomposition and let the ecosystem settle into its equilibrium, or steady-state under 'normal conditions'. Pools are increased by the same degree their turnover was increased (e.g., turnover 10x faster means the pool must be 10x larger). During AD and postAD spinup we cycle over several years of input data and hold other inputs constant (e.g., atmospheric CO2 concentrations, N deposition, etc.). For transient runs these inputs also change over time.

3. transient: Transient runs are used to compare with observations, and include high frequency output that we can compare with flux tower measurements. The end of the spinup simulation is used as the initial conditions for a transient simulation, set in the user_nl_clm file.




