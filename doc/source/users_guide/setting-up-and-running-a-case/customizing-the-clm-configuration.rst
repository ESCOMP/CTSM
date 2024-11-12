.. include:: ../substitutions.rst

.. _configuring-clm:

********************************
 Customizing CLM's Configuration
********************************

The section of the |cesmrelease| Quickstart `CESM Create a Case <https://escomp.github.io/cesm/release-cesm2/quickstart.html#create-a-case>`_ gives instructions on creating a case. Also see a similar section in the CIME User's Guide `CIME Create a case <http://esmci.github.io/cime/users_guide/create-a-case.html>`_. What is of interest here is how to customize your use of CLM for the case that you created.

For CLM when ``preview_namelist``, ``case.build``, or ``case.run`` are called there are two steps that take place:

1. The CLM ``configure`` script is called to setup the build-time configuration for CLM (see :ref:`more-info-clm-config-script`). The env variables for ``configure`` are locked after the ``case.build`` step. So the results of the CLM ``configure`` are locked after the build has taken place.

2. The CLM ``build-namelist`` script is called to generate the run-time namelist for CLM (more information on ``build-namelist`` is given below in :ref:`def-nl-items-and-defaults`).

When customizing your case at the ``case.setup`` step you are able to modify the process by effecting either one or both of these steps. The CLM ``configure`` and ``build-namelist`` scripts are both available in the ``$CTSMROOT/bld`` directory in the distribution. Both of these scripts have a ``-help`` option that is useful to examine to see what types of options you can give either of them.

There are five different types of customization for the configuration that we will discuss: |version| in |cesmrelease| build-time options, |version| in |cesmrelease| run-time options, User Namelist, other noteworthy |cesmrelease| configuration items, the CLM ``configure`` script options, and the CLM ``build-namelist`` script options.

Information on all of the CLM script, configuration, build and run items is found under ``$CTSMROOT/cime_config/config_component.xml``. See `CLM CASEROOT Variable Definitions <http://www.cesm.ucar.edu/models/cesm2/component_settings/clm4_5_input.html>`_.

================================
 CLM Script configuration items
================================

Below we list each of the CESM configuration items that are specific to CLM. All of these are available in your: ``env_build.xml`` and ``env_run.xml`` files.
::

   CLM_ACCELERATED_SPINUP
   CLM_BLDNML_OPTS
   CLM_CO2_TYPE
   CLM_CONFIG_OPTS
   CLM_CPPDEFS
   CLM_FORCE_COLDSTART
   CLM_NAMELIST_OPTS
   CLM_NML_USE_CASE
   CLM_USRDAT_NAME
   COMP_LND

For the precedence of the different options to ``build-namelist`` see the section on precedence below.

The first item ``CLM_CONFIG_OPTS`` has to do with customizing the CLM build-time options for your case, the rest all have to do with generating the namelist.

CLM_CONFIG_OPTS
  The option ``CLM_CONFIG_OPTS`` is all about passing command line arguments to the CLM ``configure`` script. It is important to note that some compsets, may already put a value into the ``CLM_CONFIG_OPTS`` variable. You can still add more options to your ``CLM_CONFIG_OPTS`` but make sure you add to what is already there rather than replacing it. Hence, we recommend using the ``-append`` option to the xmlchange script. In :ref:`more-info-clm-config-script` below we will go into more details on options that can be customized in the CLM ``configure`` script. It's also important to note that the ``$CTSMROOT/cime_config/buildnml`` script may already invoke certain CLM ``configure`` options and as such those command line options are NOT going to be available to change at this step (nor would you want to change them). The options to CLM ``configure`` are given with the ``-help`` option which is given in :ref:`more-info-clm-config-script`... note:: ``CLM_CONFIG_OPTS`` is locked after the ``case.build`` script is run. If you want to change something in ``CLM_CONFIG_OPTS`` you'll need to clean the build and rerun ``case.build``. The other env variables can be changed at run-time so are never locked.

CLM_NML_USE_CASE
  ``CLM_NML_USE_CASE`` is used to set a particular set of conditions that set multiple namelist items, all centering around a particular usage of the model. (See :ref:`precedence-of-opts` for the precedence of this option relative to the others.) To list the valid options do the following:
  ::

     > cd $CTSMROOT
     > ./bld/build-namelist -use_case list

  The output of the above command is:
  ::

     CLM build-namelist - use cases: 1850-2100_rcp2.6_glacierMEC_transient 1850-2100_rcp2.6_transient  \
     1850-2100_rcp4.5_glacierMEC_transient 1850-2100_rcp4.5_transient  \
     1850-2100_rcp6_glacierMEC_transient 1850-2100_rcp6_transient  \
     1850-2100_rcp8.5_glacierMEC_transient 1850-2100_rcp8.5_transient 1850_control  \
     1850_glacierMEC_control 2000-2100_rcp8.5_transient 2000_control 2000_glacierMEC_control  \
     20thC_glacierMEC_transient 20thC_transient glacierMEC_pd stdurbpt_pd
     Use cases are:...

     1850-2100_rcp2.6_glacierMEC_transient = Simulate transient land-use, and aerosol deposition changes  \
     with historical data from 1850 to 2005 and then with the RCP2.6 scenario from IMAGE

     1850-2100_rcp2.6_transient = Simulate transient land-use, and aerosol deposition changes with  \
     historical data from 1850 to 2005 and then with the RCP2.6 scenario from IMAGE

     1850-2100_rcp4.5_glacierMEC_transient = Simulate transient land-use, and aerosol deposition changes  \
     with historical data from 1850 to 2005 and then with the RCP4.5 scenario from MINICAM

     1850-2100_rcp4.5_transient = Simulate transient land-use, and aerosol deposition changes with  \
     historical data from 1850 to 2005 and then with the RCP4.5 scenario from MINICAM

     1850-2100_rcp6_glacierMEC_transient = Simulate transient land-use, and aerosol deposition changes  \
     with historical data from 1850 to 2005 and then with the RCP6 scenario from AIM

     1850-2100_rcp6_transient = Simulate transient land-use, and aerosol deposition changes with  \
     historical data from 1850 to 2005 and then with the RCP6 scenario from AIM

     1850-2100_rcp8.5_glacierMEC_transient = Simulate transient land-use, and aerosol deposition changes  \
     with historical data from 1850 to 2005 and then with the RCP8.5 scenario from MESSAGE

     1850-2100_rcp8.5_transient = Simulate transient land-use, and aerosol deposition changes with  \
     historical data from 1850 to 2005 and then with the RCP8.5 scenario from MESSAGE

     1850_control = Conditions to simulate 1850 land-use
     1850_glacierMEC_control = Running an IG case for 1850 conditions with the ice sheet model glimmer
     2000-2100_rcp8.5_transient = Simulate transient land-use, and aerosol deposition changes with  \
     historical data from 2000 to 2005 and then with the RCP8.5 scenario from MESSAGE

     2000_control = Conditions to simulate 2000 land-use
     2000_glacierMEC_control = Running an IG case for 2000 conditions with the ice sheet model glimmer
     20thC_glacierMEC_transient = Simulate transient land-use, and aerosol deposition changes from 1850  \
     to 2005
     20thC_transient = Simulate transient land-use, and aerosol deposition changes from 1850 to 2005
     glacierMEC_pd = Running an IG case with the ice sheet model glimmer
     stdurbpt_pd = Standard Urban Point Namelist Settings

CLM_BLDNML_OPTS
  The option CLM_BLDNML_OPTS is for passing options to the CLM ``build-namelist`` script. As with the CLM ``configure`` script the CLM $CTSMROOT/cime_config/buildnml may already invoke certain options and as such those options will NOT be available to be set here. The best way to see what options can be sent to the ``build-namelist`` script is to do
  ::

     > cd $CTSMROOT/bld
     > ./build-namelist -help

  Here is the output from the above.
  ::

     ./SYNOPSIS
     build-namelist [options]

     Create the namelist for CLM
     OPTIONS
     -[no-]chk_res            Also check [do NOT check] to make sure the resolution and
                              land-mask is valid.
     -clm_demand "list"       List of variables to require on clm namelist besides the usuals.
                              "-clm_demand list" to list valid options.
                              (can include a list member "null" which does nothing)
     -clm_startfile "file"    CLM restart file to start from.
     -clm_start_type "type"   Start type of simulation
                              (default, cold, arb_ic, startup, continue, or branch)
                              (default=do the default type for this configuration)
                              (cold=always start with arbitrary initial conditions)
                              (arb_ic=start with arbitrary initial conditions if
                               initial conditions don't exist)
                              (startup=ensure that initial conditions are being used)
     -clm_usr_name     "name" Dataset resolution/descriptor for personal datasets.
                              Default: not used
                              Example: 1x1pt_boulderCO_c090722 to describe location,
                                       number of pts, and date files created
     -co2_type "value"        Set CO2 the type of CO2 variation to use.
     -co2_ppmv "value"        Set CO2 concentration to use when co2_type is constant (ppmv).
     -config "filepath"       Read the given CLM configuration cache file.
                              Default: "config_cache.xml".
     -csmdata "dir"           Root directory of CESM input data.
                              Can also be set by using the CSMDATA environment variable.
     -d "directory"           Directory where output namelist file will be written
                              Default: current working directory.
     -drydep                  Produce a drydep_inparm namelist for testing that will go into the
                              "drv_flds_in" file for the driver to pass dry-deposition to the atm.
                              This populates the namelist with valid drydep settings for testing.
                              Default: -no-drydep
                              Note: Can always add drydep fields to user_nl_clm even with --no-drydep
                              (Note: buildnml copies the file for use by the driver)
     -fire_emis               Produce a fire_emis_nl namelist for testing that will go into the
                              "drv_flds_in" file for the driver to pass fire emissions to the atm.
                              This populates the namelist with valid fire-emiss settings for testing.
                              Note: Can always add fire_emis fields to user_nl_clm even with --no-fire_emis
                              (Note: buildnml copies the file for use by the driver)
     -glc_grid "grid"         Glacier model grid and resolution when glacier model,
                              Only used if glc_nec > 0 for determining fglcmask
                              Default:  gland5UM
                              (i.e. gland20, gland10 etcetera)
     -glc_nec <name>          Glacier number of elevation classes [0 | 3 | 5 | 10 | 36]
                              (default is 0) (standard option with land-ice model is 10)
     -glc_smb <value>         Only used if glc_nec > 0
                              If .true., pass surface mass balance info to GLC
                              If .false., pass positive-degree-day info to GLC
                              Default: true
     -help [or -h]            Print usage to STDOUT.
     -ignore_ic_date          Ignore the date on the initial condition files
                              when determining what input initial condition file to use.
     -ignore_ic_year          Ignore just the year part of the date on the initial condition files
                                 when determining what input initial condition file to use.
     -infile "filepath"       Specify a file (or list of files) containing namelists to
                              read values from.

                              If used with a CLM build with multiple ensembles (ninst_lnd>1)
                              and the filename entered is a directory to files of the
                              form filepath/filepath and filepath/filepath_$n where $n
                              is the ensemble member number. the "filepath/filepath"
                              input namelist file is the master input namelist file
                              that is applied to ALL ensemble members.

                              (by default for CESM this is setup for files of the
                               form $CASEDIR/user_nl_clm/user_nl_clm_????)
     -inputdata "filepath"    Writes out a list containing pathnames for required input datasets in

                                 file specified.
     -irrig "value"           If .true. turn irrigation on with namelist logical irrigate (for |version| physics)
                              (requires crop to be on in the clm configuration)
                              Seek surface datasets with irrigation turned on.  (for CLM4.0 physics)
                              Default: .false.
     -l_ncpl "LND_NCPL"       Number of CLM coupling time-steps in a day.
     -lnd_frac "domainfile"   Land fraction file (the input domain file)
     -mask "landmask"         Type of land-mask (default, navy, gx3v5, gx1v5 etc.)
                              "-mask list" to list valid land masks.
     -namelist "namelist"     Specify namelist settings directly on the commandline by supplying
                              a string containing FORTRAN namelist syntax, e.g.,
                                 -namelist "&clm_inparm dt=1800 /"
     -no-megan                DO NOT PRODUCE a megan_emis_nl namelist for testing that will go into the
                              "drv_flds_in" file for the driver to pass VOCs to the atm.
                              MEGAN (Model of Emissions of Gases and Aerosols from Nature)
                              This removes setting default values for testing MEGAN fields
                              Note: Can always add megan fields to user_nl_clm even with --no-megan
                              (Note: buildnml copies the file for use by the driver)
     -[no-]note               Add note to output namelist  [do NOT add note] about the
                              arguments to build-namelist.
     -rcp "value"             Representative concentration pathway (rcp) to use for
                              future scenarios.
                              "-rcp list" to list valid rcp settings.
     -res "resolution"        Specify horizontal grid.  Use nlatxnlon for spectral grids;
                              dlatxdlon for fv grids (dlat and dlon are the grid cell size
    			      in degrees for latitude and longitude respectively)
                              "-res list" to list valid resolutions.
     -s                       Turns on silent mode - only fatal messages issued.
     -sim_year "year"         Year to simulate for input datasets
                              (i.e. 1850, 2000, 1850-2000, 1850-2100)
                              "-sim_year list" to list valid simulation years
     -bgc_spinup "on|off"     CLM 4.5 Only. For CLM 4.0, spinup is controlled from configure.
                              Turn on given spinup mode for BGC setting of CN
                              on : Turn on Accelerated Decomposition (spinup_state = 1)
                              off : run in normal mode (spinup_state = 0)

                              Default is off.

                              Spinup is now a two step procedure. First, run the model
                              with spinup = "on". Then run the model for a while with
                              spinup = "off". The exit spinup step happens automatically
                              on the first timestep when using a restart file from spinup
                              mode.

                              The spinup state is saved to the restart file.
                              If the values match between the model and the restart
                              file it proceeds as directed.

                              If the restart file is in spinup mode and the model is in
                              normal mode, then it performs the exit spinup step
                              and proceeds in normal mode after that.

                              If the restart file has normal mode and the model is in
                              spinup, then it enters spinup. This is useful if you change
                              a parameter and want to rapidly re-equilibrate without doing
                              a cold start.

     -test                    Enable checking that input datasets exist on local filesystem.
     -verbose [or -v]         Turn on verbose echoing of informational messages.
     -use_case "case"         Specify a use case which will provide default values.
                              "-use_case list" to list valid use-cases.
     -version                 Echo the SVN tag name used to check out this CLM distribution.

     Note: The precedence for setting the values of namelist variables is (highest to lowest):
      1. namelist values set by specific command-line options, like, -d, -sim_year
             (i.e.  CLM_BLDNML_OPTS env_run variable)
      2. values set on the command-line using the -namelist option,
             (i.e. CLM_NAMELIST_OPTS env_run variable)
      3. values read from the file(s) specified by -infile,
             (i.e.  user_nl_clm files)
      4. datasets from the -clm_usr_name option,
             (i.e.  CLM_USRDAT_NAME env_run variable)
      5. values set from a use-case scenario, e.g., -use_case
             (i.e.  CLM_NML_USE_CASE env_run variable)
      6. values from the namelist defaults file.

The ``$CTSMROOT/cime_config/buildnml`` script already sets the resolution and mask as well as the CLM ``configure`` file, and defines an input namelist and namelist input file, and the output namelist directory, and sets the: start-type (from ``RUN_TYPE``); namelist options (from ``CLM_NAMELIST_OPTS``); ``co2_ppmv`` (from ``CCSM_CO2_PPMV``); ``co2_type`` (from ``CLM_CO2_TYPE``); ``lnd_frac`` (from ``LND_DOMAIN_PATH`` and ``LND_DOMAIN_FILE``); ``l_ncpl`` (from ``LND_NCPL``); ``glc_grid``, ``glc_smb``, ``glc_nec`` (from ``GLC_GRID``, ``GLC_SMB``, and ``GLC_NEC``); and ``clm_usr_name`` (to ``CLM_USRDAT_NAME``). Hence only the following different options can be set:

#. ``-bgc_spinup``

#. ``-chk_res``

#. ``-clm_demand``

#. ``-drydep``

#. ``-fire_emis``

#. ``-ignore_ic_date``

#. ``-ignore_ic_year``

#. ``-irrig``

#. ``-no-megan``

#. ``-note``

#. ``-rcp``

#. ``-sim_year``

#. ``-verbose``

``-bgc_spinup`` is an option only available for |version| for any configuration when CN is turned on (so either CLMCN or CLMBGC). It can be set to "on" or "off". If "on" the model will go into Accelerated Decomposition mode, while for "off" (the default) it will have standard decomposition rates. If you are starting up from initial condition files the model will check what mode the initial condition file is in and do the appropriate action on the first time-step to change the Carbon pools to the appropriate spinup setting. See :ref:`spinning-up-clm-bgc` for an example using this option.

.. todo::
    Update the above.

``-chk_res`` ensures that the resolution chosen is supported by CLM. If the resolution is NOT supported it will cause the CLM ``build-namelist`` to abort when run. So when either ``preview_namelist``, ``case.build`` or ``case.run`` is executed it will abort early. Since, the CESM scripts only support certain resolutions anyway, in general this option is NOT needed in the context of running CESM cases.

``-clm_demand`` asks the ``build-namelist`` step to require that the list of variables entered be set. Typically, this is used to require that optional filenames be used and ensure they are set before continuing. For example, you may want to require that fpftdyn be set to get dynamically changing vegetation types. To do this you would do the following.
::

   > ./xmlchange CLM_BLDNML_OPTS="-clm_demand fpftdyn"

To see a list of valid variables that you could set do this:
::

   > cd $CTSMROOT/doc
   > ../bld/build-namelist -clm_demand list

.. note:: Using a 20th-Century transient compset or the ``20thC_transient`` use-case using ``CLM_NML_USE_CASE`` would set this as well, but would also use dynamic nitrogen and aerosol deposition files, so using ``-clm_demand`` would be a way to get *just* dynamic vegetation types and NOT the other files as well.

``-drydep`` adds a dry-deposition namelist for testing to the driver. This is a driver namelist, but adding the option here has CLM ``build-namelist`` create the ``drv_flds_in`` file that the driver will copy over and use. Invoking this option does have an impact on performance even for I compsets and will slow the model down. It's also only useful when running with an active atmosphere model that makes use of this information.

``-ignore_ic_date`` ignores the Initial Conditions (IC) date completely for finding initial condition files to startup from. Without this option or the ``-ignore_ic_year`` option below, the date of the file comes into play.

``-ignore_ic_year`` ignores the Initial Conditions (IC) year for finding initial condition files to startup from. The date is used, but the year is ignored. Without this option or the ``-ignore_ic_date`` option below, the date and year of the file comes into play.

When ``-irrig on`` is used ``build-namelist`` will try to find surface datasets that have the irrigation model enabled (when running with Sattellitte Phenology). When running with the prognostic crop model on, ``-irrig on`` will turn crop irrigation on, while ``-irrig off`` will manage all crop areas as rain-fed without irrigation.

``no-megan`` means do NOT add a MEGAN model Biogenic Volatile Organic Compounds (BVOC) testing namelist to the driver. This namelist is created by default, so normally this WILL be done. This is a driver namelist, so unless ``no-megan`` is specified the CLM ``build-namelist`` will create the ``drv_flds_in`` file that the driver will copy over and use. (If you are running with CAM and CAM produces this file as well, its file will have precedence).

``-note`` adds a note to the bottom of the namelist file, that gives the details of how ``build-namelist`` was called, giving the specific command-line options given to it.

``-rcp`` is used to set the representative concentration pathway for the future scenarios you want the data-sets to simulate conditions for, in the input datasets. To list the valid options do the following:
::

   > cd $CTSMROOT/doc
   > ../bld/build-namelist -rcp list

``-sim_year`` is used to set the simulation year you want the data-sets to simulate conditions for in the input datasets. The simulation ``year`` can also be a range of years in order to do simulations with changes in the dataset values as the simulation progresses. To list the valid options do the following:
::

   > cd $CTSMROOT/doc
   > ../bld/build-namelist -sim_year list

``CLM_NAMELIST_OPTS``
  passes namelist items into one of the CLM namelists. (See :ref:`precedence-of-opts` for the precedence of this option relative to the others.)

  .. note:: For character namelist items you need to use "&apos;" as quotes for strings so that the scripts don't get confused with other quotes they use.

  Example, you want to set ``hist_dov2xy`` to ``.false.`` so that you get vector output to your history files. To do so edit ``env_run.xml`` and add a setting for ``hist_dov2xy``. So do the following:
  ::

     > ./xmlchange CLM_NAMELIST_OPTS="hist_dov2xy=.false."

  Example, you want to set ``hist_fincl1`` to add the variable 'HK' to your history files. To do so edit ``env_run.xml`` and add a setting for ``hist_fincl1``. So do the following:
  ::

     > ./xmlchange CLM_NAMELIST_OPTS="hist_fincl1=&apos;HK&apos;"

  For lists of the history fields available see :ref:`customizing_section`.

``CLM_FORCE_COLDSTART``
   when set to on, *requires* that your simulation do a cold start from arbitrary initial conditions. If this is NOT set, it will use an initial condition file if it can find an appropriate one, and otherwise do a cold start. ``CLM_FORCE_COLDSTART`` is a good way to ensure that you are doing a cold start if that is what you want to do.

.. todo::
    Update the below, as ``queryDefaultNamelist.pl`` no longer exists.

``CLM_USRDAT_NAME``
   Provides a way to enter your own datasets into the namelist. The files you create must be named with specific naming conventions outlined in :ref:`creating-your-own-singlepoint-dataset`. To see what the expected names of the files are, use the ``queryDefaultNamelist.pl`` to see what the names will need to be. For example if your ``CLM_USRDAT_NAME`` will be "1x1_boulderCO", with a "navy" land-mask, constant simulation year range, for 1850, the following will list what your filenames should be:
   ::

      > cd $CTSMROOT/bld
      > queryDefaultNamelist.pl -usrname "1x1_boulderCO" -options mask=navy,sim_year=1850,sim_year_range="constant"  -csmdata $CSMDATA

   An example of using ``CLM_USRDAT_NAME`` for a simulation is given in Example :numref:`creating-your-own-singlepoint-dataset`.

``CLM_CO2_TYPE``
   sets the type of input CO2 for either "constant", "diagnostic" or prognostic". If "constant" the value from ``CCSM_CO2_PPMV`` will be used. If "diagnostic" or "prognostic" the values MUST be sent from the atmosphere model.

===============
 User Namelist
===============

``CLM_NAMELIST_OPTS`` as described above allows you to set any extra namelist items you would like to appear in your namelist. However, it only allows you a single line to enter namelist items, and strings must be quoted with &apos; which is a bit awkward. If you have a long list of namelist items you want to set (such as a long list of history fields) a convenient way to do it is to add to the ``user_nl_clm`` that is created after the ``case.setup`` command runs. The file needs to be in valid FORTRAN namelist format (with the exception that the namelist name ``&namelist`` and the end of namelist marker ``/`` are excluded. The ``preview_namelist`` or ``case.run`` step will abort if there are syntax errors. All the variable names must be valid and the values must be valid for the datatype and any restrictions for valid values for that variable. Here's an example ``user_nl_clm`` namelist that sets a bunch of history file related items, to create output history files monthly, daily, every six and 1 hours.

----------------------------------
Example: user_nl_clm namelist file
----------------------------------

::

   !----------------------------------------------------------------------------------
   ! Users should add all user specific namelist changes below in the form of
   ! namelist_var = new_namelist_value
   !
   ! EXCEPTIONS:
   ! Set use_cndv           by the compset you use and the CLM_BLDNML_OPTS -dynamic_vegetation setting
   ! Set use_vichydro       by the compset you use and the CLM_BLDNML_OPTS -vichydro           setting
   ! Set use_cn             by the compset you use and CLM_BLDNML_OPTS -bgc  setting
   ! Set use_crop           by the compset you use and CLM_BLDNML_OPTS -crop setting
   ! Set spinup_state       by the CLM_BLDNML_OPTS -bgc_spinup      setting
   ! Set co2_ppmv           with CCSM_CO2_PPMV                      option
   ! Set fatmlndfrc         with LND_DOMAIN_PATH/LND_DOMAIN_FILE    options
   ! Set finidat            with RUN_REFCASE/RUN_REFDATE/RUN_REFTOD options for hybrid or branch cases
   !                        (includes $inst_string for multi-ensemble cases)
   !                        or with CLM_FORCE_COLDSTART to do a cold start
   !                        or set it with an explicit filename here.
   ! Set maxpatch_glc       with GLC_NEC                            option
   ! Set glc_do_dynglacier  with GLC_TWO_WAY_COUPLING               env variable
   !----------------------------------------------------------------------------------
   hist_fincl2    = 'TG','TBOT','FIRE','FIRA','FLDS','FSDS',
                    'FSR','FSA','FGEV','FSH','FGR','TSOI',
		    'ERRSOI','BUILDHEAT','SABV','SABG',
		    'FSDSVD','FSDSND','FSDSVI','FSDSNI',
		    'FSRVD','FSRND','FSRVI','FSRNI',
		    'TSA','FCTR','FCEV','QBOT','RH2M','H2OSOI',
                    'H2OSNO','SOILLIQ','SOILICE',
                    'TSA_U', 'TSA_R',
                    'TREFMNAV_U', 'TREFMNAV_R',
                    'TREFMXAV_U', 'TREFMXAV_R',
                    'TG_U', 'TG_R',
                    'RH2M_U', 'RH2M_R',
                    'QRUNOFF_U', 'QRUNOFF_R',
                    'SoilAlpha_U',
                    'Qanth', 'SWup', 'LWup', 'URBAN_AC', 'URBAN_HEAT'
   hist_fincl3 = 'TG:I', 'FSA:I', 'SWup:I', 'URBAN_AC:I', 'URBAN_HEAT:I',
                 'TG_U:I', 'TG_R:I',
   hist_fincl4 = 'TG', 'FSA', 'SWup', 'URBAN_AC', 'URBAN_HEAT'
   hist_mfilt  = 1, 30,  28, 24
   hist_nhtfrq = 0, -24, -6, -1

**Note:** The comments at the top are some guidance given in the default ``user_nl_clm`` and just give some guidance on how to set variables and use the file.

**Note:** You do NOT need to specify the namelist group that the variables are in because the CLM ``build-namelist`` knows the namelist that specific variable names belong to, and it puts them there.

Obviously, all of this would be difficult to put in the CLM_NAMELIST_OPTS variable, especially having to put &apos; around all the character strings. For more information on the namelist variables being set here and what they mean, see the section on CLM namelists below, as well as the namelist definition that gives details on each variable.

.. _precedence-of-opts:

---------------------
Precedence of Options
---------------------

Note: The precedence for setting the values of namelist variables with the different ``env_build.xml``, ``env_run.xml`` options is (highest to lowest):

1. Namelist values set by specific command-line options, like, ``-d``, ``-sim_year`` (i.e. ``CLM_BLDNML_OPTS`` ``env_build.xml`` variable)

#. Values set on the command-line using the ``-namelist`` option, (i.e. ``CLM_NAMELIST_OPTS`` ``env_run.xml`` variable)

#. Values read from the file specified by ``-infile``, (i.e. ``user_nl_clm`` file)

#. Datasets from the ``-clm_usr_name`` option, (i.e. ``CLM_USRDAT_NAME`` ``env_run.xml`` variable)

#. Values set from a use-case scenario, e.g., ``-use_case`` (i.e. ``CLM_NML_USE_CASE`` ``env_run.xml`` variable)

#. Values from the namelist defaults file.

Thus a setting in ``CLM_BLDNML_OPTS`` will override a setting for the same thing given in a use case with ``CLM_NML_USE_CASE``. Likewise, a setting in ``CLM_NAMELIST_OPTS`` will override a setting in ``user_nl_clm``.

.. _setting-initial-conditions:

------------------------------------
Setting Your Initial Conditions File
------------------------------------

Especially with CLMBGC and CLMCN starting from initial conditions is very important. Even with CLMSP it takes many simulation years to get the model fully spunup. There are a couple different ways to provide an initial condition file.

- :ref:`doing-a-hybrid-sim-for-init-conds`
- :ref:`doing-a-branch-sim-for-init-conds`
- :ref:`providing-finidat-in-usernlclm`
- :ref:`adding-finidat-to-xml`

  **Note:** Your initial condition file MUST agree with the surface dataset you are using to run the simulation. If the two files do NOT agree you will get a run-time about a mis-match in PFT weights, or in the number of PFT's or columns.  To get around this you'll need to add the ``use_init_interp=T`` namelist flag in your namelist so that the initial conditions will be interpolated on startup.**

.. _doing-a-hybrid-sim-for-init-conds:

-------------------------------------------------------
Doing a hybrid simulation to provide initial conditions
-------------------------------------------------------

The first option is to setup a hybrid simulation and give a ``RUN_REFCASE`` and ``RUN_REFDATE`` to specify the reference case simulation name to use. When you setup coupled cases (assuming a CESM checkout), at the standard resolution of "f09" it will already do this for you. For example, if you run a "B1850" compset at "f09_g17_gl4" resolution the following settings will already be done for you.

``./xmlchange RUN_TYPE=hybrid,RUN_REFCASE=b.e20.B1850.f09_g17.pi_control.all.297,RUN_REFDATE=0130-01-01,GET_REFCASE=TRUE``

Setting the ``GET_REFCASE`` option to ``TRUE`` means it will copy the files from the RUN_REFDIR usually under: ``$DIN_LOC_ROOT/cesm2_init/$RUN_REFCASE/$RUN_REFDATE`` directory. Note, that the ``RUN_REFCASE`` and ``RUN_REFDATE`` variables are expanded to get the directory name above. If you do NOT set ``GET_REFCASE`` to ``TRUE`` then you will need to have placed the file in your run directory yourself. In either case, the file is expected to be named: ``$RUN_REFCASE.clm2.r.$RUN_REFDATE-00000.nc`` with the variables expanded of course.

.. _doing-a-branch-sim-for-init-conds:

-------------------------------------------------------
Doing a branch simulation to provide initial conditions
-------------------------------------------------------

The setup for running a branch simulation is essentially the same as for a hybrid. With the exception of setting ``RUN_TYPE`` to branch rather than hybrid. A branch simulation runs the case essentially as restarting from it's place before to exactly reproduce it (but possibly output more or different fields on the history files). While a hybrid simulation allows you to change the configuration or run-time options, as well as use a different code base than the original case that may have fewer fields on it than a full restart file. The ``GET_REFCASE`` option works similarly for a branch case as for a hybrid.

.. _providing-finidat-in-usernlclm:

-------------------------------------------------
Providing a finidat file in your user_nl_clm file
-------------------------------------------------

Setting up a branch or hybrid simulation requires the initial condition file to follow a standard naming convention, and a standard input directory if you use the ``GET_REFCASE`` option. If you want to name your file willy nilly and place it anywhere, you can set it in your ``user_nl_clm`` file as in this example.
::

   finidat    = '/glade/home/$USER/myinitdata/clmi_I1850CN_f09_g17_gl4_0182-01-01.c120329.nc'

Note, if you provide an initial condition file -- you can NOT set ``CLM_FORCE_COLDSTART`` to ``TRUE``.

.. _adding-finidat-to-xml:

-------------------------------------------
 Adding a finidat file to the XML database
-------------------------------------------

Like other datasets, if you want to use a given initial condition file to be used for all (or most of) your cases you'll want to put it in the XML database so it will be used by default. The initial condition files, are resolution dependent, and dependent on the number of PFT's and other variables such as GLC_NEC or if irrigation is on or off. See Chapter 3 for more information on this.

------------------------------------
Other noteworthy configuration items
------------------------------------

For running "I" cases there are several other noteworthy configuration items that you may want to work with. Most of these involve settings for the DATM, but one ``CCSM_CO2_PPMV`` applies to all models. The list of DATM settings is `here <http://esmci.github.io/cime/data_models/data-atm.html>`_. If you are running a B, E, or F case that doesn't use the DATM obviously the DATM_* settings will not be used. All of the settings below are in your ``env_build.xml`` and ``env_run.xml`` files
::

   CCSM_CO2_PPMV
   CCSM_BGC
   DATM_MODE
   DATM_PRESAERO
   DATM_YR_ALIGN
   DATM_YR_START
   DATM_YR_END
   DATM_CPLHIST_CASE

``CCSM_CO2_PPMV``
   Sets the mixing ratio of CO2 in parts per million by volume for ALL CESM components to use. Note that most compsets already set this value to something reasonable. Also note that some compsets may tell the atmosphere model to override this value with either historic or ramped values. If the ``CCSM_BGC`` variable is set to something other than "none" the atmosphere model will determine CO2, and CLM will listen and use what the atmosphere sends it. On the CLM side the namelist item ``co2_type`` tells CLM to use the value sent from the atmosphere rather than a value set on it's own namelist.

``DATM_MODE``
   Sets the mode that the DATM model should run in this determines how data is handled as well as what the source of the data will be. Many of the modes are setup specifically to be used for ocean and/or sea-ice modeling. The modes that are designed for use by CLM are (CLM_QIAN, CLMCRUNCEP, CLMCRUNCEPv7, CLMGSWP3v1 and CLM1PT):
   ::

     CLMCRUNCEP
     CLMCRUNCEPv7
     CLMGSWP3v1
     CLM_QIAN
     CLM1PT
     CPLHISTForcing

``CLMCRUNCEP``
   The standard mode for CLM4.5 of using global atmospheric data that was developed by CRU using NCEP data from 1901 to 2010 (version 4 of this series). See :ref:`clmcruncep-and-its-datm` for more information.

``CLMCRUNCEPv7``
   Version 7 of the CRUNCEP data from 1901 to 2016. See :ref:`clmcruncep-and-its-datm` for more information.

``CLMGSWP3v1``
   GSWP3 version 1 forcing data based on NCEP reanalysis with bias corrections by GSWP3 from 1901 to 2010.

``CLM_QIAN``
   The standard mode for CLM4.0 of using global atmospheric data that was developed by Qian et. al. for CLM using NCEP data from 1948 to 2004. See :ref:`clmqian-and-its-datm` for more information. 

``CLM1PT``
   This is for the special cases where we have single-point tower data for particular sites. Right now we only have data for three urban locations: Mexico City Mexico, Vancouver Canada, and the urban-c alpha site. We also have data for the US-UMB AmeriFlux tower site for University of Michigan Biological Station. See :ref:`clm1pt-and-its-datm` for more information.

``CPLHISTForcing``
   This is for running with atmospheric forcing from a previous CESM simulation. See :ref:`cplhistforcing` for more information.

``DATM_PRESAERO``
  sets the prescribed aerosol mode for the data atmosphere model. The list of valid options include:

  ``clim_1850`` = constant year 1850 conditions

  ``clim_2000`` = constant year 2000 conditions

  ``trans_1850-2000`` = transient 1850 to year 2000 conditions

  ``rcp2.6`` = transient conditions for the rcp=2.6 W/m2 future scenario

  ``rcp4.5`` = transient conditions for the rcp=4.5 W/m2 future scenario

  ``rcp6.0`` = transient conditions for the rcp=6.0 W/m2 future scenario

  ``rcp8.5`` = transient conditions for the rcp=8.5 W/m2 future scenario

  ``pt1_pt1`` = read in single-point or regional datasets

DATM_YR_START
  ``DATM_YR_START`` sets the beginning year to cycle the atmospheric data over for ``CLM_QIAN`` or ``CLMCRUNCEP`` or ``CPLHISTForcing`` modes.

DATM_YR_END
  ``DATM_YR_END`` sets the ending year to cycle the atmospheric data over for ``CLM_QIAN`` or ``CLMCRUNCEP`` or ``CPLHISTForcing`` modes.

DATM_YR_ALIGN
  ``DATM_YR_START`` and ``DATM_YR_END`` determine the range of years to cycle the atmospheric data over, and ``DATM_YR_ALIGN`` determines which year in that range of years the simulation will start with.

DATM_CPLHIST_CASE
  ``DATM_CPLHIST_CASE`` sets the casename to use for the ``CPLHISTForcing`` mode.

-----------------------------
Downloading DATM Forcing Data
-----------------------------

In Chapter One of the `CESM User's Guide <link-to-CESM-UG>`_ there is a section on "Downloading input data". The normal process of setting up cases will use the "scripts/ccsm_utils/Tools/check_input_data" script to retrieve data from the CESM subversion inputdata repository. This is true for the standard `CLM_QIAN` forcing as well.

The `CLMCRUNCEP` data is uploaded into the subversion inputdata repository as well -- but as it is 1.1 Terabytes of data downloading it is problematic (*IT WILL TAKE SEVERAL DAYS TO DOWNLOAD THE ENTIRE DATASET USING SUBVERSION*). Because of its size you may also need to download it onto a separate disk space. We have done that on derecho for example where it resides in ``$ENV{CESMROOT}/lmwg`` while the rest of the input data resides in ``$ENV{CESMDATAROOT}/inputdata``. The data is also already available on: janus, franklin, and hopper. If you download the data, we recommend that you break your download into several chunks, by setting up a case and setting the year range for ``DATM_YR_START`` and ``DATM_YR_END`` in say 20 year sections over 1901 to 2010, and then use ``check_input_data`` to export the data.

The ``CPLHISTForcing`` DATM forcing data is unique -- because it is large compared to the rest of the input data, and we only have a disk copy on derecho. The DATM assumes the path for derecho of ``/glade/p/cesm/shared_outputdata/cases/ccsm4/$DATM_CPLHIST_CASE`` for the data. So you will need to change this path in order to run on any other machine.

--------------------------------------
Customizing via the build script files
--------------------------------------

The final thing that the user may wish to do before ``case.setup`` is run is to edit the build script files which determine the configuration and namelist. The variables in ``env_build.xml`` or ``env_run.xml`` typically mean you will NOT have to edit build script files. But, there are rare instances where it is useful to do so. The build script files are copied to your case directory and are available under ``Buildconf``. The list of build script files you might wish to edit are:

``clm.buildexe.csh``
``$CTSMROOT/cime_config/buildnml``
``datm.buildexe.csh``
``datm.buildnml.csh``

.. _more-info-clm-config-script:

--------------------------------------------
More information on the CLM configure script
--------------------------------------------

The CLM ``configure`` script defines the details of a clm configuration and summarizes it into a ``config_cache.xml`` file. The ``config_cache.xml`` will be placed in your case directory under ``Buildconf/clmconf``. The `config_definition_ctsm.xml <https://github.com/ESCOMP/CTSM/blob/master/bld/config_files/config_definition_ctsm.xml>`_ in ``$CTSMROOT/bld/config_files`` gives a definition of each CLM configuration item, it is viewable in a web-browser. Many of these items are things that you would NOT change, but looking through the list gives you the valid options, and a good description of each.

Help on CLM configure
---------------------

Coupling this with looking at the options to CLM ``configure`` with ``-help`` as below will enable you to understand how to set the different options.
::

   > cd $CTSMROOT/bld
   > configure -help

The output to the above command is as follows:
::

   SYNOPSIS
     configure [options]

     Configure CLM in preparation to be built.
   OPTIONS
     User supplied values are denoted in angle brackets (<>).  Any value that contains
     white-space must be quoted.  Long option names may be supplied with either single
     or double leading dashes.  A consequence of this is that single letter options may
     NOT be bundled.

     -bgc <name>            Build CLM with BGC package [ none | cn | cndv ]
                            (default is none).
     -cache <file>          Name of output cache file (default: config_cache.xml).
     -cachedir <file>       Name of directory where output cache file is written
                            (default: CLM build directory).
     -clm4me <name>         Turn Methane model: [on | off]
                              Requires bgc=cn/cndv (Carbon Nitrogen model)
                            (ONLY valid for |version|!)
     -clm_root <dir>        Root directory of clm source code
                            (default: directory above location of this script)
     -cppdefs <string>      A string of user specified CPP defines.  Appended to
                            Makefile defaults.  e.g. -cppdefs '-DVAR1 -DVAR2'
     -vichydro <name>       Turn VIC hydrologic parameterizations : [on | off] (default is off)
     -crop <name>           Toggle for prognostic crop model. [on | off] (default is off)
                            (can ONLY be turned on when BGC type is CN or CNDV)
     -comp_intf <name>      Component interface to use (default ESMF, currently no other option)
     -defaults <file>       Specify full path to a configuration file which will be used
                            to supply defaults instead of the defaults in bld/config_files.
                            This file is used to specify model configuration parameters only.
                            Parameters relating to the build which are system dependent will
                            be ignored.
     -exlaklayers <name>    Turn on extra lake layers (25 layers instead of 10) [on | off]
                            (ONLY valid for |version|!)
     -help [or -h]          Print usage to STDOUT.
     -nofire                Turn off wildfires for BGC setting of CN
                            (default includes fire for CN)
     -noio                  Turn history output completely off (typically for testing).
     -phys <name>           Value of clm4_0 or |version| (default is clm4_0)
     -silent [or -s]        Turns on silent mode - only fatal messages issued.
     -sitespf_pt <name>     Setup for the given site specific single-point resolution.
     -snicar_frc <name>     Turn on SNICAR radiative forcing calculation. [on | off]
                            (default is off)
     -spinup <name>         CLM 4.0 Only. For CLM 4.5, spinup is controlled from  build-namelist.
                            Turn on given spinup mode for BGC setting of CN		  (level)
                              AD            Turn on Accelerated Decomposition from	      (2)
                                            bare-soil
                              exit          Jump directly from AD spinup to normal mode	      (1)
                              normal        Normal decomposition ("final spinup mode")	      (0)
                                            (default)
                            The recommended sequence is 2-1-0
     -usr_src <dir1>[,<dir2>[,<dir3>[...]]]
                            Directories containing user source code.
     -verbose [or -v]       Turn on verbose echoing of settings made by configure.
     -version               Echo the SVN tag name used to check out this CLM distribution.
     -vsoilc_centbgc <name> Turn on vertical soil Carbon profile, CENTURY model decomposition, \

                            split Nitrification/de-Nitrification into two mineral
                            pools for NO3 and NH4 (requires clm4me Methane model), and
                            eliminate inconsistent duplicate soil hydraulic
                            parameters used in soil biogeochem.
                            (requires either CN or CNDV)
                            (ONLY valid for |version|!)
                            [on,off or colon delimited list of no options] (default off)
                              no-vert     Turn vertical soil Carbon profile off
                              no-cent     Turn CENTURY off
                              no-nitrif   Turn the Nitrification/denitrification off
                            [no-vert,no-cent,no-nitrif,no-vert:no-cent]

We've given details on how to use the options in ``env_build.xml`` and ``env_run.xml`` to interact with the CLM ``configure`` and ``build-namelist`` scripts, as well as giving a good understanding of how these scripts work and the options to them. In the next section we give further details on the CLM namelist. You could customize the namelist for these options after ``case.setup`` is run.
