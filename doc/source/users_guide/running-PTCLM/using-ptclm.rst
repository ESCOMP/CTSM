.. _using-ptclm.rst:

**************************
Using PTCLM
**************************

There are three types of options to PTCLM1: required, setup/run-time, and dataset generation options. 
The three required options are the three settings that MUST be specified for PTCLM to work at all. The other settings have default values that will default to something useful. The setup/run-time options control how the simulation will be setup and run. The dataset generation options control the generation of datasets needed when PTCLM is run. Most options use a double dash "--" "longname" such as "--list", but the most common options also have a short-name with a single dash (such as -m instead of --machine).

The required options to PTCLM are: inputdata directory (-d), machine (-m) and site-name (-s). 
Inputdata directory is the directory where you have the CESM inputdata files, you need to have write access to this directory, so if you are running on a machine that you do NOT have write access to the standard inputdata location (such as NCAR cheyenne or LBNL hopper) you need to link the standard files to a location you do have control over. We recommend using the ``scripts/link_dirtree`` tool to do that. "machine" is the scripts name for the machine/compiler you will be using for your case. And finally site-name is the name of the site that you want to run for. Site-name can either be a valid supported dataset name or a Fluxnet site name from the list of sites you are running on (see the --sitegroupname for more information about the site lists).

After PTCLM is run a case directory where you can then setup, build and run your CESM case as normal. 
It also creates a ``README.PTCLM`` in that directory that documents the commandline options to PTCLM that were used to create it.

After "help" the "list" option is one of the most useful options for getting help on using PTCLM. 
This option gives you information about some of the other options to PTCLM. To get a list of the machine, sites, and compsets that can be used for PTCLM use the "--list" option as follows.
::

   > cd scripts/ccsm_utils/Tools/lnd/clm/PTCLM
   > ./PTCLM.py --list

The output to the above command is as follows:
::

   /bin/sh: line 1: PTCLM.py: command not found

Steps in running PTCLM
=========================

1. Setup Inputdata directory with write access (use link_dirtree script)

   You need to setup an inputdata directory where you have write access to it. 
   Normally, for NCAR machines the data is on an inputdata where the user does NOT have write access to it. 
   A way that you can get around this is to use the **link_dirtree** script to create softlinks from the normal location to a location you have write access to. 
   So for example on cheyenne:
   ::

      > setenv CSMDATA $CESMDATAROOT/inputdata
      > setenv MYCSMDATA $HOME/inputdata
      > mkdir $MYCSMDATA
      > cd scripts
      > ./link_dirtree $CSMDATA $MYCSMDATA

   See `the Section called Managing Your Own Data-files in Chapter 3 <CLM-URL>`_ for more information on this.

2. Build the CLM tools
   Next you need to make sure all the CLM FORTRAN tools are built.
   ::

      > cd $CTSMROOT/tools/+|version|/mksurfdata_map
      > gmake
      > gmake clean
      > cd ../../../../../../tools/mapping/gen_domain_files/src
      > ../../../../scripts/ccsm_utils/Machines/configure -mach cheyenne -compiler intel
      > gmake
      > gmake clean

3. Run PTCLM
   Next you actually run PTCLM1 which does the different things listed below:

   a. PTCLM names your case based on your input  
      ::

	 [Prefix_]SiteCode_Compset[_QIAN]

      Where: 
	 ``Prefix`` is from the caseidprefix option (or blank if not used).

	 ``SiteCode`` is the site name you entered with the -s option.

	 ``Compset`` is the compset name you entered with the -c option.

	 ``_QIAN`` is part of the name only if the useQIAN is used.
   
      For example, the casename for the following will be:
      ::

	 > cd scripts
	 > ./PTCLM.py -m cheyenne_intel -s US-UMB -d $MYCSMDATA -c ICRUCLM45BGC --use QIAN "US-UMB_I_2000_CN_QIAN"

   b. PTCLM creates datasets for you
      It will populate $MYCSMDATA with new datasets it creates using the CLM tools.

   c. If a transient compset and PTCLM1 finds a _dynpftdata.txt file
      If you are running a transient compset (such as the "I_1850-2000_CN" compset) AND you there is a file in the PTCLM_sitedata directory under the PTCLM directory called $SITE_dynpftdata.txt it will use this file for the land-use changes. 
      Otherwise it will leave land-use constant, unless you use the pftgrid option so it uses the global dataset for landuse changes. 
      See the Section called Dynamic Land-Use Change Files for use by PTCLM for more information on this. 
      There is a sample transient dataset called US-Ha1_dynpftdata.txt. 
      Transient compsets, are compsets that create transient land-use change and forcing conditions such as: 'I_1850-2000', 'I_1850-2000_CN', 'I_RCP8.5_CN', 'I_RCP6.0_CN', 'I_RCP4.5_CN', or 'I_RCP2.6_CN'.

   d. PTCLM creates a pft-physiology for you
      PTCLM1 will create a local copy of the pft-physiology specific for your site that you could then customize with changes specific for that site.

   e. PTCLM creates a README.PTCLM for you
      PTCLM1 will create a simple text file with the command line for it in a file called README.PTCLM in the case directory it creates for you.

4. Customize, setup, build and run case as normal
   You then customize your case as you would normally. See the Chapter 1 chapter for more information on doing this.

PTCLM options
=========================

Next we discuss the setup and run-time options, dividing them up into setup, initial condition (IC), and run-time options.

Configure options include:

- --compset=MYCOMPSET
- --caseidprefix=MYCASEID
- --cesm_root=BASE_CESM
- --namelist=NAMELIST
- --rmold
- --scratchroot=SCRATCHROOT
- --sitegroupname=SITEGROUP
- --QIAN_tower_yrs
- --useQIAN

``--compset``
  The "-c" option is the most commonly used option after the required options, as it specifies the CESM scripts component set to use with PTCLM1. 
  The default compset is the "ICN" compset with CN on for present day conditions.

``--caseidprefix``
  This option gives a prefix to include in the casename when the case is created, in case you want to customize your casenames a bit. 
  By default, casenames are figured out based on the other options. The argument to this option can either be a name to prefix casenames with and/or a pathname to include. 
  Hence, if you want cases to appear in a specific directory you can give the pathname to that directory with this option.

``--cesm_root``
  This option is for running PTCLM1 with a different root directory to CESM than the version PTCLM exists in. Normally you do NOT need to use this option.

``--namelist``
  This option adds any items given into the CLM user_nl_clm namelist. This allows you to add customizations to the namelist before the clm.buildnml.csh file is created for the case.

``--rmold``
  This option will remove an old case directory of the same name if one exists. Otherwise, if an old case directory already exists and you try to run PTCLM it will return with an error.

``--scratchroot``
  This option is ONLY valid when using one of the generic machines (the -m option). This passed onto **create_newcase** and gives the location where cases will be built and run.

``--sitegroupname``
  In the PTCLM directory there is a subdirectory "PTCLM_sitedata" that contains files with the site, PFT and soil data information for groups of sites. 
  These site groups are all separate ASCII files with the same prefix followed by a "_*data.txt" name. See `the Section called PTCLM Group Site Lists <CLM-URL>`_ for more information on these files. By default we have provided three different valid group names:

EXAMPLE
-------
AmeriFlux

Fluxnet-Canada

The EXAMPLE is the group used by default and ONLY includes the US-UMB site as that is the only site we have data provided for. 
The other two site groups include the site information for all of both the AmeriFlux and Fluxnet-Canada sites. 
You can use the "sitegroupname" option to use one of the other lists, or you can create your own lists using the EXAMPLE file as an example. 
Your list of sites could be real world locations or could be theoretical "virtual" sites given to exercise CLM on differing biomes for example. 
Note, see `the Section called Converting AmeriFlux Data for use by PTCLM <CLM-URL>`_ with permission information to use the US-UMB data.

``--useQIAN``
  This option says to use the standard CLM global Qian T62 atmospheric forcing rather than any tower site forcing data available. Otherwise, PTCLM will try to find tower forcing data for the specific site entered.

``--QIAN_tower_yrs``
  This option is used with the "useQIAN" option to set the years to cycle over for the Qian data. In this case Qian atmospheric forcing will be used, but the simulation will run over the same years that tower site is available for this site.

**IC options include:**

- --coldstart
- --finidat=FINIDAT

The coldstart option says to startup with OUT an initial condition file, while the finidat option explicitly gives the initial condition file to use. Obviously, the coldstart and finidat options can NOT be used together.

``--coldstart``
  This option ensures that a cold-start will be done with arbitrary initial conditions.

``--finidat``
  This option sets the initial condition file to startup the simulation from.

**Run-time options include:**

- --debug
- --run_n=MYRUN_N
- --run_units=MYRUN_UNITS
- --stdurbpt
- --debug

This option tells PTCLM to echo what it would do if it were run, but NOT actually run anything. So it will show you the dataset creation commands it would use. It does however, run **create_newcase**, but then it only displays the **xmlchange** commands and changes that it would do. Also note that if you give the "--rmold" option it won't delete the case directory beforehand. Primarily this is intended for debugging the operation of PTCLM.

``--run_n``
  This option along with run_units is used to set the length for the simulation. "run_n" is the number of units to use. The default run length depends on the site, compset, and configuration.

``--run_units``
  This option is the units of time to use for the length of the simulation. It is used along with "run_n" to set the length of the simulation. The default run length depends on the site, compset, and configuration.

``--stdurbpt``
  This option turns on the "stdurbpt_pd" use-case for CLM_NML_USE_CASE. This option can NOT be used for compsets that set the use-case to something besides present-day.

**The dataset generation options are:**

- --pftgrid
- --soilgrid
- --nopointdata
- --owritesrfaer

The options that with a "grid" suffix all mean to create datasets using the global gridded information rather than using the site specific point data. By default the site specific point data is used. The "nopointdata" and "owritesrfaer" options have to do with file creation.

Because supported single-point datasets already have the data created for them, you MUST use the "nopointdata" and "ndepgrid" options when you are using a supported single-point site. You must use "ndepgrid" even for a compset without CN. You also can NOT use the options: "soilgrid", "pftgrid", "aerdepgrid", or "owritesrfaer".

``--pftgrid``
  This option says to use the PFT values provided on the global dataset rather than using the specific site based values from the PTCLM_sitedata/\*_pftdata.txt file when creating the surface dataset.
  This option must NOT be used when you you are using a site that is a supported single point dataset.

``--soilgrid``
  This option says to use the soil values provided on the global dataset rather than using the specific site based values from the PTCLM_sitedata/\*_soildata.txt file when creating the surface dataset.
  This option must NOT be used when you you are using a site that is a supported single point dataset.

``--nopointdata``
  This option says to NOT create any input datasets -- assume this step has already been done. 
  If datasets weren't already created, your case will fail when you try to run it. 
  In general the first time you run PTCLM for a new site you want it to generate new datasets, but the next time and future times you want to use this option so that it doesn't waste a lot of time rebuilding datasets over again.

  .. note:: This option is required when you you are using a site that is a supported single point dataset.

``--owritesrfaer``
  This option says to overwrite any surface and/or aerosol deposition datasets that were already created. 
  Otherwise, the creation of these files will be skipped if a file is already found (but it WILL create files if they don't exist).
  This option must NOT be used when you you are using a site that is a supported single point dataset.
