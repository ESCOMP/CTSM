.. include:: ../substitutions.rst

.. _using-ptclm.rst:

**************************
Using PTCLMmkdata
**************************

There are two types of options to PTCLMmkdata: required and optional. The three required options are the three settings that MUST be specified for PTCLMmkdata to work at all. The other settings have default values that will default to something useful. Most options use a double dash "--" "longname" such as "--list", but the most common options also have a short-name with a single dash.

The required options to PTCLMmkdata are: inputdata directory (-d) and site-name (-s). Inputdata directory is the directory where you have the CESM inputdata files. Finally site-name is the name of the site that you want to run for. Site-name is a Fluxnet site name from the list of sites you are running on (see the --sitegroupname for more information about the site lists).

After PTCLMmkdata is run you can run **create_newcase** to setup a case to use the datasets created. It also creates a ``README.PTCLM`` in that directory that documents the commandline options to PTCLMmkdata that were used to create it.

After "help" the "list" option is one of the most useful options for getting help on using PTCLMmkdata. This option gives you information about some of the other options to PTCLMmkdata. To get a list of the sites that can be used for PTCLMmkdata use the "--list" option as follows.
::

   > cd $CTSMROOT/tools/PTCLM
   > ./PTCLMmkdata --list

The output to the above command is as follows:
::

   /bin/sh: line 1: PTCLMmkdata: command not found

Steps in running PTCLMmkdata
============================

1. Build the CLM tools Next you need to make sure all the CLM FORTRAN tools are built.
   ::

      > cd $CTSMROOT/tools/PTCLM
      > ./buildtools
      > gmake clean

2. Run PTCLMmkdata Next you actually run PTCLMmkdata which does the different things listed below:

   a. PTCLMmkdata names your output file directory based on your input
      ::

	 [Prefix_]SiteCode

      Where:
	 ``Prefix`` is from the caseidprefix option (or blank if not used).

	 ``SiteCode`` is the site name you entered with the -s option.

      For example, the casename for the following will be:
      ::

	 > cd scripts
	 > ./PTCLMmkdata -s US-UMB -d $MYCSMDATA

   b. PTCLMmkdata creates datasets for you It will populate $MYCSMDATA with new datasets it creates using the CLM tools.

   c. If a transient compset and PTCLMmkdata finds a _dynpftdata.txt file If you are running a transient compset (such as the "I_1850-2000_CN" compset) AND you there is a file in the PTCLM_sitedata directory under the PTCLMmkdata directory called $SITE_dynpftdata.txt it will use this file for the land-use changes. Otherwise it will leave land-use constant, unless you use the pftgrid option so it uses the global dataset for landuse changes. See the Section called Dynamic Land-Use Change Files for use by PTCLMmkdata for more information on this. There is a sample transient dataset called US-Ha1_dynpftdata.txt. Transient compsets, are compsets that create transient land-use change and forcing conditions such as: 'I_1850-2000', 'I_1850-2000_CN', 'I_RCP8.5_CN', 'I_RCP6.0_CN', 'I_RCP4.5_CN', or 'I_RCP2.6_CN'.

   d. PTCLMmkdata creates a pft-physiology for you PTCLMmkdata will create a local copy of the pft-physiology specific for your site that you could then customize with changes specific for that site.

   e. PTCLMmkdata creates a README.PTCLM for you PTCLMmkdata will create a simple text file with the command line for it in a file called README.PTCLM in the case directory it creates for you.

3. Run create_newcase pointing to the directory created

4. Customize, setup, build and run case as normal You then customize your case as you would normally. See the Chapter 1 chapter for more information on doing this.

PTCLMmkdata options
=========================

Next we discuss the setup and run-time options, dividing them up into setup, initial condition (IC), and run-time options.

Configure options include:

- --cesm_root=BASE_CESM
- --sitegroupname=SITEGROUP
- --donot_use_tower_yrs

``--cesm_root``
  This option is for running PTCLMmkdata with a different root directory to CESM than the version PTCLMmkdata exists in. Normally you do NOT need to use this option.

``--sitegroupname``
  In the PTCLMmkdata directory there is a subdirectory "PTCLM_sitedata" that contains files with the site, PFT and soil data information for groups of sites. These site groups are all separate ASCII files with the same prefix followed by a "_*data.txt" name. See :ref:`adding-ptclm-site-data` for more information on these files. By default we have provided three different valid group names:

EXAMPLE
-------
AmeriFlux

Fluxnet-Canada

The EXAMPLE is the group used by default and ONLY includes the US-UMB site as that is the only site we have data provided for. The other two site groups include the site information for all of both the AmeriFlux and Fluxnet-Canada sites. You can use the "sitegroupname" option to use one of the other lists, or you can create your own lists using the EXAMPLE file as an example. Your list of sites could be real world locations or could be theoretical "virtual" sites given to exercise CLM on differing biomes for example. Note, see :ref:`converting-ameriflux-for-ptclmmkdata` with permission information to use the US-UMB data.

``--donot_use_tower_yrs``
  This option is used with the "useQIAN" option to set the years to cycle over for the Qian data. In this case Qian atmospheric forcing will be used, but the simulation will run over the same years that tower site is available for this site.

**Run-time options include:**

- --debug

This option tells PTCLMmkdata to echo what it would do if it were run, but NOT actually run anything. So it will show you the dataset creation commands it would use. It does however, run **create_newcase**, but then it only displays the **xmlchange** commands and changes that it would do. Also note that if you give the "--rmold" option it won't delete the case directory beforehand. Primarily this is intended for debugging the operation of PTCLMmkdata.

**The dataset generation options are:**

- --pftgrid
- --soilgrid

The options that with a "grid" suffix all mean to create datasets using the global gridded information rather than using the site specific point data. By default the site specific point data is used. The "nopointdata" and "owritesrfaer" options have to do with file creation.

Because supported single-point datasets already have the data created for them, you MUST use the "nopointdata" and "ndepgrid" options when you are using a supported single-point site. You must use "ndepgrid" even for a compset without CN. You also can NOT use the options: "soilgrid", "pftgrid", "aerdepgrid", or "owritesrfaer".

``--pftgrid``
  This option says to use the PFT values provided on the global dataset rather than using the specific site based values from the PTCLM_sitedata/\*_pftdata.txt file when creating the surface dataset. This option must NOT be used when you you are using a site that is a supported single point dataset.

``--soilgrid``
  This option says to use the soil values provided on the global dataset rather than using the specific site based values from the PTCLM_sitedata/\*_soildata.txt file when creating the surface dataset. This option must NOT be used when you you are using a site that is a supported single point dataset.

