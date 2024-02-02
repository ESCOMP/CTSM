.. _introduction-to-ptclm.rst:

.. include:: ../substitutions.rst

.. _what-is-ptclm:

=====================
 What is PTCLMmkdata?
=====================

PTCLMmkdata (pronounced Pee-Tee Cee-L-M make data is a Python script to help you set up PoinT CLM simulations.

It runs the CLM tools for you to get datasets set up, and copies them to a location you can use them including the changes needed for a case to use the dataset with namelist and XML changes.

Then you run **create_newcase** and point to the directory so that the namelist and XML changes are automatically applied.

PTCLMmkdata has a simple ASCII text file for storing basic information for your sites.

We also have complete lists for AmeriFlux and Fluxnet-Canada sites, although we only have the meteorology data for one site.

For other sites you will need to obtain the meteorology data and translate it to a format that the CESM datm model can use.

But, even without meteorology data PTCLMmkdata is useful to setup datasets to run with standard ``CLM_QIAN`` data.

The original authors of PTCLMmkdata are: Daniel M. Ricciuto, Dali Wang, Peter E. Thornton, Wilfred M. Post all at Environmental Sciences Division, Oak Ridge National Laboratory (ORNL) and R. Quinn Thomas at Cornell University. It was then modified fairly extensively by Erik Kluzek at NCAR. We want to thank all of these individuals for this contribution to the CESM effort. We also want to thank the folks at University of Michigan Biological Stations (US-UMB) who allowed us to use their Fluxnet station data and import it into our inputdata repository, especially Gil Bohrer the PI on record for this site.

.. _details-of-ptclm:

=======================
 Details of PTCLMmkdata
=======================

To get help on PTCLM2_180611 use the "--help" option as follows.
::

   > cd $CTSMROOT/tools/PTCLM
   > ./PTCLMmkdata --help

The output to the above command is as follows:
::

   Usage: PTCLM.py [options] -d inputdatadir -m machine -s sitename

   Python script to create cases to run single point simulations with tower site data.

   Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit

  Required Options:
    -d CCSM_INPUT, --csmdata=CCSM_INPUT
                        Location of CCSM input data
    -m MYMACHINE, --machine=MYMACHINE
                        Machine, valid CESM script machine (-m list to list valid
                        machines)
    -s MYSITE, --site=MYSITE
                        Site-code to run, FLUXNET code or CLM1PT name (-s list to list
                        valid names)

  Configure and Run Options:
    -c MYCOMPSET, --compset=MYCOMPSET
                        Compset for CCSM simulation (Must be a valid 'I' compset [other
                        than IG compsets], use -c list to list valid compsets)
    --coldstart         Do a coldstart with arbitrary initial conditions
    --caseidprefix=MYCASEID
                        Unique identifier to include as a prefix to the case name
    --cesm_root=BASE_CESM
                        Root CESM directory (top level directory with models and scripts
                        subdirs)
    --debug             Flag to turn on debug mode so won't run, but display what would
                        happen
    --finidat=FINIDAT   Name of finidat initial conditions file to start CLM from
    --list              List all valid: sites, compsets, and machines
    --namelist=NAMELIST
                        List of namelist items to add to CLM namelist (example:
                        --namelist="hist_fincl1='TG',hist_nhtfrq=-1"
    --QIAN_tower_yrs    Use the QIAN forcing data year that correspond to the tower
                        years
    --rmold             Remove the old case directory before starting
    --run_n=MYRUN_N     Number of time units to run simulation
    --run_units=MYRUN_UNITS
                        Time units to run simulation (steps,days,years, etc.)
    --quiet             Print minimul information on what the script is doing
    --sitegroupname=SITEGROUP
                        Name of the group of sites to search for you selected site in
                        (look for prefix group names in the PTCLM_sitedata directory)
    --stdurbpt          If you want to setup for standard urban namelist settings
    --useQIAN           use QIAN input forcing data instead of tower site meterology data
    --verbose           Print out extra information on what the script is doing

  Input data generation options:
    These are options having to do with generation of input datasets.  Note: When
    running for supported CLM1PT single-point datasets you can NOT generate new
    datasets.  For supported CLM1PT single-point datasets, you MUST run with the
    following settings: --nopointdata And you must NOT set any of these: --soilgrid
    --pftgrid --owritesrf

    --nopointdata       Do NOT make point data (use data already created)
    --owritesrf         Overwrite the existing surface datasets if they exist (normally
                        do NOT recreate them)
    --pftgrid           Use pft information from global gridded file (rather than site
                        data)
    --soilgrid          Use soil information from global gridded file (rather than site
                        data)

  Main Script Version Id: $Id: PTCLM.py 47576 2013-05-29 19:11:16Z erik $ Scripts URL: $HeadURL: https://svn-ccsm-models.cgd.ucar.edu/PTCLM/trunk_tags/PTCLM1_130529/PTCLM.py $:

Here we give a simple example of using PTCLMmkdata for a straightforward case of running at the US-UMB Fluxnet site on cheyenne where we already have the meteorology data on the machine. Note, see :ref:`converting-ameriflux-for-ptclmmkdata` for permission information to use this data.

Example 6-1. Example of running PTCLMmkdata for US-UMB on cheyenne
------------------------------------------------------------------
::

   > setenv CSMDATA   $CESMDATAROOT/inputdata
   > setenv MYDATAFILES `pwd`/mydatafiles
   > setenv SITE      US-UMB
   > setenv MYCASE    testPTCLM

   # Next build all of the clm tools you will need
   > cd $CTSMROOT/tools/PTCLM
   > buildtools
   # next run PTCLM (NOTE -- MAKE SURE python IS IN YOUR PATH)
   > cd $CTSMROOT/tools/PTCLM
   # Here we run it using qcmd so that it will be run on a batch node
   > qcmd --  ./PTCLMmkdata --site=$SITE --csmdata=$CSMDATA --mydatadir=$MYDATAFILES >& ptclmrun.log &
   > cd $CIMEROOT/scripts
   > ./create_newcase --user-mods-dir $MYDATAFILES/1x1pt_$SITE --case $MYCASE --res CLM_USRDAT --compset I1PtClm50SpGs
   # Next setup, build and run as normal
   > cd $MYCASE
   > ./case.setup

PTCLMmkdata includes a README file that gives some extra details and a simple example.

.. include:: ../../../../tools/PTCLM/README
   :literal:
