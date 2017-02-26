.. _quickstart:

============
 Quickstart
============

Running the CLM requires a suite of UNIX utilities and programs and you should make sure you have all of these available before trying to go forward with using it. 
If you are missing one of these you should contact the systems administrator for the machine you wish to run on and make sure they are installed.

List of utilities required for CESM in the "CESM1.2.0 Software/Operating System Prerequisites" section in `http://www.cesm.ucar.edu/models/cesm1.2//cesm/doc/usersguide/book1.html <CLM-URL>`_
- UNIX bash shell (for some of the CLM tools scripts)
- NCL (for some of the offline tools for creating/modifying CLM input datasets see `Chapter 2 <CLM-URL>`_ for more information on NCL)
- Python (optional, needed for PTCLM)
- xsltproc, docbook and docbook utilities (optional, needed to build the Users-Guide)

Before working with CLM4.5 read the QuickStart Guide in the `CESM1.2.0 Scripts User's Guide <CLM-URL>`_. Once you are familiar with how to setup cases for any type of simulation with CESM you will want to direct your attention to the specifics of using CLM.

For some of the details of setting up cases for CLM4.5 read the README and text files available from the "models/lnd/clm/doc" directory (see the "CLM Web pages" section for a link to the list of these files). Here are the important ones that you should be familiar with.

1. `README file <CLM-URL>`_ describing the directory structure.

2. `Quickstart.userdatasets <CLM-URL>`_ file describing how to use your own datasets in the model (also see `the Section called Creating your own single-point/regional surface datasets in Chapter 5 <CLM-URL>`_).

3. `models/lnd/clm/doc/KnownBugs <CLM-URL>`_ file describing known problems in CLM4.5 (that we expect to eventually fix).

4. `models/lnd/clm/doc/KnownLimitationss <CLM-URL>`_ file describing known limitations in CLM4.5 and workarounds that we do NOT expect to fix.

The IMPORTANT_NOTES file talks about important things for users to know about using the model scientifically. It content is given in the next chapter on `"What is scientifically validated and functional in CLM4.5 in CESM1.2.0?" <CLM-URL>`_.

The ChangeLog/ChangeSum talk about advances in different versions of CLM. The content of these files is largely explained in the previous chapter on `"What is new with CLM4.5 in CESM1.2.0 since previous public releases?" <CLM-URL>`_.

Note other directories have README files that explain different components and tools used when running CLM and are useful in understanding how those parts of the model work and should be consulted when using tools in those directories. For more details on configuring and customizing a case with CLM see `Chapter 1 <CLM-URL>`_.

The Quickstart.GUIDE (which can be found in ``models/lnd/clm/doc``) is repeated here.
::

          Quick-Start to Using cpl7 Scripts for clm4_5

	  Assumptions: You want to use yellowstone with clm4_5 BGC
              to do a clm simulation with data atmosphere and the
              latest CRUNCEP atm forcing files and settings. You also want to cycle
              the CRUNCEP atm data between 1901 to 1920 and you want to run at
              0.9x1.25 degree resolution.

	  Process:

	  # Create the case

	  cd scripts

	  ./create_newcase -case <testcase> -mach yellowstone_intel -res f09_g16 -compset  I1850CRUCLM45BGC
	  (./create_newcase -help -- to get help on the script)

	  # Setup the case

   cd <testcase>
   ./xmlchange id1=val1,id2=val2  # to make changes to any settings in the env_*.xml files
   ./cesm_setup
   (./cesm_setup -help -- to get help on the script, this creates the ./<testcase>.run  \ 
   script)

   # Add any namelist changes to the user_nl_* files

   $EDITOR user_nl_*

   # Compile the code

   ./<testcase>.build

   # Submit the run

   ./<testcase>.submit

Information on Compsets:

     "I" compsets are the ones with clm and datm7 without ice and ocean. They 
     specify either CLM4.0 physics or CLM4.5 physics.
     Most of the "I" compsets for CLM4.0 use the CLM_QIAN data with solar following
     the cosine of solar zenith angle, precipitation constant, and other
     variables linear interpolated in time (and with appropriate time-stamps on
     the date).  Useful "I" compsets for CLM4.5 use the CRUNCEP data in place
     of the CLM_QIAN data.

     To list all the compsets use:
     ./create_newcase -list compsets

     Some of the CLM4.5 I compsets are:

     Alias               Description
     1850CRUCLM45        CLM4.5 to simulate year=1850 with CLMN45SP (Satellite Phenology)
     I1850CRUCLM45BGC    CLM4.5 to simulate year=1850 with CLM45BGC biogeophysics model (BGC)
     I20TRCRUCLM45BGC    CLM4.5 with BGC on with transient PFT over 1850-2000

     While some of the CLM4 I compsets are:

     Alias               Description
     ICN                 CLM4.0 to simulate year=2000 with Carbon-Nitrogen BGC model (CN)
     I1850CN             CLM4.0 to simulate year=1850 with Carbon-Nitrogen BGC model (CN)
     I20TRCN             CLM4.0 with CN on with transient PFT over 1850-2000
     IRCP26CN            CLM4.0 with CN on with transient PFT over 1850-2100 for RCP=2.6 scenario
     IRCP45CN            CLM4.0 with CN on with transient PFT over 1850-2100 for RCP=4.5 scenario
     IRCP60CN            CLM4.0 with CN on with transient PFT over 1850-2100 for RCP=6.0 scenario
     IRCP85CN            CLM4.0 with CN on with transient PFT over 1850-2100 for RCP=8.5 scenario

Automatically resubmitting jobs:

   After doing a short simulation that you believe is correct

   ./xmlchange CONTINUE_RUN=TRUE

   # Change RESUBMIT to number greater than 0, and CONTINUE_RUN to TRUE...

   ./<testcase>.submit

