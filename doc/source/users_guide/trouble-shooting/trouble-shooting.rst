.. _trouble-shooting:

*********************
Trouble Shooting 
*********************

In this chapter we give some guidance on what to do when you encounter some of the most common problems. We can't cover all the problems that a user could potentially have, but we will try to help you recognize some of the most common situations. And we'll give you some suggestions on how to approach the problem to come up with a solution.

In general you will run into one of three type of problems:

1. *setup-time*
#. *build-time*
#. *run-time*


Setup Problems
==============

The first type of problem happens when you invoke the **case.setup** command. 
This indicates there is something wrong with your input datasets, or the details of what you are trying to setup the model to do. 
There's also a trouble-shooting chapter in the `|cesmrelease| Scripts User's Guide <CLM-URL>`_. 
Many of the problems with configuration can be resolved with the guidelines given there. 
Here we will restrict ourselves to problems from the input files.

Example: Missing datasets
----------------------------------------------------------------
::

   > ./create_newcase -case ne60rcp6 -res ne60_g16 -compset IRCP60CN \
   -mach cheyenne_intel
   > ./case.setup

The following is what is displayed to the screen.
::

   .
   .
   .
   Running preview_namelist script 
   CLM configure done.
   CLM adding use_case 1850-2100_rcp6_transient defaults for var clm_demand with val fpftdyn 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var clm_start_type with val startup 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var model_year_align_ndep with val 1850 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var rcp with val 6 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var sim_year with val 1850 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var sim_year_range with val 1850-2100 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var stream_year_first_ndep with val 1850 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var stream_year_last_ndep with val 2100 
   CLM adding use_case 1850-2100_rcp6_transient defaults for var use_case_desc with val Simulate transient land-use, aerosol and Nitrogen deposition changes
   with historical data from 1850 to 2005 and then with the RCP6 scenario from AIM

   build-namelist - No default value found for fpftdyn.
	       Are defaults provided for this resolution and land mask?
   ERROR: clm.buildnml.csh failed
   ERROR: /Users/erik/clm_cesm1_1_1_rel/scripts/ne60rcp6/preview_namelists failed: 25344

The important thing to note here is the line:
::

   ERROR: clm.buildnml.csh failed

which tells us that the problem is in the land **clm.buildnml.csh**. It may also indicate problems in one of the other buildnml.csh files (atm, cesm, cpl, glc, ice, or ocn), in which case you should consult the appropriate model user's guide.

In the example, the error is that the CLM XML database does NOT have a ``finidat`` for the given resolution, rcp scenario and ocean mask. That means you will need to create the file and then supply the file into your case. See `Chapter 2 <CLM-URL>`_ for more information on creating files, and see `Chapter 3 <CLM-URL>`_ for more information on adding files to the XML database. Alternatively, you can provide the file to your case by creating a user namelist as shown in `the Section called User Namelist in Chapter 1 <CLM-URL>`_.

.. note:: The two most common problems from your **clm.buildnml.csh** will be errors from the CLM **configure** or **build-namelist**. For more information on these scripts see: `the Section called More information on the CLM configure script in Chapter 1 <CLM-URL>`_ and `the section on CLM_BLDNML_OPTS <CLM-URL>`_.


Build problems
================

The following is an example of running the build for a case and having it fail in the land model build. 
As you can see it lists which model component is being built and the build log for that component.
::

    CCSM BUILDEXE SCRIPT STARTING
    - Build Libraries: mct pio csm_share 
      Sat Jun 19 21:21:19 MDT 2010 /ptmp/erik/test_build/mct/mct.bldlog.100619-212107
      Sat Jun 19 21:22:18 MDT 2010 /ptmp/erik/test_build/pio/pio.bldlog.100619-212107
      Sat Jun 19 21:23:18 MDT 2010
      /ptmp/erik/test_build/csm_share/csm_share.bldlog.100619-212107
      Sat Jun 19 21:24:00 MDT 2010 /ptmp/erik/test_build/run/cpl.bldlog.100619-212107
      Sat Jun 19 21:24:00 MDT 2010 /ptmp/erik/test_build/run/atm.bldlog.100619-212107
      Sat Jun 19 21:24:06 MDT 2010 /ptmp/erik/test_build/run/lnd.bldlog.100619-212107
      ERROR: clm.buildexe.csh failed, see /ptmp/erik/test_build/run/lnd.bldlog.100619-212107
      ERROR: cat /ptmp/erik/test_build/run/lnd.bldlog.100619-212107

You can then examine the build log that failed and see what went wrong. Most compilers will give the full filepath and line number for the file that filed to compile.

Run Time Problems
=================

Tracking down problems while the model is running is much more difficult to do than setup or build problems. 
In this section we will give some suggestions on how to find run time problems. 
Below we show the log file results of a job that aborted while running.
::

    CCSM PRESTAGE SCRIPT HAS FINISHED SUCCESSFULLY
    Sun Jun 20 18:24:06 MDT 2010 -- CSM EXECUTION BEGINS HERE
    Sun Jun 20 18:24:35 MDT 2010 -- CSM EXECUTION HAS FINISHED
    Model did not complete - see /ptmp/erik/test_run/run/cpl.log.100620-182358

In the next section we will talk about using the different log files to track down problems, and find out where the problem is coming from. In the section after that we give some general advice on debugging problems and some suggestions on ideas that may be helpful to track the problem down. Some of the examples below are from the `$CTSMROOT/doc/KnownBugs <CLM-URL>`_ file.

Tracking Problems by Querying Log Files
---------------------------------------

The first thing to do when tracking down problems is to query the different log files to see if you can discover where the problem occurs, and any error messages about it. 
It's important to figure out if the problem comes in at initialization or in the run phase of the model, and in which model component the problem happens. 
There are different log files for the different major components, and they all end with the date and time in YYMMDD-HHMMSS format (2-digit: year, month, day, hour minute and second). 
When the model runs to completion the log files will be copied to the logs directory in the script directory, but when the model fails they will remain in the run directory. 
Here's an example list of log files from an "I" case where the model dies in the land model initialization. 
For "I" cases the sea-ice and ocean components are just stubs and don't create log files (and unless running with the active land-ice model "glc" log files won't be created either).
::

   atm.log.100620-182358
   cesm.log.100620-182358
   cpl.log.100620-182358
   lnd.log.100620-182358

The coupler log file
--------------------

The first log file to check is the coupler log file so that you can see where the model dies and which model component it fails in. When the model dies at initialization the last model component listed is the component that failed.

Example of a case that fails in the CLM land model initialization.
::

   (seq_timemgr_clockPrint)     Prev Time   = 00001201   00000
   (seq_timemgr_clockPrint)     Next Time   = 99991201   00000
   (seq_timemgr_clockPrint)     Intervl yms =     9999       0           0
   
   (seq_mct_drv) : Initialize each component: atm, lnd, ocn, and ice
   (seq_mct_drv) : Initialize atm component
   (seq_mct_drv) : Initialize lnd component

The cesm log file
-----------------

The cesm log files are to some extent the "garbage collection" of log output. 
The CLM sends it's output from it's master processor, but sends other output and possibly errors to the cesm log file. 
Because, of this, often error messages are somewhere in the cesm log file. 
However, since there is so much other output it may be difficult to find. 
For example, here is some output from an older version of CESM (CESM1.0.2) where the RTM river routing file (before it was converted to NetCDF) was not provided and the error on the open statement for the file was embedded near the end of the cesm log file.
::

   NODE#  NAME
   (    0)  be1105en.ucar.edu
   "/gpfs/proj2/fis/cgd/home/erik/clm_trunk/$CTSMROOT/src/riverroute/RtmMod.F90", line
   239: 1525-155 The file name provided in the OPEN statement for unit 1 has zero length or
   contains all blanks.  The program will recover by ignoring the OPEN statement.
   "/gpfs/proj2/fis/cgd/home/erik/clm_trunk/$CTSMROOT/src/riverroute/RtmMod.F90", line
   241: 1525-001 The READ statement on the file fort.1 cannot be completed because the end
   of the file was reached.  The program will stop.

   Running: ./cesm.exe 
   Please wait...

   Memory usage for   ./cesm.exe (task #   0) is:      51696 KB. Exit status: 1. Signal: 0

Although the example is from an earlier version of the model it still serves to illustrate finding problems from the cesm log file.

When working with the cesm log file, for a run-time problem, you will need to be able to separate it's output into three categories: pre-crash, crash, and post-crash. 
The pre-crash section is everything that is normal output for good operation of the model. 
The crash section is the section where the model dies and reports on the actual problem. 
the post-crash section is the cleanup and finalization after the model dies. 
The most important part of this of course is the crash section. 
The tricky part is distinguishing it from the other sections. 
Also because the cesm log file most likely has duplicated output from multiple processors it is even more difficult to distinguish the different sections and to some extent the sections may be intertwined, as different processors reach the different sections at different times. 
Because, of this reducing the number of processors for your simulation may help you sort out the output in the file (see `the Section called Run with a smaller set of processors <CLM-URL>`_). 
Also much of the output from the cesm log file are system level information having to do with MPI multiprocessing. 
Usually you can ignore this information, but it makes it more difficult to trudge through.


Sometimes the cesm log file is the ONLY file available, because the model terminates early in initialization. 
In this case understanding the output in the cesm log file becomes even more important. 
This also indicates the model did NOT advance far enough to reach the initialization of the individual model components. 
This may mean that the initialization of the multiprocessing for MPI and/or OpenMP failed, or that the reading of the driver namelist file "drv_in" failed.


Here we show those three sections for a cesm log file where a two task job failed on reading the namelist file. 
For a typical job with many tasks similar sections of this will be repeated not just twice but for each task and hence make it harder to read.


*Pre-crash section of the cesm log file*
::

   ATTENTION: 0031-386  MP_INSTANCES setting ignored when LoadLeveler is not being used.

   ATTENTION: 0031-386  MP_INSTANCES setting ignored when LoadLeveler is not being used.
   ATTENTION: 0031-378 MP_EUIDEVICE setting ignored when LoadLeveler is not being used.  
   ATTENTION: 0031-386  MP_INSTANCES setting ignored when LoadLeveler is not being used.
      0:INFO: 0031-724  Executing program: </usr/local/lsf/7.0/aix5-64/bin/lsnrt_run>
      1:INFO: 0031-724  Executing program: </usr/local/lsf/7.0/aix5-64/bin/lsnrt_run>
      0:/contrib/bin/cesm_launch: process 401894 bound to logical CPU 0 on host be0310en.ucar.edu ...
      1:/contrib/bin/cesm_launch: process 439264 bound to logical CPU 1 on host be0310en.ucar.edu ...
      0:INFO: 0031-619  64bit(us, Packet striping on)  ppe_rmas MPCI_MSG: MPI/MPCI library was compiled on   Wed Aug  5 13:36:06 2009
      0: 
      1:LAPI version #14.26 2008/11/23 11:02:30 1.296 src/rsct/lapi/lapi.c, lapi, rsct_rpt53, rpt53s004a 09/04/29 64bit(us)  library compiled on Wed Apr 29 15:30:42 2009
      1:.
      1:LAPI is using lightweight lock.
      0:LAPI version #14.26 2008/11/23 11:02:30 1.296 src/rsct/lapi/lapi.c, lapi, rsct_rpt53, rpt53s004a 09/04/29 64bit(us)  library compiled on Wed Apr 29 15:30:42 2009
      0:.
      0:LAPI is using lightweight lock.
      0:Use health ping for failover/recovery
      1:Use health ping for failover/recovery
      0:Initial communication over instance 2.
      1:Initial communication over instance 0.
      1:IB RDMA initialization completed successfully
      1:The MPI shared memory protocol is used for the job
      0:IB RDMA initialization completed successfully
      0:LAPI job ID for this job is: 1684890719
      0:The MPI shared memory protocol is used for the job
      0:(seq_comm_setcomm)  initialize ID (  7 GLOBAL ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_setcomm)  initialize ID (  2   ATM  ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_setcomm)  initialize ID (  1   LND  ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_setcomm)  initialize ID (  4   ICE  ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_setcomm)  initialize ID (  5   GLC  ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_setcomm)  initialize ID (  3   OCN  ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_setcomm)  initialize ID (  6   CPL  ) pelist   =     0     1     1 ( npes =     2) ( nthreads =  1)
      0:(seq_comm_joincomm) initialize ID (  8 CPLATM ) join IDs =     6     2       ( npes =     2) ( nthreads =  1)
      0:(seq_comm_joincomm) initialize ID (  9 CPLLND ) join IDs =     6     1       ( npes =     2) ( nthreads =  1)
      0:(seq_comm_joincomm) initialize ID ( 10 CPLICE ) join IDs =     6     4       ( npes =     2) ( nthreads =  1)
      0:(seq_comm_joincomm) initialize ID ( 11 CPLOCN ) join IDs =     6     3       ( npes =     2) ( nthreads =  1)
      0:(seq_comm_joincomm) initialize ID ( 12 CPLGLC ) join IDs =     6     5       ( npes =     2) ( nthreads =  1)
      0:  
      0: (seq_comm_printcomms) ID layout : global pes vs local pe for each ID
      0:     gpe        LND      ATM      OCN      ICE      GLC      CPL    GLOBAL   CPLATM   CPLLND   CPLICE   CPLOCN   CPLGLC    nthrds
      0:     ---      ------   ------   ------   ------   ------   ------   ------   ------   ------   ------   ------   ------    ------
      0:       0   :       0        0        0        0        0        0        0        0        0        0        0        0        1
      1:       1   :       1        1        1        1        1        1        1        1        1        1        1        1        1
      1:  
      0: (t_initf) Read in prof_inparm namelist from: drv_in
      1: (seq_io_init) cpl_io_stride, iotasks or root out of bounds - resetting to defaults  4 0 1
      0: piolib_mod.f90 1353 1 2 1 2
      1: piolib_mod.f90 1353 1 2 1 2
      0: pio_support::pio_die:: myrank= 0 : ERROR: piolib_mod.f90: 1354 : not enough procs for the stride
      1: pio_support::pio_die:: myrank= 1 : ERROR: piolib_mod.f90: 1354 : not enough procs for the stride

*Crash section of the cesm log file*
::

   0:
   0:  Traceback:
   1:
   1:  Traceback:
   0:    Offset 0x00000c4c in procedure __pio_support_NMOD_piodie, near line 88 in file pio_support.F90.in
   1:    Offset 0x00000c4c in procedure __pio_support_NMOD_piodie, near line 88 in file pio_support.F90.in
   0:    Offset 0x00000fd0 in procedure __piolib_mod_NMOD_init, near line 1354 in file piolib_mod.F90
   1:    Offset 0x00000fd0 in procedure __piolib_mod_NMOD_init, near line 1354 in file piolib_mod.F90
   1:    Offset 0x00000398 in procedure __seq_io_mod_NMOD_seq_io_init, near line 247 in file /gpfs/proj2/fis/cgd/home/erik/clm_trunk/models/drv/shr/seq_io_mod.F90
   0:    Offset 0x00000398 in procedure __seq_io_mod_NMOD_seq_io_init, near line 247 in file /gpfs/proj2/fis/cgd/home/erik/clm_trunk/models/drv/shr/seq_io_mod.F90
   0:    Offset 0x0001aa88 in procedure ccsm_driver, near line 465 in file /gpfs/proj2/fis/cgd/home/erik/clm_trunk/models/drv/driver/ccsm_driver.F90
   0:    --- End of call chain ---
   1:    Offset 0x0001aa88 in procedure ccsm_driver, near line 465 in file /gpfs/proj2/fis/cgd/home/erik/clm_trunk/models/drv/driver/ccsm_driver.F90
   1:    --- End of call chain ---

*Post-crash section of the cesm log file*
::

      1:Communication statistics of task 1 is associated with task key: 1684890719_1
      0:Communication statistics of task 0 is associated with task key: 1684890719_0
      0:
      0:Running: ./cesm.exe 
      0:Please wait...
      0:
      0:Memory usage for   ./cesm.exe (task #   0) is:     198892 KB. Exit status: 134. Signal: 0
      1:
      1:Running: ./cesm.exe 
      1:Please wait...
      1:
      1:Memory usage for   ./cesm.exe (task #   0) is:     198572 KB. Exit status: 134. Signal: 0
      INFO: 0031-656  I/O file STDOUT closed by task 0
      INFO: 0031-656  I/O file STDERR closed by task 0
      ERROR: 0031-250  task 0: IOT/Abort trap
      INFO: 0031-656  I/O file STDOUT closed by task 1
      INFO: 0031-656  I/O file STDERR closed by task 1
      ERROR: 0031-250  task 1: IOT/Abort trap
      INFO: 0031-639  Exit status from pm_respond = 0
      ATTENTION: 0031-386  MP_INSTANCES setting ignored when LoadLeveler is not being used.
      Job  /usr/local/lsf/7.0/aix5-64/bin/poejob /contrib/bin/ccsm_launch /contrib/bin/job_memusage.exe ./cesm.exe

      TID   HOST_NAME   COMMAND_LINE            STATUS            TERMINATION_TIME
      ===== ========== ================  =======================  ===================
      00000 be0310en   /contrib/bin/ccs  Exit (134)               08/31/2010 12:32:57
      00001 be0310en   /contrib/bin/ccs  Exit (134)               08/31/2010 12:32:57

The CLM log file
----------------

Of course when you are working with and making changes to CLM, most of your focus will be on the CLM log file and the errors it shows. 
As already pointed out if you don't see errors in the ``lnd.log.*`` file you should look in the ``cesm.log.*`` to see if any errors showed up there.

Here's an example of the ``lnd.log.*`` file when running ``PTS_MODE`` with initial conditions (this is bug 1025 in the `$CTSMROOT/doc/KnownLimitationss <CLM-URL>`_ file).
::

 Successfully initialized variables for accumulation
 
 reading restart file I2000CN_f09_g17_gl4_c100503.clm2.r.0001-01-01-00000.nc                                                                                                                                                                                                              
 Reading restart dataset
 ERROR - setlatlon.F:Cant get variable dim for lat or lsmlat
 ENDRUN: called without a message string

The DATM log file
-----------------

When working with "I cases" the second most common problems after CLM problems are problems with the data atmosphere model. So examining the ``atm.log.*`` is important.

Here's an example of a problem that occurs when the wrong prescribed aerosol file is given to a ``pt1_pt1`` simulation.
::

   (datm_comp_init)  atm mode = CLMNCEP
   (shr_strdata_init)  calling shr_dmodel_mapSet for fill
   (shr_strdata_init)  calling shr_dmodel_mapSet for remap
   ('shr_map_getWts') ERROR: yd outside bounds  19.5000000000000000
   (shr_sys_abort) ERROR: ('shr_map_getWts')  ERROR yd outside 90 degree bounds
   (shr_sys_abort) WARNING: calling shr_mpi_abort() and stopping

The batch log files
-------------------

The names of the batch log files will depend on the batch system of the machine that is being used. They will normally be in the script directory. Usually, they don't contain important information, but they are a last resort place to look for error messages. On the NCAR system "cheyenne" the batch files are called with names that start with the batch submission script and then either "stderr.o" or "stdout.o", with the job number at the end.

General Advice on Debugging Run time Problems
=============================================

Here are some suggestions on how to track down a problem while running. In general if the problem still occurs for a simpler case, it will be easier to track down.

1. *Run in DEBUG mode*
#. *Run with a smaller set of processors*
#. *Run in serial mode with a single processor*
#. *Run at a lower resolution*
#. *Run a simpler case*
#. *Run with a debugger*

Run in DEBUG mode
-----------------

The first thing to try is to run in DEBUG mode so that float point trapping will be triggered as well as array bounds checking and other things the compiler can turn on to help you find problems. 
To do this edit the ``env_build.xml`` file and set DEBUG to TRUE as follows:
::

   > ./xmlchange DEBUG=TRUE


Run with a smaller set of processors
------------------------------------

Another way to simplify the system is to run with a smaller set of processors. You will need to clean the setup and edit the --env_mach_pes.xml--. For example, to run with four processors:
::

   > ./case.setup -clean
   > ./xmlchange NTASKS_ATM=4,NTASKS_LND=4,NTASKS_ICE=4,NTASKS_OCN=4,NTASKS_CPL=4,NTASKS_GLC=4
   > ./case.setup

Another recommended simplification is to run without threading, so set the NTHRDS for each component to "1" if it isn't already. Sometimes, multiprocessing problems require a certain number of processors before they occur so you may not be able to debug the problem without enough processors. But, it's always good to reduce it to as low a number as possible to make it simpler. For threading problems you may have to have threading enabled to find the problem, but you can run with 1, 2, or 3 threads to see what happens.

Run in serial mode with a single processor
------------------------------------------

Simplifying to one processor removes all multi-processing problems and makes the case as simple as possible. If you can enable ``MPILIB=mpi-serial`` you will also be able to run interactively rather than having to submit to a job queue, which sometimes makes it easier to run and debug. If you can use ``MPILIB=mpi-serial`` you can also use threading, but still run interactively in order to use more processors to make it faster if needed.
::

   > ./case.setup -clean
   # Set tasks and threads for each component to 1
   # You could also set threads to something > 1 for speed, but still
   # run interactively if threading isn't an issue.
   
   > ./xmlchange NTASKS_ATM=1,NTHRDS_ATM=1,NTASKS_LND=1,NTHRDS_LND=1,NTASKS_ICE=1,NTHRDS_ICE=1
   > ./xmlchange NTASKS_OCN=1,NTHRDS_OCN=1,NTASKS_CPL=1,NTHRDS_CPL=1,NTASKS_GLC=1,NTHRDS_GLC=1
   # set MPILIB to mpi-serial so that you can run interactively
   > ./xmlchange MPILIB=mpi-serial
   > ./case.setup  
   # Then build your case
   # And finally run, by running the *.run script interactively

Run at a lower resolution
-------------------------

If you can create a new case running at a lower resolution and replicate the problem it may be easier to solve. This of course requires creating a whole new case, and trying out different lower resolutions.

Run a simpler case
------------------

Along the same lines, you might try running a simpler case, trying another compset with a simpler setup and see if you can replicate the problem and then debug from that simpler case. Again, of course you will need to create new cases to do this.

Run with a debugger
-------------------

Another suggestion is to run the model with a debugger such as: **dbx**, **gdb**, or **totalview**. 
Often to run with a debugger you will need to reduce the number of processors as outlined above. 
Some debuggers such as **dbx** will only work with one processor, while more advanced debuggers such as **totalview** can work with both MPI tasks and OMP threads. 
Even simple debuggers though can be used to query core files, to see where the code was at when it died (for example using the **where** in **dbx** for a core file can be very helpful. 
For help in running with a debugger you will need to contact your system administrators for the machine you are running on.

