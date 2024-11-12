.. include:: ../substitutions.rst

.. _trouble-shooting:

***************
Troubleshooting
***************

In this chapter we give some guidance on what to do when you encounter some of the most common problems.

In general you may run into one of four types of problems:

1. *case-creation*
#. *setup-time*
#. *build-time*
#. *run-time*

Start with the `CIME Trouble Shooting Guide <http://esmci.github.io/cime/versions/master/html/users_guide/troubleshooting.html>`_ , especially if you encounter one of the first three types of problems. The CIME troubleshooting guide also provides some useful tips regarding run-time errors. If this doesn't identify and solve your problem, then try some of the suggestions below for run-time errors.

General Advice on Debugging Run time Problems
=============================================

The model has been run for thousands and thousands of simulation years in many different configurations, both fully-coupled and in land-only modes, without problems.  If you have modified the model in any way, by using either different input datasets or new or modified code, then that is the first place to look if you encounter an error. However, the model is not completely infallible, as noted below.

It is important to examine all of the component log files in the run directory for errors. An error in the land model may not appear in the lnd log file, it may show up in the cesm log. In a land-only simulation, errors associated with the data atmosphere model may show up in the atm log or the cesm log, or both. The two logs together may contain useful information about the error.  Frequently, the error output in the log files will include a **traceback** of code where the error occurred. Identifying the specific line of code where the error occurred is the first step in diagnosing the error and developing a solution. If a traceback doesn't appear in the log files, then try running in debug mode as noted in the CIME troubleshooting guide. An example of a traceback in the cesm log is given below
::

   398: ERROR: Carbon or Nitrogen patch negative =   -60.0630620423182
   398:  -1.49270707132601
   398: ERROR: limits =   -60.0000000000000       -6.00000000000000
   398: iam = 362: local  patch    index = 482
   398: iam = 362: global patch    index = 163723
   398: iam = 362: global column   index = 104283
   398: iam = 362: global landunit index = 32348
   398: iam = 362: global gridcell index = 13723
   398: iam = 362: gridcell longitude    =  120.0000000
   398: iam = 362: gridcell latitude     =  -70.0000000
   398: iam = 362: pft      type         = 10
   398: iam = 362: column   type         = 1
   398: iam = 362: landunit type         = 1
   398: ENDRUN:
   398: ERROR:
   398: ERROR: carbon or nitrogen state critically negative ERROR in CNPrecisionControl
   398: Mod.F90 at line 209
   398:Image              PC                Routine            Line        Source
   398:cesm.exe           000000000383B3EA  Unknown               Unknown  Unknown
   398:cesm.exe           0000000002F1E5D0  shr_abort_mod_mp_         114  shr_abort_mod.F90
   398:cesm.exe           0000000001AF22BF  abortutils_mp_end          50  abortutils.F90
   398:cesm.exe           0000000001D02677  cnprecisioncontro         693  CNPrecisionControlMod.F90
   398:cesm.exe           0000000001CFCC58  cnprecisioncontro         207  CNPrecisionControlMod.F90
   398:cesm.exe           00000000021FB4F5  cndrivermod_mp_cn         575  CNDriverMod.F90
   398:cesm.exe           0000000001D0F5C7  cnvegetationfacad         866  CNVegetationFacade.F90
   398:cesm.exe           0000000001AFEC96  clm_driver_mp_clm         925  clm_driver.F90
   398:cesm.exe           0000000001AE744B  lnd_comp_mct_mp_l         458  lnd_comp_mct.F90
   398:cesm.exe           0000000000429414  component_mod_mp_         737  component_mod.F90
   398:cesm.exe           000000000040AE4B  cime_comp_mod_mp_        2622  cime_comp_mod.F90
   398:cesm.exe           000000000042904C  MAIN__                    133  cime_driver.F90
   398:cesm.exe           0000000000408D22  Unknown               Unknown  Unknown
   398:libc.so.6          00002B8B95D306E5  __libc_start_main     Unknown  Unknown
   398:cesm.exe           0000000000408C29  Unknown               Unknown  Unknown

Here, the output is identifying the sequence of Fortran statements involved in the error, starting with line 133 in cime_driver.F90 and ending with line 114 in shr_abort_mod.F90.  In this case the run is triggering an error check in the model related to negative carbon/nitrogen at line 693 of CNPrecisionControlMod.F90. In addition, there is additional information related to the error indicating the carbon or nitrogen state is critically negative at line 209 in CNPrecisionControlMod.F90, which is
::

   call TruncateCandNStates( bounds, filter_soilp, num_soilp, cs%leafc_patch(bounds%begp:bounds%endp), &
                             ns%leafn_patch(bounds%begp:bounds%endp), &
                             pc(bounds%begp:), pn(bounds%begp:), __LINE__, &
                             num_truncatep, filter_truncatep)

So here we know that it is either leaf nitrogen (leafn) or leaf carbon (leafc) that has triggered the error.  Furthermore, from the information above and looking at the code, we see that it is leaf carbon that is triggering the error, as the value is -60.0630620423182 and the limit is -60.

At this point it is useful as a next step to identify the particular patch index and perhaps the pft type that is triggering the error. In this case, the endrun call is already written to provide this information: the patch index and pft type causing the error, along with some other information, are printed in the lines beginning with ``iam``. The ``iam`` value gives the CTSM processor number (this can be obtained in the code via the ``iam`` variable defined in ``spmdMod``). The local patch index is the value of ``p`` in the current patch loop; "local" implies that it refers to this processor's indexing. However, this same value of ``p`` may appear on other processors, since the local indexing on each processor starts with 1. So, to get the unique patch causing the problem, you either need to use the processor's ``iam`` index (there is only one patch with local index 482 on processor 362), or use the global indices printed below the local index. The "global" term here refers to the global index space across all processors (there is only one patch with a global index of 163723 across all processors). See below for how to use the ``get_global_index`` function to translate from local to global indices.

If you are writing your own ``endrun`` call, you can get this additional information by specifying the ``subgrid_index`` and ``subgrid_level`` arguments; for example:

::
  
  call endrun(subgrid_index=p, subgrid_level=subgrid_level_patch, msg=errMsg(sourcefile, __LINE__))

(The ``subgrid_level_patch`` constant, and similar constants for the other subgrid levels, are defined in ``decompMod``, so can be accessed via ``use decompMod, only : subgrid_level_patch``.)

You can get this same information without aborting the run via a call to ``write_point_context``, which is also defined in the ``abortutils`` module; e.g.:

::
  
  if (abs(carbon_patch(p)) < ccrit) then
     call write_point_context(subgrid_index=p, subgrid_level=subgrid_level_patch)
  end if

Or, if all you want is the global index of ``p`` for the sake of writing extra diagnostic prints like the example below, then you can use the ``get_global_index`` function defined in ``decompMod``, like:

::
   
   if (abs(carbon_patch(p)) < ccrit) then
      write(iulog,*) 'carbon patch significantly negative at local, global p = ', &
           p, get_global_index(subgrid_index=p, subgrid_level=subgrid_level_patch)
   end if

In all of these cases, the output will appear in either the cesm or lnd log file. In the above example, we see that the local patch index is 482 on processor 362 and the global patch index is 163723. From there, one can use this patch index to write out variables that are used in updating leafc, for example, leafc is updated a number of times in CNCStateUpdate1Mod.F90.

There are two equivalent methods to write a conditional statement to provide more output for the problem patch within a loop over all patches. The first method is to translate the local index to a global index:

::
   
   use decompMod, only : get_global_index, subgrid_level_patch
   ...
   if (get_global_index(p, subgrid_level_patch) == 163723) then
      write(iulog,*)'CNCStateUpdate1Mod leafc: ',cs_veg%leafc_patch(p)
      write(iulog,*)'CNCStateUpdate1Mod +leafc_xfer_to_leafc: ',cf_veg%leafc_xfer_to_leafc_patch(p)*dt
   end if

The second method is to use the local index along with the processor number:

::
   
   use spmdMod, only : iam
   ...
   if (p == 482 .and. iam == 362) then
      write(iulog,*)'CNCStateUpdate1Mod leafc: ',cs_veg%leafc_patch(p)
      write(iulog,*)'CNCStateUpdate1Mod +leafc_xfer_to_leafc: ',cf_veg%leafc_xfer_to_leafc_patch(p)*dt
   end if

By placing these write statements in the code, one can get a sense of how leafc is evolving toward a negative state and why. This is a very complex example of troubleshooting. To make a long story short, as described `here <https://github.com/ESCOMP/CTSM/issues/1163>`_, the error turned out to be caused by a few lines in the phenology code that weren't handling a 20 minute time step properly, thus an actual bug in the code. This was also a good example of where a much less computationally expensive land-only simulation was able to be used for debugging instead of the orginal expensive fully-coupled simulation.

Another method of troubleshooting is to use the ``point_of_interest`` module.

Use the point_of_interest module
--------------------------------

It is common, when debugging, to want to print the values of various variables for all patches or columns of certain landunit types within a certain grid cell of interest. For example, one might be able to identify a certain grid cell with an erroneous value for a particular history field variable (e.g., GPP) using for example ncview. Once the latitude and longitude of this grid cell has been determined, the point_of_interest module (``src/utils/point_of_interest.F90``) helps create the logical functions needed to do this. This module is compiled into every CTSM build, but is not invoked by default. To use it

(1) Enter in the latitude/longitude of the point of interest in the function ``at_poi`` in ``point_of_interest.F90`` by setting the variables ``poi_lat`` and ``poi_lon``.

(2) You may customize the ``point_of_interest.F90`` code by changing the example function (``poi_c``) and/or adding new functions. Look for comments about "Customize" to see what to customize.

(3) Add calls to these functions in the CTSM code

The example function in ``point_of_interest.F90`` is ``poi_c``. It finds columns with a given landunit type (in this case, the natural vegetated landunit). That function can be used in a column-level loop to find columns with that landunit within the grid cell of interest. Its typical use in CTSM code is
::
   
   do fc = 1, num_nolakec
      c = filter_nolakec(fc)
      ! Various code here, maybe setting foo and bar variables
      if (poi_c(c)) then
         write(iulog,*) 'DEBUG: foo, bar = ', foo(c), bar(c)
      end if
   end do

You will also need a ``use`` statement in the module from which you are calling ``poi_c``
::

   use point_of_interest, only : poi_c

Here are some other suggestions on how to track down a problem encountered while running. These involve setting up and running simpler cases.  In general if the problem still occurs for a simpler case, it will be easier to track down. However, we note that most errors are specific to the case being run.

#. *Run with a smaller set of processors*
#. *Run in serial mode with a single processor*
#. *Run at a lower resolution*
#. *Run a simpler case*
#. *Run with a debugger*

Run with a smaller set of processors
------------------------------------

One way to simplify the system is to run with a smaller set of processors. You will need to clean the setup and edit ``env_mach_pes.xml``. For example, to run with four processors:
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

Another suggestion is to run the model with a debugger such as: ``ddt``, ``dbx``, ``gdb``, or ``totalview``. Often to run with a debugger you will need to reduce the number of processors as outlined above. Some debuggers such as ``dbx`` will only work with one processor, while more advanced debuggers such as ``totalview`` can work with both MPI tasks and OMP threads. Even simple debuggers though can be used to query core files, to see where the code was at when it died (for example using the ``where`` in ``dbx`` for a core file can be very helpful. For help in running with a debugger you will need to contact your system administrators for the machine you are running on.

