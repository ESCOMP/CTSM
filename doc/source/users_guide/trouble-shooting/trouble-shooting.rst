.. include:: ../substitutions.rst

.. _trouble-shooting:

*********************************************
Introduction
*********************************************

In this chapter we give some guidance on what to do when you encounter some of the most common problems.

In general you may run into one of four types of problems:

1. *case-creation*
2. *setup-time*
3. *build-time*
4. *run-time*

Start with the `CIME Trouble Shooting Guide <https://esmci.github.io/cime/versions/master/html/ccs/troubleshooting.html>`_ , especially if you encounter one of the first three types of problems. The CIME troubleshooting guide also provides some useful tips regarding run-time errors. If this doesn't identify and solve your problem, then try some of the suggestions below for run-time errors.

*********************************************
General Advice on Debugging Run time Problems
*********************************************

Release versions of the model have been run for thousands and thousands of simulation years in many different configurations, both fully-coupled and in land-only modes, without problems.  If you have modified the model in any way, by using either different input datasets or new or modified code, then that is the first place to look if you encounter an error. In particular, see the following section :numref:`List of common problems to watch out for when developing / reviewing code` for advice on problems to watch out for when developing or reviewing code. 

.. _List of common problems to watch out for when developing / reviewing code:

=====================================================================================================
List of common problems to watch out for when developing / reviewing code
=====================================================================================================

-----------------------------------------------------------------------------------------------------
Possible sources of model crashes or non-physical results
-----------------------------------------------------------------------------------------------------

Potential for divide-by-zero
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Solution: put code in a conditional that checks for 0, and handles 0 values specially

Potential for other floating point exceptions (e.g,. raising 0 to a negative power, etc.)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Solution: put code in a conditional that handles mathematically impossible cases specially

Potential for a quantity that should be non-negative to go negative, due to rounding errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Solutions:
::

  foo = max(foo, 0._r8)

or use the ``truncate_small_values`` subroutine in ``NumericsMod``.

-----------------------------------------------------------------------------------------------------
Possible sources of science bugs and maintainability problems
-----------------------------------------------------------------------------------------------------

Block of code is copy & pasted, or effectively duplicated
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is one of the most common issues we encounter. It causes both maintainability problems and future bugs (when one block of code gets modified and thus becomes out of sync with the other block). Sometimes a block of code is copied and pasted entirely; other times this issue is more subtle, with logic effectively duplicated in two places even though the code looks somewhat different.

This problem can be identified by asking yourself: If I make a change in place A, will I need to make a corresponding change in place B?

Possible solutions:

- Introduce a subroutine that holds the common code, with parameters to allow for any differences in behavior between the different locations

- Save the result of a calculation in a variable that can be reused in multiple places

Variable is missing from the restart file, or is on the restart file when it doesn't need to be
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Variables need to be written to and read from the restart file if their value persists from one time step to the next. This is mainly the case for the model's fundamental state variables. However, in order to minimize restart file size (reasons of disk usage, performance, and understandability), we try to avoid adding variables to the restart file unless they're truly needed.

Some good rules of thumb to use are:

- If a variable's value at time ``t`` directly depends on its value at time ``t-1``, often through evolution equations of the form ``x = x + flux*dtime``, then it likely needs to be on the restart file.

- Imagine setting a variable to 0 at the start of each time step. Would this cause incorrect results? If so, it likely needs to be on the restart file.

- However, if a variable can be recalculated based on other variables, then it probably should *not* be on the restart file. Instead, its initial value can be calculated in model initialization, after reading the restart file.

There are also some cases where a variable is needed on the restart file because its value is referenced (i.e., it appears on the right-hand side of an equation) earlier in the driver loop than where it is set. This is a subtle issue, and needs to be kept in mind when developing and reviewing code. For example, the relevant parts of the driver loop could look like this:
::

  foo = bar*2

  (more code here)

  bar = ...

In reality, these lines will be in subroutines called from the main driver loop, so an understanding is needed of the calling order of the model's subroutines.

An ideal solution in this case is to reorder the code so that ``bar`` is calculated before it is used. The next most ideal solution is to recalculate ``bar`` from other variables in initialization after the restart file is read. However, if neither of these are possible, then ``bar`` needs to be added to the restart file.

We have many restart tests in the automated test suite; these catch many problems with variables being absent from the restart file that should be there. However, these tests cannot catch all problems - particularly if a variable's value only needs to persist between time steps in certain circumstances (such as if the restart is done mid-day). In addition, these tests cannot catch problems with variables being added to the restart file unnecessarily.

Variable is used on right-hand side of an equation before it has its "final" value
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example:
::

   bar = ...

   foo = bar + 1._r8

   (more code here)

   bar = bar + 1._r8

In this case, it's possible that the ``foo`` assignment should really have happened after the increment to ``bar``.

Solution: Check code carefully for assignments to variables that are used on the right-hand side of equations. Ideally, this search would only need to be done within the current subroutine. But in practice, variables in CTSM are sometimes updated in multiple subroutines, so you should extend this search to make sure your new code happens in the correct place in the driver loop. (i.e., make sure that there aren't subroutines called later in the driver that update the quantity that you're using on the right-hand side of the equation.)

------------------------------------------------------------------------------------------------------------
Possible sources of answer changes with changing processor count (PEM test failures, also seen in ERP tests)
------------------------------------------------------------------------------------------------------------

Problems specific to parallelization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

PEM and ERP tests are designed to catch problems specific to parallelization. In CTSM, these aren't the most common sources of errors with these tests, but we'll start with a few parallel-specific reasons that these tests could fail.

Missing or incorrect broadcast for a namelist variable
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Namelist variables are read on the master proc and then should be broadcast to all other processors. If a broadcast statement is missing or incorrect for a namelist variable, then the namelist value could be wrong on all other processors. This will lead to answer changes with changing processor count because changing processor count will change which grid cells are on the master proc (with the correct namelist value) vs. other processors.

Processor count dependence of a parallel algorithm
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

An obvious source of answer changes with changing processor count is processor count dependence of a parallel algorithm. A common issue here is an MPI reduction that depends on the processor count – e.g., a sum across multiple processors, which could depend on the order in which the sum is taken.

However, we do not have many parallel algorithms in CTSM; these would mainly apply in cases where there is communication between grid cells.

Incorrect indexing of a subgrid variable
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One common cause of processor count dependence is the incorrect indexing of a subgrid variable. For example, if a patch variable is indexed by `g` instead of `p`, then it will access the wrong index. This will be picked up in a PEM (or ERP) test because exactly which point it accesses is processor count-dependent.

Here are some `git grep` commands that can help find this problem:

.. code:: shell

   > git grep -i '_patch \*( \*[glc] \*)'
   > git grep -i '_col \*( \*[glp] \*)'
   > git grep -i '_lun \*( \*[gcp] \*)'
   > git grep -i '_grc \*( \*[lcp] \*)'

However, since we often strip the suffix in associate statements, you cannot rely on these grep commands to detect this issue.

Scalar variable used before it is set in a given loop iteration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Another common cause of answer changes with changing processor counts is a scalar variable being used before it is set in a given loop iteration. This means that its value depends on a previous loop iteration, or possibly the value that was set in an earlier loop in this subroutine. Changing processor counts changes which grid cell is operated on first in a loop for a given processor, and also which grid cell is the previous loop iteration for a given grid cell.

There are a few common specific ways that this appears in CTSM code, as noted below:

Missing setting of `c`, `l` or `g` in a loop
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Often, a loop over one subgrid level will access variables in arrays at a coarser subgrid level. For example, a loop over patches will access column and gridcell-level variables. This requires settings like `c = patch%column(p)`. Sometimes there is a bug where a given loop is missing one of these needed settings; instead its setting comes from the previous loop in that subroutine. In this case, all patches – on all grid cells – will use the same `c` value.

Scalar variable set in a conditional but accessed outside that conditional
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Sometimes a scalar variable is set inside a conditional but is accessed outside that conditional. There may be multiple branches of the conditional with the intent that the scalar is set for all cases, but there may be a missing branch, so in some situations the scalar doesn't end up getting set for a particular point. The value will then be taken from the previous loop iteration.

------------------------------------------------------------------------------------------------------------
Possible sources of threading bugs
------------------------------------------------------------------------------------------------------------

Whole-array assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

  foo(:) = 0._r8

should be replaced by the following, assuming ``foo`` is a column-level array:
::

  foo(bounds%begc:bounds%endc) = 0._r8

or, better, initialize ``foo`` within a loop over the appropriate filter, or a loop over bounds%begc to bounds%endc, ideally subset by active points.

Forgetting to index into a variable for assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
E.g., `variable_patch = 0._r8` instead of `variable_patch(p) = 0._r8`.

Improper argument passing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See [this page on the old wiki](https://wiki.ucar.edu/display/ccsm/CLM+Coding+Conventions#CLMCodingConventions-ArgumentPassingargPass)

Identify specific loops with issues by turning off threading for loops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One way to identify threaded loops with issues is to turn off more and more loops until you identify ones with issues.

For example a subroutine with a loop like this...
::

  !$OMP PARALLEL DO PRIVATE (nc,bounds_clump)
  do nc = 1,nclumps
     call get_clump_bounds(nc, bounds_clump)
     ...
  end do
  !$OMP END PARALLEL DO

Incorrect list of private variables in threaded loops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For threading each processor thread needs to have its own version of temporary or loop indexing variables.  For longer loops the list needed can be quite long, if you don't have the right list of variables in the private list, different threads will share these variables and result in strange behavior.

So in the above example, if either nc, or bounds_clump weren't in the private list the loop would not be able to function correctly.

Turn off OpenMP parallelism by removing the first line or adding an extra comment character (!) to comment it out.

-----------------------------------------------------------------------------------------------------
Possible sources of performance problems
-----------------------------------------------------------------------------------------------------

Doing an operation on all array elements rather than just active points
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Active points are (generally, but not entirely) ones with > 0 weight on the grid cell.

Solution: Use the filters, which only include active points

Nested loops in the wrong order for cache-friendliness and/or vectorizability
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-----------------------------------------------------------------------------------------------------
Other things to look for in a PR review
-----------------------------------------------------------------------------------------------------

Binary files committed without the use of git lfs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If someone commits a binary file without git lfs enabled, it will actually be committed directly. The same thing will happen even if they have git lfs enabled if the file has an extension that isn't currently tracked by git lfs. Look for any such binary files when looking through the list of changed files. These will often appear in the doc directory.

See the .gitattributes file at the top level of the repository for files typically handled by git lfs).

.. _A specific troubleshooting example:

=====================================================================================================
A specific troubleshooting example
=====================================================================================================

First of all, it is important to examine all of the component log files in the run directory for errors. An error in the land model may not appear in the lnd log file, it may show up in the cesm log. In a land-only simulation, errors associated with the data atmosphere model may show up in the atm log or the cesm log, or both. The two logs together may contain useful information about the error.  Frequently, the error output in the log files will include a **traceback** of code where the error occurred. Identifying the specific line of code where the error occurred is the first step in diagnosing the error and developing a solution. If a traceback doesn't appear in the log files, then try running in debug mode as noted in the CIME troubleshooting guide. An example of a traceback in the cesm log is given below (this specific example is from cesm2_2_beta05)
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

Here, the output is identifying the sequence of Fortran statements involved in the error, starting with line 133 in cime_driver.F90 and ending with line 114 in shr_abort_mod.F90.  In this case the run is triggering an error check in the model related to negative carbon/nitrogen at line 693 of CNPrecisionControlMod.F90. In addition, there is additional information related to the error indicating the carbon or nitrogen state is critically negative at line 207 in CNPrecisionControlMod.F90, which is a subroutine call
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

Another method of troubleshooting is to :ref:`Use the point_of_interest module`.

.. _Use the point_of_interest module:

=====================================================================================================
Use the point_of_interest module
=====================================================================================================

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

=====================================================================================================
Run with a smaller set of processors
=====================================================================================================

One way to simplify the system is to run with a smaller set of processors. You will need to clean the setup and edit ``env_mach_pes.xml``. For example, to run with four processors:
::

   > ./case.setup -clean
   > ./xmlchange NTASKS_ATM=4,NTASKS_LND=4,NTASKS_ICE=4,NTASKS_OCN=4,NTASKS_CPL=4,NTASKS_GLC=4
   > ./case.setup

Another recommended simplification is to run without threading, so set the NTHRDS for each component to "1" if it isn't already. Sometimes, multiprocessing problems require a certain number of processors before they occur so you may not be able to debug the problem without enough processors. But, it's always good to reduce it to as low a number as possible to make it simpler. For threading problems you may have to have threading enabled to find the problem, but you can run with 1, 2, or 3 threads to see what happens.

=====================================================================================================
Run in serial mode with a single processor
=====================================================================================================

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

=====================================================================================================
Run at a lower resolution
=====================================================================================================

If you can create a new case running at a lower resolution and replicate the problem it may be easier to solve. This of course requires creating a whole new case, and trying out different lower resolutions.

=====================================================================================================
Run a simpler case
=====================================================================================================

Along the same lines, you might try running a simpler case, trying another compset with a simpler setup and see if you can replicate the problem and then debug from that simpler case. Again, of course you will need to create new cases to do this.

=====================================================================================================
Run with a debugger
=====================================================================================================

Another suggestion is to run the model with a debugger such as: ``ddt``, ``dbx``, ``gdb``, or ``totalview``. Often to run with a debugger you will need to reduce the number of processors as outlined above. Some debuggers such as ``dbx`` will only work with one processor, while more advanced debuggers such as ``totalview`` can work with both MPI tasks and OMP threads. Even simple debuggers though can be used to query core files, to see where the code was at when it died (for example using the ``where`` in ``dbx`` for a core file can be very helpful. For help in running with a debugger you will need to contact your system administrators for the machine you are running on.

