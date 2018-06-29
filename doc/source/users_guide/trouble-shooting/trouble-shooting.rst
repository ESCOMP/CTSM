.. _trouble-shooting:

.. include:: ../substitutions.rst

*********************
Trouble Shooting 
*********************

In this chapter we give some guidance on what to do when you encounter some of the most common problems. We can't cover all the problems that a user could potentially have, but we will try to help you recognize some of the most common situations. And we'll give you some suggestions on how to approach the problem to come up with a solution.

In general you will run into one of three type of problems:

1. *case-creation*
#. *setup-time*
#. *build-time*
#. *run-time*

See the `CIME Trouble Shooting Guide <http://esmci.github.io/cime/users_guide/troubleshooting.html>`_ for some help on the first three.


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

Another suggestion is to run the model with a debugger such as: **ddt**, **dbx**, **gdb**, or **totalview**. 
Often to run with a debugger you will need to reduce the number of processors as outlined above. 
Some debuggers such as **dbx** will only work with one processor, while more advanced debuggers such as **totalview** can work with both MPI tasks and OMP threads. 
Even simple debuggers though can be used to query core files, to see where the code was at when it died (for example using the **where** in **dbx** for a core file can be very helpful. 
For help in running with a debugger you will need to contact your system administrators for the machine you are running on.

