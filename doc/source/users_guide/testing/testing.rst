.. _testing:

*******
Testing
*******

Technically, you could use the customization we gave in `Chapter 1 <CLM-URL>`_ to test various configuration and namelist options for CLM. 
Sometimes, it's also useful to have automated tests though to test that restarts give exactly the same results as without a restart. 
It's also useful to have automated tests to run over a wide variety of configurations, resolutions, and namelist options. 
To do that we have several different types of scripts set up to make running comprehensive testing of CLM easy. 
There are two types of testing scripts for CLM. 
The first are the CESM test scripts, which utilize the **create_newcase** scripts that we shown how to use in this User's Guide. 
The second are a set of stand-alone scripts that use the CLM **configure** and **build-namelist** scripts to build and test the model as well as testing the CLM tools as well. 
Below we will go into further details of how to use both methods.


CIME Testing scripts
====================

We first introduce the test scripts that work for all CESM components. 
The CIME script **create_test** runs a specific type of test, at a given resolution, for a given compset using a given machine. 
There is a list of different tests, but the "ERI" tests do several things at once, running from startup, as well as doing exact branch and restart tests. 
So to run "ERI" testing at 2-degree with the I1850CRUCLM45 compset on cheyenne_intel you do the following.
::

   > cd scripts
   > ./create_test -testname ERI.f19_g17_gl4.I1850CRUCLM45.cheyenne_intel
   > cd ERI.f19_g17_gl4.I1850CRUCLM45.cheyenne_intel.$id
   > ./ERI.f19_g17_gl4.I1850CRUCLM45.cheyenne_intel.$id.build
   > ERI.f19_g17_gl4.I1850CRUCLM45.cheyenne_intel.$id.submit

When the test is done it will update the file TestStatus with either a PASS or FAIL message.

We already have a standard list of tests for clm (the "aux_clm" list of tests). To run the CLM cheyenne intel compiler test list, for the same machine and compiler you would do the following:
::

   > cd scripts
   > ./create_test -xml_mach cheyenne -xml_compiler intel -xml_category aux_clm  -mach cheyenne -compiler intel
   # Normally it will submit the jobs as they are ready, but if it's interrupted you 
   # may need to submit by hand as follows...
   # Submit the suite of tests (note $id refers to the integer job number for this job)
   > ./cs.submit.$id.cheyenne
   # Later check the tests with...
   > ./cs.status.$id
   # The above will give a PASS or FAIL message for each test.

For more information on doing testing with the CESM scripts see the `|cesmrelease| User's Guide <CLM-URL>`_ on testing.

Testing PTCLM
=============

There is a simple test script for PTCLM called ``testcases.csh`` in the PTCLM directory (``scripts/ccsm_utils/Tools/lnd/clm/PTCLM``). 
The test script is setup to run on the machines: cheyenne, frankfurt, yong, and titan. 
You simply run the script interactively. 
The script will write out the status of tests to a file called: ``tc.job#.status``.

There are a few environment variables that can be used with ``testcases.csh`` to change it's operation.

``CESM_ROOT``: To test with a separate root to CESM code set this env variable to the root directory to use.
``CLM_SOFF``: If set to ``TRUE`` - stop on first failed test rather than continuing to run.
``CLM_RETAIN_FILES``: If set to ``FALSE`` - cleanup tools build first.
``DEBUG``: If set to ``TRUE`` - setup cases, but do not build or run.
