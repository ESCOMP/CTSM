.. include:: ../substitutions.rst

.. _testing:

*******
Testing
*******

Technically, you could use the customization we gave in :ref:`customizing_section` to test various configuration and namelist options for CLM. Sometimes, it's also useful to have automated tests though to test that restarts give exactly the same results as without a restart. It's also useful to have automated tests to run over a wide variety of configurations, resolutions, and namelist options. To do that we have several different types of scripts set up to make running comprehensive testing of CLM easy. There are four types of testing scripts for CLM. The first is the ``./run_sys_tests`` script at the top level that runs the CESM test scripts, which utilize the CIME test infrastructure. The second are the python unit and system tests under the python directory. The third are the PF unit-tests for the Fortran code, that are also run as part of the system tests in "1". The fourth is a stand-alone script that use the CLM ``bld/build-namelist`` scripts to build and test namelist creation for the model. Below we will go into further details of how to use each of these methods.

#. run_sys_tests
#. python tests
#. PFunit Fortran unit tests
#. namelist testing

CIME Testing scripts
====================

We first introduce the test scripts that work for all CESM components. The CIME script ``cime/scripts/create_test`` runs a specific type of test, at a given resolution, for a given compset using a given machine. See `CIME Chapter on Testing <http://esmci.github.io/cime/users_guide/testing.html>`_ for how to use it to run single tests as well as lists of tests. The standard testname for CLM is "aux_clm" for derecho with intel and gnu compilers as well as the CGD machine izumi for intel, nag, and pgi compilers. There's also a shorter test list called "clm_short". Also see the `CTSM Wiki on Testing <https://github.com/ESCOMP/ctsm/wiki/System-Testing-Guide>`_.

run_sys_tests for CTSM test suites
==================================
run_sys_tests is a script at the top level of the CTSM repository that runs the CIME test scripts for a test suite specified in the testlist XML file. To use it, simply do:

::

   > /run_sys_tests -s aux_clm --skip_compare --skip-generate

Which will run the ``aux_clm`` test suite, skipping the compare and generate steps. Other common CTSM test suites include: fates, clm_short, and ctsm_sci.

CTSM Python Tests
=================
The CTSM python tests are located in the ``python/ctsm/test`` directory. These include unit tests for the python code as well as system tests that run the python code to test various configurations and namelist options. To run all of the python tests, simply do "make" in the python directory.

::

   > cd python
   > make

CTSM Fortran Unit Tests
=======================

.. include:: ../../../../src/README.unit_testing
   :literal:

CTSM Build-namelist Tests
=========================

Test the namelist build script by running the following:

::

   > cd bld/unit_testers
   > ./build-namelist_test.pl 1>namelist_test.log 2>&1

When that's complete, inspect ``namelist_test.log`` (e.g., with ``less namelist_test.log``). If you see ``Successfully ran all testing for build-namelist`` but nothing like ``# Looks like you failed 4 tests of 1999.``, then everything went fine.

If something went wrong, you can find the failing tests like so:

::

   > grep -E "^[0-9]+/[0-9]+ < [a-zA-Z]+" namelist_test.log | grep -v "PASS"
