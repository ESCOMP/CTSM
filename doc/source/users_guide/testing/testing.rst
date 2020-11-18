.. _testing:

.. include:: ../substitutions.rst

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
See `CIME Chapter on Testing <http://esmci.github.io/cime/users_guide/testing.html>`_ for how to use it to run single
tests as well as lists of tests. The standard testname for CLM is "aux_clm" for cheyenne with intel and gnu compilers as
well as the CGD machine hobart for intel, nag, and pgi compilers.  There's also a shorter test list called "clm_short". Also
see the `CTSM Wiki on Testing <https://github.com/ESCOMP/ctsm/wiki/System-Testing-Guide>`_.

CTSM Tools Testing
==================

.. include:: ../../../../test/tools/README
   :literal:

CTSM Fortran Unit Tests
=======================

.. include:: ../../../../src/README.unit_testing
   :literal:

CTSM Build-namelist Tests
=========================

Run the following perl tester that 

::
   > cd bld/unit_testers
   > ./build-namelist_test.pl


Testing PTCLM
=============

.. include:: ../../../../tools/PTCLM/README
   :literal:

To run on cheyenne, you do the following:


.. include:: ../../../../tools/PTCLM/test/README.run_cheyenne
   :literal:
