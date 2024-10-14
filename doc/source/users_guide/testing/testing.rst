.. include:: ../substitutions.rst

.. _testing:

*******
Testing
*******

Technically, you could use the customization we gave in :ref:`customizing_section` to test various configuration and namelist options for CLM. Sometimes, it's also useful to have automated tests though to test that restarts give exactly the same results as without a restart. It's also useful to have automated tests to run over a wide variety of configurations, resolutions, and namelist options. To do that we have several different types of scripts set up to make running comprehensive testing of CLM easy. There are two types of testing scripts for CLM. The first are the CESM test scripts, which utilize the ``cime/scripts/create_newcase`` scripts that we shown how to use in this User's Guide. The second are a set of stand-alone scripts that use the CLM ``configure`` and ``bld/build-namelist`` scripts to build and test the model as well as testing the CLM tools as well. Below we will go into further details of how to use both methods.

.. todo::
    Does ``configure`` script still exist?

CIME Testing scripts
====================

We first introduce the test scripts that work for all CESM components. The CIME script ``create_test`` runs a specific type of test, at a given resolution, for a given compset using a given machine. See `CIME Chapter on Testing <http://esmci.github.io/cime/users_guide/testing.html>`_ for how to use it to run single tests as well as lists of tests. The standard testname for CLM is "aux_clm" for cheyenne with intel and gnu compilers as well as the CGD machine hobart for intel, nag, and pgi compilers. There's also a shorter test list called "clm_short". Also see the `CTSM Wiki on Testing <https://github.com/ESCOMP/ctsm/wiki/System-Testing-Guide>`_.

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
