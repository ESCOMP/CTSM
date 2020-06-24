.. _setting-ctsm-runtime-options:

.. highlight:: shell

==============================
 Setting CTSM runtime options
==============================

Overview
========

This section describes the process for creating the runtime input text files for CTSM and
LILAC. These files, which are in Fortran namelist format, have hard-coded file
names. These files must exist with the expected names in the directory from which the
model is run:

- ``lnd_in``: This is the main namelist input file for CTSM

- ``lnd_modelio.nml``: This sets CTSM's PIO (parallel i/o library) configuration settings

- ``lilac_in``: This namelist controls the operation of LILAC

.. note::

   There are a number of other runtime input files to both CTSM and LILAC, in NetCDF
   format. The paths to these other files are specified in either ``lnd_in`` or
   ``lilac_in``.

