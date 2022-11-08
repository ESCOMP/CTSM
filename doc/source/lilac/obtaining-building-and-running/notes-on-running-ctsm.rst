.. highlight:: shell

.. _notes-on-running-ctsm:

=======================
 Notes on running CTSM
=======================

.. _runtime-environment-variables:

Environment variables that may need to be set at runtime
========================================================

With the MPT MPI library (which is the default MPI library on NCAR's cheyenne machine), it is important to set the environment variable ``MPI_TYPE_DEPTH`` to 16 when running CTSM (this setting is required by the Parallel IO library). Typically you should set this variable in your job submission script, using either:

.. code-block:: Bash

    export MPI_TYPE_DEPTH=16

or:

.. code-block:: Tcsh

    setenv MPI_TYPE_DEPTH 16

prior to running the model.
