.. highlight:: shell

.. _notes-on-running-ctsm:

=======================
 Notes on running CTSM
=======================

.. _runtime-environment-variables:

Environment variables that may need to be set at runtime
========================================================

Currently none. This is only a placeholder. Typically you should set this variable in your job submission script, using either:

.. code-block:: Bash

    export <NONE>=<value>

or:

.. code-block:: Tcsh

    setenv <NONE> <value>

prior to running the model.
