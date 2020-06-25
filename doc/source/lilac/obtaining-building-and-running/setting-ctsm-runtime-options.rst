.. highlight:: shell

.. _setting-ctsm-runtime-options:

==============================
 Setting CTSM runtime options
==============================

Overview and quick start
========================

This section describes the process for creating the runtime input text files for CTSM and
LILAC. These files, which are in Fortran namelist format, have hard-coded file
names. These files must exist with the expected names in the directory from which the
model is run:

- ``lnd_in``: This is the main namelist input file for CTSM

- ``lnd_modelio.nml``: This sets CTSM's PIO (parallel i/o library) configuration settings

- ``lilac_in``: This namelist controls the operation of LILAC

.. note::

   There are a number of other required runtime input files to both CTSM and LILAC, in
   NetCDF format. The paths to these other files are specified in either ``lnd_in`` or
   ``lilac_in``.

The basic process for creating the necessary input files is the following; this process is
also illustrated in :numref:`Figure ctsm_lilac_runtime_file_workflow`:

#. Run the ``build_ctsm`` script described in section
   :numref:`obtaining-and-building-ctsm`. In addition to building CTSM, this also stages
   the necessary files in the ``runtime_inputs`` subdirectory of your specified build
   directory. Then ``cd`` to this ``runtime_inputs`` subdirectory to do the following
   steps (it is fine to do these steps even while CTSM is still building).

#. Modify the ``ctsm.cfg`` file to set high-level options to CTSM. (A few options need to
   be set; most can be left at their default values or changed if desired.) Optionally,
   also set specific namelist values in ``user_nl_ctsm``.

#. Run the script, ``make_runtime_inputs``. (This creates the files ``lnd_in`` and
   ``clm.input_data_list``.)

#. Modify ``lilac_in`` as needed. (Typically you will only need to set values for
   ``atm_mesh_filename`` and ``lnd_mesh_filename``; other variables can typically be kept
   at their default values.)

#. Run the script, ``download_input_data`` to download any of CTSM's standard input files
   that are needed based on settings in ``lnd_in`` and ``lilac_in``. (This step may be
   unnecessary if all of the needed input data already exists. However, it doesn't hurt to
   run it in this case.)

#. Copy ``lnd_in``, ``lnd_modelio.nml`` and ``lilac_in`` to the directory from which you
   will be running the model.

.. _Figure ctsm_lilac_runtime_file_workflow:

.. figure:: ctsm_lilac_runtime_file_workflow.*

   CTSM/LILAC runtime file workflow. Files in black can be (and in some cases must be)
   edited before running the next step. Files in blue should **not** be edited directly.

