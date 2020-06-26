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

   CTSM/LILAC runtime file workflow. Files in blue can be (and in some cases must be)
   edited before running the next step. Files in purple (with italicized names) should
   **not** be edited directly.

More details on these steps are given in the following subsections.

Creating initial runtime inputs with build_ctsm
===============================================

The ``build_ctsm`` script, which is described in detail in section
:numref:`obtaining-and-building-ctsm`, creates initial runtime input files in addition to
building the model. This script creates a number of files in the ``runtime_inputs``
subdirectory of the specified build directory. For a few variables in these runtime input
files, ``build_ctsm`` sets initial values based on options provided to this
script. Important options for these runtime inputs include ``--no-pnetcdf``,
``--inputdata-path`` and ``--max-mpitasks-per-node``. (Run ``build_ctsm`` with the ``-h``
or ``--help`` option for more information.)

Once this script creates and populates the ``runtime_inputs`` subdirectory, it is safe to
proceed with the following steps, even if CTSM has not finished building.

For the following steps, you should ``cd`` to this ``runtime_inputs`` subdirectory.

Modifying ctsm.cfg and user_nl_ctsm
===================================

CTSM has hundreds of runtime parameters. Most of these parameters can be set individually,
but in many cases it makes more sense to think in terms of high-level options. These
high-level options set groups of parameters, creating configurations that the core CTSM
developers feel are useful - and these standard configurations are generally tested both
from a scientific and software perspective.

The two text files, ``ctsm.cfg`` and ``user_nl_ctsm``, together with CTSM's scripting
infrastructure and XML database controlled by the ``make_runtime_inputs`` script, work
together to allow you to configure CTSM's runtime parameters at both a high level and
individually.

ctsm.cfg
--------

``ctsm.cfg`` controls high-level options that, in many cases, set the default values for
multiple individual runtime parameters. All of the available high-level options appear in
this file; you can change the values of variables, but cannot add or remove any variables
in this file.

The first set of options in this file specifies key file names:

- ``lnd_domain_file`` must be specified. This file specifies CTSM's grid and land
  mask. The general process for creating this file is described in section
  :numref:`creating-domain-files`.

- ``fsurdat`` also must be specified. This file specifies a variety of spatially-varying
  properties. This file is grid-specific, but can be created from grid-independent files
  using CTSM's toolchain described in section :numref:`creating-surface-datasets`.

- ``finidat`` should generally be specified, although it's not absolutely essential. This
  file specifies CTSM's initial conditions. If this isn't specified, the model will use a
  standard set of initial conditions, interpolated to your grid. However, particularly for
  NWP / prediction applications, you will typically want a customized initial condition
  file. The process for generating this file will depend on your atmosphere model and
  workflow, but an example for WRF is given in section
  :numref:`wrf-create-input-namelists-for-ctsm-and-lilac`.

The remainder of this file specifies a variety of high-level options, each of which sets
the default values for a number of CTSM's runtime parameters. The default values should be
reasonable starting points, but you may want to configure these. Details on these options
and allowed values are given in comments in ``ctsm.cfg``.

user_nl_ctsm
------------

This file allows you to override individual CTSM namelist variables. This includes
variables whose default values are set based on settings in ``ctsm.cfg`` and others. The
file is initially populated with some settings controlling CTSM's diagnostic (history)
file output. These pre-populated settings can be changed, and additional settings can be
added to this file.

There is some documentation of these settings in section :numref:`customizing-a-case`, and
in the `CESM release documentation
<http://www.cesm.ucar.edu/models/cesm2/settings/current/clm5_0_nml.html>`_, but note that
the latter is slightly out of date with respect to the latest version of CTSM. An easy way
to see the list of available variables is to run ``make_runtime_inputs`` in order to
generate an initial ``lnd_in`` file; most of the variables given in that file can be
specified in ``user_nl_ctsm``, and then ``make_runtime_inputs`` can be rerun. **As noted
below, it is better NOT to edit the** ``lnd_in`` **file directly, instead using the
workflow documented here.**

Running make_runtime_inputs
===========================

Once you have made the modifications you want to ``ctsm.cfg`` and ``user_nl_ctsm``, run
the script ``make_runtime_inputs`` from the ``runtime_inputs`` directory. This takes
``ctsm.cfg`` and ``user_nl_ctsm`` as inputs, and generates two output files: ``lnd_in``
and ``clm.input_data_list``. ``lnd_in`` will be read by CTSM. ``clm.input_data_list`` is
an automatic extraction of a subset of ``lnd_in`` specifying the paths of various other
input files that will be needed by CTSM; this is used by the ``download_input_data``
script to automatically download the relevant files.

It is safe to rerun ``make_runtime_inputs`` as often as you want, incrementally changing
``ctsm.cfg`` and/or ``user_nl_ctsm``.

.. important::

   We recommend that you do NOT modify ``lnd_in`` directly. Instead, to make changes to
   the ``lnd_in`` file, you should modify ``user_nl_ctsm`` and rerun
   ``make_runtime_inputs``. There are a few reasons for following this workflow:

   - Hand edits to ``lnd_in`` will be lost if you later rerun ``make_runtime_inputs``,
     whereas edits to ``user_nl_ctsm`` will be maintained.

   - ``make_runtime_inputs`` performs various validations of the contents of
     ``user_nl_ctsm``; these validations would be bypassed if you edited ``lnd_in``
     directly.

   - If you change any file paths, ``make_runtime_inputs`` will ensure that
     ``clm.input_data_list`` remains in sync with ``lnd_in``.

Modifying lilac_in
==================

Unlike ``lnd_in``, the ``lilac_in`` file can be hand-edited. Most of the settings in this
file can be left at their default values, but there are two variables whose values you
must set (as indicated by their default values, ``FILL_THIS_IN``):

- ``atm_mesh_filename``: This should specify the path to an ESMF mesh file describing the
  atmosphere model's grid.

- ``lnd_mesh_filename``: This should specify the path to an ESMF mesh file describing the
  land model's grid. If the land model is running on the same grid as the atmosphere
  model (which is typical), this can be the same file as ``atm_mesh_filename``.

Other settings you may want to change are:

- Settings in ``lilac_history_input``: ``lilac_histfreq_option`` and
  ``lilac_histfreq_n``. Together, these specify the output frequency from LILAC
  itself. Note that this is separate from CTSM's output: LILAC's output contains
  instantaneous snapshots of the fields passed from the atmosphere to CTSM and vice
  versa, whereas CTSM's output is much more extensive. For many purposes, it's fine to
  leave LILAC's output turned off (as is the default). Allowable options for
  ``lilac_histfreq_option`` are ``never``, ``nsteps``, ``nseconds``, ``nminutes``,
  ``nhours``, ``ndays``, ``nmonths`` and ``nyears``.

- Settings in ``atmaero_stream``: These specify a dataset containing atmospheric aerosols,
  for the (typical) case where the atmosphere model is not sending these aerosols itself.

Running download_input_data
===========================

CTSM requires a variety of runtime input files in NetCDF format. These files are listed in
the ``lnd_in`` file, and are consolidated in the file ``clm.input_data_list`` (which is
produced by ``make_runtime_inputs``). In addition, a few other NetCDF files are listed in
``lilac_in``, of which the file listed in ``atmaero_stream`` is typically a standard input
file (as opposed to one that you, the user, has provided).

**Idea of inputdata directory... this can be a symlink, or it can be an actual directory,
depending on how build_ctsm was run.**

**As a convenience, can run download_input_data.**

**There will likely be errors about some files not being able to be obtained, particularly
for user-provided files.**

Copying the necessary files to the model's run directory
========================================================

.. todo::

   TODO: Fill this section in
