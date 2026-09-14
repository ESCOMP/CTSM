.. include:: ../substitutions.rst

.. _what-are-the-clm-tools:

========================
 What are the CLM tools
========================

There are tools provided with CLM that allow you, for example, to create your own input datasets at resolutions you choose or to compare CLM history files between different cases. Tools are available in the ``$CTSMROOT/tools`` directory. Some tools are FORTRAN stand-alone programs in their own directory. There is a suite of NCL scripts in the ``$CTSMROOT/tools/unsupported`` directory. Some of the tools also call the ESMF regridding program.

The tools produce files that can be used with |version|. If you need files for earlier versions of the model, you will likely need to use the tools present in the earlier versions.

The list of scripts and programs in ``$CTSMROOT/tools`` is as follows:

1. ``mksurfdata_esmf`` to create surface datasets from gridded datasets that we refer to as raw datasets (ctsm5_2 and newer versions).

#. ``crop_calendars`` to regrid and process GGCMI sowing and harvest date files for use in CTSM

#. ``site_and_regional`` contains scripts to handle input datasets for site and regional cases; these scripts help with creation of datasets using the standard process, or subsetting existing datasets, and also overwriting aspects for specific cases

#. ``modify_input_files`` contains scripts to modify CTSM input files, in particular surface datasets and mesh files

#. ``unsupported`` contains miscellaneous useful unsupported tools contributed by users; these tools may or may not work

#. ``cprnc`` to compare two NetCDF files. This is not in ``$CIMEROOT/tools`` and information is available in section :numref:`comparing-history-files`.

Subsequent sections provide details about these tools, while the following ``$CTSMROOT/tools/README.md`` goes through the complete process for creating input files needed to run the CLM:

.. include:: ../../../../tools/README.md
   :code: markdown

The ``$CTSMROOT/tools/unsupported/README.md`` covers what you need to know about the unsupported tools:

.. include:: ../../../../tools/unsupported/README.md
   :code: markdown

