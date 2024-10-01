.. include:: ../substitutions.rst

.. _what-are-the-clm-tools:

========================
 What are the CLM tools
========================

.. todo::
    Remove references to mkprocdata_map?

There are several tools provided with CLM that allow you to create your own input datasets at resolutions you choose, or to interpolate initial conditions to a different resolution, or used to compare CLM history files between different cases. The tools are all available in the ``$CTSMROOT/tools`` directory. Most of the tools are FORTRAN stand-alone programs in their own directory, but there is also a suite of NCL scripts in the ``$CTSMROOT/tools//ncl_scripts`` directory, and some of the tools are scripts that may also call the ESMF regridding program. Some of the NCL scripts are very specialized and not meant for general use, and we won't document them here. They still contain documentation in the script itself and the README file in the tools directory.

The tools produce files that can be used for CLM4.5 and |version|. They do **NOT** produce files that can be used for CLM4.0. If you need files for CLM4.0, you'll need to use a previous version of CLM.

The list of generally important scripts and programs are as follows.

1. ``./mkmapgrids`` to create SCRIP grid data files from old CLM format grid files that can then be used to create new CLM datasets (deprecated). There is also a NCL script (``./mkmapgrids/mkscripgrid.ncl``) to create SCRIP grid files for regular latitude/longitude grids.

#. ``./mkmapdata`` to create SCRIP mapping data file from SCRIP grid files (uses ESMF).

#. ``mksurfdata_esmf`` to create surface datasets from grid datasets (clm4_0 and |version| versions).

.. todo::
    Update the below, as domain files aren't needed with nuopc.
#. *$CIMEROOT/tools/mapping/gen_domain_files/gen_domain* to create a domain file for datm from a mapping file. The domain file is then used by BOTH datm AND CLM to define the grid and land-mask.

#. ``$CIMEROOT/tools/mapping/gen_domain_files/gen_domain`` to create a domain file for datm from a mapping file. The domain file is then used by BOTH datm AND CLM to define the grid and land-mask.

#. ``$CIMEROOT/tools/cprnc`` to compare two NetCDF files.

In the sections to come we will go into detailed description of how to use each of these tools in turn. First, however we will discuss the common environment variables and options that are used by all of the FORTRAN tools. Second, we go over the outline of the entire file creation process for all input files needed by CLM for a new resolution, then we turn to each tool. In the last section we will discuss how to customize files for particular observational sites.

The FORTRAN tool (``mksurfdata_esmf``) runs, with a namelist and has a namelist builder for it.

In the following sections, we will outline how to make these files available for build-namelist so that you can easily create simulations that include them. In the chapter on single-point and regional datasets we also give an alternative way to enter new datasets without having to edit files.

------------------------------------
Running FORTRAN tools with namelists
------------------------------------

``mksurfdata_esmf`` runs with a namelist that is read from standard input. Hence, you create a namelist and then run them by redirecting the namelist file into standard input as follows:
::

   ./program < namelist

**mksurfdata_esmf** also has a script that creates the namelist and runs the program for you. The namelist values are also documented along with the other namelists in the:
::

   $CTSMROOT/tools/mksurfdata_esmf/gen_mksurfdata_namelist.xml`` file

-----------------------------------------------
Running FORTRAN tools with command line options
-----------------------------------------------

.. todo::
    Update the below, as domain files aren't needed with nuopc.

**gen_domain** and **cprnc** run with command line arguments. The detailed sections below will give you more information on the command line arguments specific to each tool. Also running the tool without any arguments will give you a general synopsis on how to run the tool.
