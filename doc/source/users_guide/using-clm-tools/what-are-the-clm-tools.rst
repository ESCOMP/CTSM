.. include:: ../substitutions.rst

.. _what-are-the-clm-tools:

========================
 What are the CLM tools
========================

There are several tools provided with CLM that allow you to create your own input datasets at resolutions you choose, or to interpolate initial conditions to a different resolution, or used to compare CLM history files between different cases. The tools are all available in the ``$CTSMROOT/tools`` directory. Most of the tools are FORTRAN stand-alone programs in their own directory, but there is also a suite of NCL scripts in the ``$CTSMROOT/tools//ncl_scripts`` directory, and some of the tools are scripts that may also call the ESMF regridding program. Some of the NCL scripts are very specialized and not meant for general use, and we won't document them here. They still contain documentation in the script itself and the README file in the tools directory.

The tools produce files that can be used for CLM4.5 and |version|. They do **NOT** produce files that can be used for CLM4.0. If you need files for CLM4.0, you'll need to use a previous version of CLM.

The list of generally important scripts and programs are as follows.

1. *./mkmapgrids* to create SCRIP grid data files from old CLM format grid files that can then be used to create new CLM datasets (deprecated). There is also a NCL script (``./mkmapgrids/mkscripgrid.ncl`` to create SCRIP grid files for regular latitude/longitude grids.

#. *./mkmapdata* to create SCRIP mapping data file from SCRIP grid files (uses ESMF).

#. *mksurfdata_esmf* to create surface datasets from grid datasets (clm4_0 and |version| versions).

#. *./mkprocdata_map* to interpolate output unstructured grids (such as the CAM HOMME dy-core "ne" grids like ne30np4) into a 2D regular lat/long grid format that can be plotted easily. Can be used by either clm4_0 or |version|.

#. *$CIMEROOT/tools/mapping/gen_domain_files/gen_domain* to create a domain file for datm from a mapping file. The domain file is then used by BOTH datm AND CLM to define the grid and land-mask.

#. *$CIMEROOT/tools/cprnc* to compare two NetCDF files.

In the sections to come we will go into detailed description of how to use each of these tools in turn. First, however we will discuss the common environment variables and options that are used by all of the FORTRAN tools. Second, we go over the outline of the entire file creation process for all input files needed by CLM for a new resolution, then we turn to each tool. In the last section we will discuss how to customize files for particular observational sites.

The FORTRAN tools (mksurfdata_esmf and mkprocdata_map) run, with a namelist (mksurfdata_esmf) to provide options, or with command line arguments (mkprocdata_map).

In the following sections, we will outline how to make these files available for build-namelist so that you can easily create simulations that include them. In the chapter on single-point and regional datasets we also give an alternative way to enter new datasets without having to edit files.

------------------------------------
Running FORTRAN tools with namelists
------------------------------------

**mksurfdata_esmf** runs with a namelist that is read from standard input. Hence, you create a namelist and then run them by redirecting the namelist file into standard input as follows:
::

   ./program < namelist

There is a sample namelist called ``$CTSMROOT/tools/mksurfdata_esmf/mksurfdata_esmf.namleist`` that shows you what the namelist should look like. **mksurfdata_esmf** also has a script that creates the namelist and runs the program for you. Namelists that you create should be similar to the example namelist. The namelist values are also documented along with the other namelists in the:
::

   $CTSMROOT/bld/namelist_files/namelist_definition.xml`` file
        and default values in the:
   $CTSMROOT/bld/namelist_files/namelist_defaults_clm_tools.xml`` file.

-----------------------------------------------
Running FORTRAN tools with command line options
-----------------------------------------------

**gen_domain**, mkprocdata_map, and **cprnc** run with command line arguments. The detailed sections below will give you more information on the command line arguments specific to each tool. Also running the tool without any arguments will give you a general synopsis on how to run the tool.

-----------------------------------------
Running FORTRAN tools built with SMP=TRUE
-----------------------------------------

When you enable ``SMP=TRUE`` on your build of one of the tools that make use of it, you are using OpenMP for shared memory parallelism (SMP). In SMP loops are run in parallel with different threads run on different processors all of which access the same memory (called on-node). Thus you can only usefully run up to the number of processors that are available on a single-node of the machine you are running on. For example, on the NCAR machine cheyenne there are 36 processors per node, so you can use up to 36 processors.

.. _using-ncl:

---------
Using NCL
---------

In the tools directory ``$CTSMROOT/tools/ncl_scripts`` and in a few other locations there are scripts that use NCAR Command Language (NCL). Unlike the FORTRAN tools, you will need to get a copy of NCL in order to use them. You also won't have to build an executable in order to use them, hence no Makefile is provided. NCL is provided for free download as either binaries or source code from: `http://www.ncl.ucar.edu/ <http://www.ncl.ucar.edu/>`_. The NCL web-site also contains documentation on NCL and it's use. These scripts are stand-alone and at most use environment variables to control how to use them. In some cases there are perl scripts with command line arguments that call the NCL scripts to control what they do.
