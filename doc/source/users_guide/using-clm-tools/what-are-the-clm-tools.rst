.. _what-are-the-clm-tools:

========================
 What are the CLM tools
========================

There are several tools provided with CLM that allow you to create your own input datasets at resolutions you choose, or to interpolate initial conditions to a different resolution, or used to compare CLM history files between different cases. 
The tools are all available in the ``$CTSMROOT/tools`` directory. 
Most of the tools are FORTRAN stand-alone programs in their own directory, but there is also a suite of NCL scripts in the ``./ncl_scripts`` directory, and some of the tools are scripts that may also call the ESMF regridding program. 
Some of the NCL scripts are very specialized and not meant for general use, and we won't document them here. 
They still contain documentation in the script itself and the README file in the tools directory.

The tools are divided into three directories for three categories: clm4_0, +|version|, and shared. 
The first two are of course for tools that are designed to work with either the CLM4.0 or +|version| versions of the model. 
The last one are shared utilities that can be used by either, or have a "-phys" option so you can specify which version you want to use.

The list of generally important scripts and programs are as follows.

1. *tools/cprnc* (relative to top level directory) to compare NetCDF files with a time axis.

#. *./mkmapgrids* to create SCRIP grid data files from old CLM format grid files that can then be used to create new CLM datasets (deprecated). There is also a NCL script (``./mkmapgrids/mkscripgrid.ncl`` to create SCRIP grid files for regular latitude/longitude grids.

#. *./mkmapdata* to create SCRIP mapping data file from SCRIP grid files (uses ESMF).

#. *./gen_domain* to create a domain file for datm from a mapping file. The domain file is then used by BOTH datm AND CLM to define the grid and land-mask.

#. *mksurfdata_map* to create surface datasets from grid datasets (clm4_0 and +|version| versions).

#. *./mkprocdata_map* to interpolate output unstructured grids (such as the CAM HOMME dy-core "ne" grids like ne30np4) into a 2D regular lat/long grid format that can be plotted easily. Can be used by either clm4_0 or +|version|.

In the sections to come we will go into detailed description of how to use each of these tools in turn. 
First, however we will discuss the common environment variables and options that are used by all of the FORTRAN tools. 
Second, we go over the outline of the entire file creation process for all input files needed by CLM for a new resolution, then we turn to each tool. 
In the last section we will discuss how to customize files for particular observational sites.

The tools run either one of two ways, with a namelist to provide options, or with command line arguments (and NOT both). 
**gen_domain** and **cprnc** run with command line arguments, and the other tools run with namelists.

In the following sections, we will outline how to make these files available for build-namelist so that you can easily create simulations that include them. 
In the chapter on single-point and regional datasets we also give an alternative way to enter new datasets without having to edit files.

------------------------------------
Running FORTRAN tools with namelists
------------------------------------

**mksurfdata_map** and **mkmapgrids** run with namelists that are read from standard input. 
Hence, you create a namelist and then run them by redirecting the namelist file into standard input as follows:
::

   ./program < namelist

For programs with namelists there is at least one sample namelist with the name "program".namelist (i.e. 
``mksurfdata_map.namelist`` for the **mksurfdata_map** program). 
There may also be other sample namelists that end in a different name besides "namelist". 
Namelists that you create should be similar to the example namelist. 
The namelist values are also documented along with the other namelists in the: 
::

   $CTSMROOT/bld/namelist_files/namelist_definition.xml`` file 
        and default values in the: 
   $CTSMROOT/bld/namelist_files/namelist_defaults_clm_tools.xml`` file.

-----------------------------------------------
Running FORTRAN tools with command line options
-----------------------------------------------

**gen_domain**, and **cprnc** run with command line arguments. 
The detailed sections below will give you more information on the command line arguments specific to each tool. 
Also running the tool without any arguments will give you a general synopsis on how to run the tool. 

-----------------------------------------
Running FORTRAN tools built with SMP=TRUE
-----------------------------------------

When you enable ``SMP=TRUE`` on your build of one of the tools that make use of it, you are using OpenMP for shared memory parallelism (SMP). 
In SMP loops are run in parallel with different threads run on different processors all of which access the same memory (called on-node). 
Thus you can only usefully run up to the number of processors that are available on a single-node of the machine you are running on. 
For example, on the NCAR machine cheyenne there are 16 processors per node, but the SMT hardware on the machine allows you to submit twice as many threads or 32 threads. 


---------
Using NCL
---------

In the tools directory ``$CTSMROOT/tools/ncl_scripts`` and in a few other locations there are scripts that use NCAR Command Language (NCL). 
Unlike the FORTRAN tools, you will need to get a copy of NCL in order to use them. 
You also won't have to build an executable in order to use them, hence no Makefile is provided. 
NCL is provided for free download as either binaries or source code from: `http://www.ncl.ucar.edu/ <http://www.ncl.ucar.edu/>`_. 
The NCL web-site also contains documentation on NCL and it's use. These scripts are stand-alone and at most use environment variables to control how to use them. In some cases there are perl scripts with command line arguments that call the NCL scripts to control what they do.
