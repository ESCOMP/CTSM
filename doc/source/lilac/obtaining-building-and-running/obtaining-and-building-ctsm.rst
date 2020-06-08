.. _obtaining-and-building-ctsm:

.. highlight:: shell

=======================================
 Obtaining and building CTSM and LILAC
=======================================

This section describes the process for obtaining and building the CTSM library and its
dependencies, and linking to these libraries in an atmosphere model's build.

Quick start example
===================

The basic process for obtaining and building CTSM is the following:

Obtain CTSM by running::

  git clone https://github.com/ESCOMP/CTSM.git
  cd CTSM
  ./manage_externals/checkout_externals

Then build CTSM and its dependencies. On a machine that has been ported to CIME, the
command will look like this (example given for NCAR's ``cheyenne`` machine)::

  ./build_ctsm /glade/scratch/$USER/ctsm_build_dir --compiler intel --machine cheyenne

On a machine that has *not* been ported to CIME, you will need to provide some additional
information. Run ``build_ctsm -h`` for details, but the basic command will look like
this::

  ./build_ctsm ~/ctsm_build_dir --os Darwin --compiler gnu --netcdf-path /usr/local --esmf-lib-path /Users/sacks/ESMF/esmf8.0.0/lib/libO/Darwin.gfortranclang.64.mpich3.default

Further details on these commands are given below.

.. _building-ctsm-and-lilac-prerequisites:

Prerequisites
=============

Building CTSM requires:

- a Unix-like operating system (Linux, AIX, OS X, etc.)

- git version 1.8 or newer

- python3

  - The default version of python (when you run python without specifying 2 vs. 3) should
    be python3

- perl version 5

- a GNU version of the make tool (gmake)

- CMake

- Fortran and C compilers

  - See https://github.com/escomp/cesm#details-on-fortran-compiler-versions for
    information on compiler versions known to work with CESM, and thus CTSM.

- LAPACK and BLAS libraries

- a NetCDF library version 4.3 or newer built with the same compiler you will use for CTSM

  - a PnetCDF library is optional

- a functioning MPI environment

  - typically, this includes compiler wrappers like ``mpif90`` and ``mpicc``

Obtaining CTSM
==============

CTSM and its dependencies (excluding the :ref:`prerequisites noted
above<building-ctsm-and-lilac-prerequisites>`) can be obtained with::

  git clone https://github.com/ESCOMP/CTSM.git
  cd CTSM
  ./manage_externals/checkout_externals

By default, this will put you on the ``master`` branch of CTSM, which is the main
development branch. You can checkout a different branch or tag using ``git checkout``;
**be sure to rerun** ``./manage_externals/checkout_externals`` **after doing so.**

For more details, see
https://github.com/ESCOMP/CTSM/wiki/Quick-start-to-CTSM-development-with-git

Building CTSM and its dependencies
==================================

Overview
--------

CTSM provides a build script, ``build_ctsm``, for building CTSM and its dependencies. (The
dependencies built with this build script include various libraries that are packaged with
CIME_. This does *not* build the :ref:`prerequisites noted
above<building-ctsm-and-lilac-prerequisites>`: it is assumed that those are already built
on your machine.)

There are two possible workflows for building CTSM and its dependencies. The first works
if you are using a machine that has been ported to CIME_; the second works if you are
using a machine that has *not* been ported to CIME_. Both workflows are described
below. If you are using a machine that has not been ported to CIME, it is possible to do a
complete CIME port and then use the first workflow (by following the `CIME porting guide
<http://esmci.github.io/cime/versions/master/html/users_guide/porting-cime.html>`_), but
unless you need to do so for other reasons (such as running CESM, or running CTSM in a
land-only configuration forced by a data atmosphere, using the CIME_ scripting
infrastructure), it is generally simpler to use the second workflow below: A full CIME
port requires many settings that are not needed for just building CTSM.

There is a third usage where you simply want to rebuild after making some source code
changes to CTSM. This is also documented below.

All of these workflows use CIME's build system behind the scenes. Typically, you will not
need to be aware of any of those details, but if problems arise, you may want to consult
the `CIME documentation`_.

Building on a CIME-supported machine
------------------------------------

If you are using a machine that has been ported to CIME_ (for example, NCAR's ``cheyenne``
machine), then you do not need to specify much information to ``build_ctsm``. In addition,
in this case, CIME will load the appropriate modules and set the appropriate environment
variables at build time, so you do not need to do anything to set up your environment
ahead of time.

To build CTSM and its dependencies in this case, run::

  ./build_ctsm /PATH/TO/CTSM/BUILD --machine MACHINE --compiler COMPILER

where you should fill in the capitalized arguments with appropriate values for your
machine.

.. note::

   The given directory (``/PATH/TO/CTSM/BUILD``) must *not* exist. This directory is
   created for you by the build script.

Some other options to ``build_ctsm`` are supported in this case (but many are not, since
they are only applicable to the non-CIME-supported machine workflow); run ``build_ctsm
-h`` for details.

Besides the build files themselves, there are two key files that are needed for the build
of the atmosphere model:

1. ``/PATH/TO/CTSM/BUILD/ctsm.mk``: This Makefile-formatted file gives variables that
   should be set in the atmosphere model's build. :ref:`See below for information on how
   to use this file<including-ctsm-in-the-atmosphere-model-build>`.

2. ``/PATH/TO/CTSM/BUILD/ctsm_build_environment.sh`` or
   ``/PATH/TO/CTSM/BUILD/ctsm_build_environment.csh``: These files specify the build
   environment that CIME used to build CTSM and its dependencies. **Before building the
   atmosphere model, you should source the appropriate file** (based on your shell - use
   the ``.sh`` file for bash and similar shells, and the ``.csh`` file for tcsh and
   similar shells). **This will ensure that the atmosphere model is built with the same
   compiler and library versions as CTSM.** For example, with bash: ``source
   /PATH/TO/CTSM/BUILD/ctsm_build_environment.sh``.

Rebuilding after changing CTSM source code
------------------------------------------

To rebuild after changing CTSM source code, you should follow one of the above workflows,
but the ``build_ctsm`` command will simply be::

  ./build_ctsm /PATH/TO/CTSM/BUILD --rebuild

where ``/PATH/TO/CTSM/BUILD`` should point to the same directory you originally used.

.. _including-ctsm-in-the-atmosphere-model-build:

Including CTSM in the atmosphere model's build
==============================================

.. todo::

   TODO: Fill this section in

.. _CIME: http://esmci.github.io/cime
.. _CIME documentation: http://esmci.github.io/cime
