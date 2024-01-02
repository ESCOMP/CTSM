.. highlight:: shell

.. _obtaining-and-building-ctsm:

=======================================
 Obtaining and building CTSM and LILAC
=======================================

This section describes the process for obtaining and building the CTSM library and its
dependencies, and linking to these libraries in an atmosphere model's build.

.. important::

   This documentation only applies to the process where you are building CTSM with LILAC
   for use in an atmosphere model that has *not* been integrated with CESM or CIME. If you
   are using CTSM within CESM, or running CTSM in land-only mode with a data atmosphere,
   then you should refer to the :ref:`general CTSM user's guide<users-guide>` as well as
   the `CIME documentation`_.

Quick start example / overview
==============================

The basic process for obtaining and building CTSM is the following:

Obtain CTSM by running::

  git clone https://github.com/ESCOMP/CTSM.git
  cd CTSM
  ./manage_externals/checkout_externals

Then build CTSM and its dependencies. On a machine that has been ported to CIME, the
command will look like this (example given for NCAR's ``cheyenne`` machine)::

  ./lilac/build_ctsm /glade/scratch/$USER/ctsm_build_dir --machine cheyenne --compiler intel

and then, before building the atmosphere model::

  source /glade/scratch/$USER/ctsm_build_dir/ctsm_build_environment.sh

On a machine that has *not* been ported to CIME, you will need to provide some additional
information. Run ``./lilac/build_ctsm -h`` for details, but the basic command will look
like this::

  ./lilac/build_ctsm ~/ctsm_build_dir --os Darwin --compiler gnu --netcdf-path /usr/local --esmf-mkfile-path /Users/sacks/ESMF/esmf8.0.0/lib/libO/Darwin.gfortranclang.64.mpich3.default/esmf.mk --max-mpitasks-per-node 4 --no-pnetcdf

In both cases, you will then need to include the necessary information in the include and
link lines of the atmosphere model's build. For a Makefile-based build, this can be done
by including the file ``/PATH/TO/CTSM/BUILD/ctsm.mk`` in the atmosphere model's build
scripts, then adding ``CTSM_INCLUDES`` to the include line and ``CTSM_LIBS`` to the link
line.

Further details on these steps are given below.

.. _building-ctsm-and-lilac-prerequisites:

Prerequisites
=============

Building CTSM requires:

- a Unix-like operating system (Linux, AIX, OS X, etc.)

- git version 1.8 or newer

- python3

  - Note that some scripts in the workflow look for 'python3' and others look for
    'python'. So python should be available under both of these names (although it is okay
    for ``python`` to refer to version 2.7.x).

- perl version 5

- a GNU version of the make tool (gmake)

- CMake

- Fortran and C compilers

  - See https://github.com/escomp/cesm#details-on-fortran-compiler-versions for
    information on compiler versions known to work with CESM, and thus CTSM.

- LAPACK and BLAS libraries

- a NetCDF library version 4.3 or newer built with the same compiler you will use for CTSM

  - a PNetCDF library is optional

- a functioning MPI environment

  - typically, this includes compiler wrappers like ``mpif90`` and ``mpicc``

- ESMF version 8 or later

  - **ESMF is not needed in general for CTSM, but is needed for LILAC**

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

CTSM provides a build script, ``lilac/build_ctsm``, for building CTSM and its dependencies. (The
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

.. _building-on-a-cime-supported-machine:

Building on a CIME-supported machine
------------------------------------

If you are using a machine that has been ported to CIME_ (for example, NCAR's ``cheyenne``
machine), then you do not need to specify much information to ``build_ctsm``. In addition,
in this case, CIME will load the appropriate modules and set the appropriate environment
variables at build time, so you do not need to do anything to set up your environment
ahead of time. **Building CTSM with LILAC requires ESMF. ESMF is currently an optional
CIME dependency, so many CIME-ported machines do not provide information on an ESMF
installation. NCAR's cheyenne machine DOES provide ESMF, but for other machines, you may
need to add this to your CIME port.**

To build CTSM and its dependencies in this case, run::

  ./lilac/build_ctsm /PATH/TO/CTSM/BUILD --machine MACHINE --compiler COMPILER

where you should fill in the capitalized arguments with appropriate values for your
machine.

.. note::

   The given directory (``/PATH/TO/CTSM/BUILD``) must *not* exist. This directory is
   created for you by the build script.

Some other options to ``build_ctsm`` are supported in this case (but many are not, since
they are only applicable to the non-CIME-supported machine workflow); run
``./lilac/build_ctsm -h`` for details.

.. important::

   If PNetCDF (parallel NetCDF) is not available on this machine, you will need to add the
   option ``--no-pnetcdf``.

   If you plan to run with OpenMP threading-based parallelization, or hybrid MPI/OpenMP,
   then it is important to add ``--build-with-openmp``.

Besides the build files themselves, ``build_ctsm`` creates the following important files
that are needed for the build of the atmosphere model:

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

Building on a machine that has not been ported to CIME
------------------------------------------------------

If you are using a machine thata has not been ported to CIME_, then you need to specify
additional information to ``build_ctsm`` that is needed by the build system. Before
building CTSM, you should load any modules and/or set any environment variables required
by the atmosphere model or CTSM builds, including all of the :ref:`prerequisites noted
above<building-ctsm-and-lilac-prerequisites>`.

The minimal amount of information needed is given by the following::

  ./lilac/build_ctsm /PATH/TO/CTSM/BUILD --os OS --compiler COMPILER --netcdf-path NETCDF_PATH --esmf-mkfile-path ESMF_MKFILE_PATH --max-mpitasks-per-node MAX_MPITASKS_PER_NODE --pnetcdf-path PNETCDF_PATH

where you should fill in the capitalized arguments with appropriate values for your
machine. Run ``./lilac/build_ctsm -h`` for details on these arguments, as well as documentation
of additional, optional arguments. Some of these optional arguments may be needed for
successful compilation, while others (such as ``--pnetcdf-path``) may be needed for good
model performance.

.. note::

   The given directory (``/PATH/TO/CTSM/BUILD``) must *not* exist. This directory is
   created for you by the build script.

.. important::

   If PNetCDF (parallel NetCDF) is not available on your machine/compiler, you should use
   the option ``--no-pnetcdf`` instead of ``--pnetcdf-path``. You must specify exactly one
   of those two options.

   If you plan to run with OpenMP threading-based parallelization, or hybrid MPI/OpenMP,
   then it is important to add ``--build-with-openmp``.

Example usage for a Mac (a simple case) is::

  ./lilac/build_ctsm ~/ctsm_build_dir --os Darwin --compiler gnu --netcdf-path /usr/local --esmf-mkfile-path /Users/sacks/ESMF/esmf8.0.0/lib/libO/Darwin.gfortranclang.64.mpich3.default/esmf.mk --max-mpitasks-per-node 4 --no-pnetcdf

Example usage for NCAR's ``cheyenne`` machine (a more complex case) is::

  module purge
  module load ncarenv/1.3 python/3.7.9 cmake intel/19.1.1 esmf_libs mkl
  module use /glade/p/cesmdata/cseg/PROGS/modulefiles/esmfpkgs/intel/19.1.1/
  module load esmf-8.2.0b23-ncdfio-mpt-O mpt/2.22 netcdf-mpi/4.8.0 pnetcdf/1.12.2 ncarcompilers/0.5.0

  ./lilac/build_ctsm /glade/scratch/$USER/ctsm_build_dir --os linux --compiler intel --netcdf-path '$ENV{NETCDF}' --pio-filesystem-hints gpfs --pnetcdf-path '$ENV{PNETCDF}' --esmf-mkfile-path '$ENV{ESMFMKFILE}' --max-mpitasks-per-node 36 --extra-cflags '-xCORE_AVX2 -no-fma' --extra-fflags '-xCORE_AVX2 -no-fma'

(It's better to use the :ref:`alternative process for a CIME-supported
machine<building-on-a-cime-supported-machine>` in this case, but the above illustrates
what would be needed for a machine similar to this that has not been ported to CIME.)

Besides the build files themselves, ``build_ctsm`` creates an important file that is
needed for the build of the atmosphere model: ``/PATH/TO/CTSM/BUILD/ctsm.mk``. This
Makefile-formatted file gives variables that should be set in the atmosphere model's
build. :ref:`See below for information on how to use this
file<including-ctsm-in-the-atmosphere-model-build>`.

Rebuilding after changing CTSM source code
------------------------------------------

To rebuild after changing CTSM source code, you should follow one of the above workflows,
but the ``build_ctsm`` command will simply be::

  ./lilac/build_ctsm /PATH/TO/CTSM/BUILD --rebuild

where ``/PATH/TO/CTSM/BUILD`` should point to the same directory you originally used.

.. _including-ctsm-in-the-atmosphere-model-build:

Including CTSM in the atmosphere model's build
==============================================

Once you have successfully built CTSM and its dependencies, you will need to add various
paths to the compilation and link lines when building your atmosphere model. For a
Makefile-based build system, we facilitate this by producing a file,
``/PATH/TO/CTSM/BUILD/ctsm.mk``, which you can include in your own build script. (We do
not yet produce an equivalent for CMake or other build systems.)

There are two important variables defined in this file:

- ``CTSM_INCLUDES``: This variable should be included in the compilation line for the
  atmosphere model's source files. It lists all paths that need to be included in these
  compilations so that the compiler can find the appropriate Fortran module files.

- ``CTSM_LIBS``: This variable should be included in the link line when creating the final
  executable. It lists paths and library names that need to be included in the link
  step. **Note: This may not include all of the libraries that are**
  :ref:`prerequisites<building-ctsm-and-lilac-prerequisites>`, **such as LAPACK, BLAS and
  NetCDF. If your atmosphere doesn't already require these, you may need to add
  appropriate information to your atmosphere model's link line.** However, it should
  already include all required link information for ESMF.

Other variables in this file do not need to be included directly in the atmosphere model's
build (they are just intermediate variables used to create ``CTSM_INCLUDES`` and
``CTSM_LIBS``).

For example, for the WRF build, we do the following: If building with CTSM, then we
expect that the user has set an environment variable::

  export WRF_CTSM_MKFILE=/PATH/TO/CTSM/BUILD/ctsm.mk

If that environment variable exists, then the ``configure`` script adds the following to
the Makefile-based build:

- Includes the ``ctsm.mk`` file (like ``include ${WRF_CTSM_MKFILE}``)

- Adds a CPP definition, ``-DWRF_USE_CTSM``, which is used to do conditional compilation
  of the CTSM-LILAC interface code

- Adds ``$(CTSM_INCLUDES)`` to its variable ``INCLUDE_MODULES``

- Adds ``$(CTSM_LIBS)`` to its variable ``LIB``

.. _CIME: http://esmci.github.io/cime
.. _CIME documentation: http://esmci.github.io/cime
