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

Building CTSM and its dependencies, and including CTSM in the atmosphere model's build
======================================================================================

Overview
--------

CTSM provides a build script, ``build_ctsm`` for building CTSM and its dependencies. (The
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
it is generally simpler to use the second workflow below.

There is a third usage where you simply want to rebuild after making some source code
changes to CTSM. This is also documented below.

All of these workflows use CIME's build system behind the scenes. Typically, you will not
need to be aware of any of those details, but if problems arise, you may want to consult
the `CIME documentation <http://esmci.github.io/cime>`_.

Building on a CIME-supported machine
------------------------------------



Rebuilding after changing CTSM source code
------------------------------------------

To rebuild after changing CTSM source code, you should follow one of the above workflows,
but the ``build_ctsm`` command will simply be::

  ./build_ctsm /PATH/TO/CTSM/BUILD --rebuild

where ``/PATH/TO/CTSM/BUILD`` should point to the same directory you originally used.

.. _CIME: https://github.com/esmci/cime
