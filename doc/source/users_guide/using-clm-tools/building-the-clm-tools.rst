========================
 Building the CLM tools
========================

.. include:: ../substitutions.rst

The CLM FORTRAN tools all have similar makefiles, and similar options for building. The tools **cprnc** and **gen_domain** use the CIME configure/build system which is described in the next section.

The Makefiles (for **mksurfdata_esmf** and **mkprocdata_map**) use GNU Make extensions and thus require that you use GNU make to use them. They also auto detect the type of platform you are on, using "uname -s" and set the compiler, compiler flags and such accordingly. There are also environment variables that can be set to set things that must be customized. All the tools use NetCDF and hence require the path to the NetCDF libraries and include files. On some platforms (such as Linux) multiple compilers can be used, and hence there are env variables that can be set to change the FORTRAN and/or "C" compilers used. The tools also allow finer control, by also allowing the user to add compiler flags they choose, for both FORTRAN and "C", as well as picking the compiler, linker and and add linker options. Finally the tools allow you to turn optimization on (which is off by default but on for **mksurfdata_esmf**) with the OPT flag so that the tool will run faster.

Options used by all:  **mksurfdata_esmf**

- ``LIB_NETCDF`` -- sets the location of the NetCDF library.
- ``INC_NETCDF`` -- sets the location of the NetCDF include files.
- ``USER_FC`` -- sets the name of the FORTRAN compiler.

Options used by: **mkprocdata_map**, and **mksurfdata_esmf**

- ``MOD_NETCDF`` -- sets the location of the NetCDF FORTRAN module.
- ``USER_LINKER`` -- sets the name of the linker to use.
- ``USER_CPPDEFS`` -- adds any CPP defines to use.
- ``USER_CFLAGS`` -- add any "C" compiler flags to use.
- ``USER_FFLAGS`` -- add any FORTRAN compiler flags to use.
- ``USER_LDFLAGS`` -- add any linker flags to use.
- ``USER_CC`` -- sets the name of the "C" compiler to use.
- ``OPT`` -- set to TRUE to compile the code optimized (TRUE or FALSE)
- ``SMP`` -- set to TRUE to turn on shared memory parallelism (i.e. OpenMP) (TRUE or FALSE)
- ``Filepath`` -- list of directories to build source code from.
- ``Srcfiles`` -- list of source code filenames to build executable from.
- ``Makefile`` -- customized makefile options for this particular tool.
- ``mkDepends`` -- figure out dependencies between source files, so make can compile in order..
- ``Makefile.common`` -- General tool Makefile that should be the same between all tools.

More details on each environment variable.

``LIB_NETCDF``
  This variable sets the path to the NetCDF library file (``libnetcdf.a``). If not set it defaults to ``/usr/local/lib``. In order to use the tools you need to build the NetCDF library and be able to link to it. In order to build the model with a particular compiler you may have to compile the NetCDF library with the same compiler (or at least a compatible one).

``INC_NETCDF``
  This variable sets the path to the NetCDF include directory (in order to find the include file ``netcdf.inc``). if not set it defaults to ``/usr/local/include``.

``MOD_NETCDF``
  This variable sets the path to the NetCDF module directory (in order to find the NetCDF FORTRAN-90 module file when NetCDF is used with a FORTRAN-90 **use statement**. When not set it defaults to the ``LIB_NETCDF`` value.

``USER_FC``
  This variable sets the command name to the FORTRAN-90 compiler to use when compiling the tool. The default compiler to use depends on the platform. And for example, on the AIX platform this variable is NOT used

``USER_LINKER``
  This variable sets the command name to the linker to use when linking the object files from the compiler together to build the executable. By default this is set to the value of the FORTRAN-90 compiler used to compile the source code.

``USER_CPPDEFS``
  This variable adds additional optional values to define for the C preprocessor. Normally, there is no reason to do this as there are very few CPP tokens in the CLM tools. However, if you modify the tools there may be a reason to define new CPP tokens.

``USER_CC``
  This variable sets the command name to the "C" compiler to use when compiling the tool. The default compiler to use depends on the platform. And for example, on the AIX platform this variable is NOT used

``USER_CFLAGS``
  This variable adds additional compiler options for the "C" compiler to use when compiling the tool. By default the compiler options are picked according to the platform and compiler that will be used.

``USER_FFLAGS``
  This variable adds additional compiler options for the FORTRAN-90 compiler to use when compiling the tool. By default the compiler options are picked according to the platform and compiler that will be used.

``USER_LDFLAGS``
  This variable adds additional options to the linker that will be used when linking the object files into the executable. By default the linker options are picked according to the platform and compiler that is used.

``SMP``
  This variable flags if shared memory parallelism (using OpenMP) should be used when compiling the tool. It can be set to either TRUE or FALSE, by default it is set to FALSE, so shared memory parallelism is NOT used. When set to TRUE you can set the number of threads by using the OMP_NUM_THREADS environment variable. Normally, the most you would set this to would be to the number of on-node CPU processors. Turning this on should make the tool run much faster.

.. warning:: Note, that depending on the compiler answers may be different when SMP is activated.

``OPT``
  This variable flags if compiler optimization should be used when compiling the tool. It can be set to either ``TRUE`` or ``FALSE``, by default it is set to for both **mksurfdata_esmf** and **mkprocdata_map**. Turning this on should make the tool run much faster.

.. warning:: Note, you should expect that answers will be different when ``OPT`` is activated.

``Filepath``
  All of the tools are stand-alone and don't need any outside code to operate. The Filepath is the list of directories needed to compile and hence is always simply "." the current directory. Several tools use copies of code outside their directory that is in the CESM distribution (either ``csm_share`` code or CLM source code).

``Srcfiles``
  The ``Srcfiles`` lists the filenames of the source code to use when building the tool.

``Makefile``
  The ``Makefile`` is the custom GNU Makefile for this particular tool. It will customize the ``EXENAME`` and the optimization settings for this particular tool.

``Makefile.common``
  The ``Makefile.common`` is the copy of the general GNU Makefile for all the CLM tools. This file should be identical between the different tools. This file has different sections of compiler options for different Operating Systems and compilers.

``mkDepends``
  The ``mkDepends`` is the copy of the perl script used by the ``Makefile.common`` to figure out the dependencies between the source files so that it can compile in the necessary order. This file should be identical between the different tools.

.. note:: There are several files that are copies of the original files. By having copies the tools can all be made stand-alone, but any changes to the originals will have to be put into the tool directories as well.

The *README.filecopies* (which can be found in ``$CTSMROOT/tools``) is repeated here.

.. include:: ../../../../tools/README.filecopies
   :literal:

================================================================
 Building the CLM tools that use the CIME configure/build system
================================================================

**cprnc** and *gen_domain** both use the CIME configure/build system rather than the CLM specific version described above.

See `CIME documentation on adding grids <http://esmci.github.io/cime/users_guide/grids.html?highlight=gen_domain#adding-grids>`_ for more information on adding grids, creating mapping files, and running **gen_domain**. Also see the CIME file: ``$CTSMROOT/tools/mapping/gen_domain_files/INSTALL`` for how to build **gen_domain**.

