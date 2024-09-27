========================
 Building the CLM tools
========================

.. include:: ../substitutions.rst

.. todo::
  Update the below, as domain files aren't needed with nuopc.
The tools **cprnc** and **gen_domain** use the CIME configure/build system which is described in the next section.

The only CLM FORTRAN tool is mksurfdata_esmf which has it's own build system that takes advantage of the cime build.

================================================================
 Building the CLM tools that use the CIME configure/build system
================================================================

.. todo::
    Update the below, as domain files aren't needed with nuopc.

``cprnc`` and ``gen_domain`` both use the CIME configure/build system rather than the CLM specific version described above.

See `CIME documentation on adding grids <http://esmci.github.io/cime/users_guide/grids.html?highlight=gen_domain#adding-grids>`_ for more information on adding grids, creating mapping files, and running ``gen_domain``. Also see the CIME file: ``$CTSMROOT/tools/mapping/gen_domain_files/INSTALL`` for how to build ``gen_domain``.

