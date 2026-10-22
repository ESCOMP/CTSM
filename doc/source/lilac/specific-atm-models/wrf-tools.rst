.. _wrf:

.. highlight:: shell

=============================
 WRF-CTSM Tools and Utilities
=============================

This section includes instructions on tools and utilities developed for
WRF-CTSM simulations.

Generate CTSM surface dataset for a WRF domain
----------------------------------------------

Before this step, make sure you have successfully created geo_em* files for
your specific WRF domain using WPS. Instructions on how to run ``geogrid.exe``
is described in here.

1. Create ESMF mesh file from WRF ``geo_em*`` files, using the make_mesh tool. Details in section :numref:`how-to-make-mesh`.

2. Create surface datasets in ``tools/mksurfdata_esmf``. Details in section :numref:`creating-surface-datasets`.


Merge WRF initial conditions into an existing CTSM initial condition file
--------------------------------------------------------------------------

The following procedure is if you'd wish to merge WRF inital conditions from
``wrfinput`` file into CTSM initial condition file ::

    module load ncl
    ncl transfer_wrfinput_to_ctsm_with_snow.ncl 'finidat="the_existing_finidat_file.nc"' 'wrfinput="your_wrfinput_file"' 'merged="the_merged_finidat_file.nc"'

.. todo::

 Versions of the transfer_wrfinput ncl script are available in /glade/work/slevis/git_wrf/ctsm_init/.

