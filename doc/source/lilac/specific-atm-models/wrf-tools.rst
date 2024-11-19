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

1. Create SCRIP grid file from WRF ``geo_em*`` files, using the following ncl
   script::

    ncl create_scrip_file.ncl

   This creates two files that are complements of each other only in the mask field

2. Create mapping files by using ``mkmapdata`` code under
   ``CTSM/tools/mkmapdata/``.

   Using environment variables set the following environment varibales needed
   by ``mkunitymap.ncl`` code::

    setenv GRIDFILE1 wrf2clm_ocean_noneg.nc
    setenv GRIDFILE2 wrf2clm_land_noneg.nc
    setenv MAPFILE wrf2clm_mapping_noneg.nc
    setenv PRINT TRUE

    ncl mkunitymap.ncl

.. warning::

    This will throw some git errors if not run in a repository.

3. Create ESMF mapping files by running ``regridbatch.sh``::

     qsub regridbatch.sh

4. In your ctsm repository directory, build::

     ../../../configure --macros-format Makefile --mpilib mpi-serial

.. todo::
    Update the below, as domain files aren't needed with nuopc.

5. Generate CTSM domain files using ``get_domain`` tool::

     ./gen_domain -m /glade/work/$USER/ctsm/nldas_grid/scrip/wrf2clm_mapping_noneg.nc -o wrf2clm_ocn_noneg -l wrf2clm_lnd_noneg

.. todo::
    Update the below, as ``mksurfdata.pl`` no longer exists.

6. Create surface datasets in ``tools/mksurfdata_esmf``::

     ./mksurfdata.pl -res usrspec -usr_gname "nldas" -usr_gdate "190124" -usr_mapdir "/glade/work/$USER/ctsm/nldas_grid/map" -y 2000 -exedir "/glade/u/home/$USER/src/ctsm/ctsm_surfdata/tools/mksurfdata_esmf" -no-crop

Merge WRF initial conditions into an existing CTSM initial condition file
--------------------------------------------------------------------------

The following procedure is if you'd wish to merget WRF inital conditions from
``wrfinput`` file into CTSM initial condition file ::

    ncl transfer_wrfinput_to_ctsm_with_snow.ncl 'finidat="the_existing_finidat_file.nc"' 'wrfinput="your_wrfinput_file"' 'merged="the_merged_finidat_file.nc"'

.. todo::

 Sam, can you please make the above ncl script available.

