.. include:: ../substitutions.rst

.. _creating-maps-for-mksurfdata:

*********************************************
Creating input for surface dataset generation
*********************************************

Generating SCRIP grid files
==================================

The utility ``mkmapdata.sh`` requires SCRIP format input files to describe the input and output grids that maps are generated for. CLM provides a utility, ``mkmapgrids`` that generates those files. The program converts old formats of CAM or CLM grid files to SCRIP grid format. There is also a NCL script (``mkscripgrid.ncl``) to create regular latitude longitude regional or single-point grids at the resolution the user desires.

SCRIP grid files for all the standard model resolutions and the raw surface datasets have already been done and the files are in the XML database. Hence, this step doesn't need to be done -- EXCEPT WHEN YOU ARE CREATING YOUR OWN GRIDS.

.. _using-mkocnmap:

Using mknocnmap.pl to create grid and maps for single-point regional grids
--------------------------------------------------------------------------

.. todo::
    Update the below, as domain files aren't needed with nuopc.

If you want to create a regular latitude/longitude single-point or regional grid, we suggest you use ``mknoocnmap.pl`` in ``$CTSMROOT/tools/mkmapdata`` which will create both the SCRIP grid file you need (using ``$CTSMROOT/tools/mkmapgrids/mkscripgrid.ncl``) AND an identity mapping file assuming there is NO ocean in your grid domain. If you HAVE ocean in your domain you could modify the mask in the SCRIP grid file for ocean, and then use ``ESMF_RegridWeightGen`` to create the mapping file, and ``gen_domain`` to create the domain file. Like other tools, ``./mkmapdata/mknoocnmap.pl`` has a help option with the following:
::

   SYNOPSIS
	mknoocnmap.pl [options]	    Gets map and grid files for a single land-only point.
   REQUIRED OPTIONS
	-centerpoint [or -p] <lat,lon> Center latitude,longitude of the grid to create.
	-name [-or -n] <name>	    Name to use to describe point

   OPTIONS
	-dx <number>                   Size of total grid in degrees in longitude direction
				       (default is 0.1)
	-dy <number>                   Size of total grid in degrees in latitude direction
				       (default is 0.1)
	-silent [or -s]		    Make output silent
	-help [or -h]		    Print usage to STDOUT.
	-verbose [or -v]		    Make output more verbose.
	-nx <number>                   Number of longitudes (default is 1)
	-ny <number>                   Number of latitudes  (default is 1)

See :numref:`Figure mknoocnmap.pl` for a visual representation of this process.

Creating mapping files for mksurfdata_esmf
==============================================

``mkmapdata.sh`` uses the above SCRIP grid input files to create SCRIP mapping data files (uses ESMF).

The bash shell script ``$CTSMROOT/tools/mkmapgrids/mkmapdata.sh`` uses ``ESMF_RegridWeightGen`` to create a list of maps from the raw datasets that are input to ``mksurfdata_esmf``. Each dataset that has a different grid, or land-mask needs a different mapping file for it, but many different raw datasets share the same grid/land-mask as other files. Hence, there doesn't need to be a different mapping file for EACH raw dataset---just for each raw dataset that has a DIFFERENT grid or land-mask. See :numref:`Figure mkmapdata.sh` for a visual representation of how this works. The bash script figures out which mapping files it needs to create and then runs ``ESMF_RegridWeightGen`` for each one. You can then either enter the datasets into the XML database (see Chapter :numref:`adding-new-resolutions-section`), or leave the files in place and use the ``-res usrspec -usr_gname -usr_gdate`` options to ``mksurfdata_esmf``. ``mkmapdata.sh`` has a help option with the following
::

   ../../tools/mkmapdata/mkmapdata.sh

   **********************
   usage on cheyenne:Figure mkmapdata.sh
   ./mkmapdata.sh

   valid arguments:
   [-f|--gridfile <gridname>]
	Full pathname of model SCRIP grid file to use
	This variable should be set if this is not a supported grid
	This variable will override the automatic generation of the
	filename generated from the -res argument
	the filename is generated ASSUMING that this is a supported
	grid that has entries in the file namelist_defaults_clm.xml
	the -r|--res argument MUST be specied if this argument is specified
   [-r|--res <res>]
	Model output resolution (default is 10x15)
   [-t|--gridtype <type>]
	Model output grid type
	supported values are [regional,global], (default is global)
   [-b|--batch]
	Toggles batch mode usage.
   If you want to run in batch mode
	you need to have a separate batch script for a supported machine
	that calls this script interactively - you cannot submit this
	script directory to the batch system
   [-l|--list]
	List mapping files required (use check_input_data to get them)
	also writes data to clm.input_data_list
   [-d|--debug]
	Toggles debug-only (don't actually run mkmapdata just echo what would happen)
   [-h|--help]
	Displays this help message
   [-v|--verbose]
	Toggle verbose usage -- log more information on what is happening

    You can also set the following env variables:
     ESMFBIN_PATH - Path to ESMF binaries
		    (default is /contrib/esmf-5.3.0-64-O/bin)
     CSMDATA ------ Path to CESM input data
		    (default is /glade/p/cesm/cseg/inputdata)
     MPIEXEC ------ Name of mpirun executable
		    (default is mpirun.lsf)
     REGRID_PROC -- Number of MPI processors to use
		    (default is 8)

   **pass environment variables by preceding above commands
     with 'env var1=setting var2=setting '
   **********************

.. _Figure mkmapdata.sh:

.. figure:: mkmapdata_details.jpeg

  Details of running mkmapdata.sh

Each of the raw datasets for ``mksurfdata_esmf`` needs a mapping file to map from the output grid you are running on to the grid and land-mask for that dataset. This is what ``mkmapdata.sh`` does. To create the mapping files you need a SCRIP grid file to correspond with each resolution and land mask that you have a raw data file in ``mksurfdata_esmf``. Some raw datasets share the same grid and land mask -- hence they can share the same SCRIP grid file. The output maps created here go into ``mksurfdata_esmf`` see :numref:`Figure Workflow of CLM5 Land Use Data Tool and mksurfdata_esmf Tool`.
