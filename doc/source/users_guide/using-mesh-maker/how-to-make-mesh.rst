.. include:: ../substitutions.rst

.. _how-to-make-mesh:

===============================================
 Creating an ESMF mesh file from a netCDF file
===============================================

This page includes instructions for using the  ``mesh_maker`` tool to create a mesh file from a netCDF file with valid 1D or 2D latitude and longitude coordinates. It also shows how to use ``mesh_plotter`` to visualize a mesh file.

.. note:: An **ESMF mesh file** is a netCDF file that includes the information about the grid's coordinates and their connectivity to each other in an **Unstructured Grid Format**. Additional information about ESMF mesh files is available `here <https://earthsystemmodeling.org/docs/release/ESMF_8_0_1/ESMF_refdoc/node3.html#SECTION03028200000000000000>`_.

You can check out the ``mesh_maker`` options like so:

::

   > tools/site_and_regional/mesh_maker --help
   
   |------------------------------------------------------------------|
   |---------------------  Instructions  -----------------------------|
   |------------------------------------------------------------------|
   This script creates ESMF unstructured GRID (mesh file) from a netCDF
   file with valid lats and lons. Provided lats and lons can be 1D or 2D.

   For example for running WRF-CTSM cases, the user can create a mesh
   file for their domain :
       ./mesh_maker.py --input wrfinput_d01 --output my_region
           --lat XLAT --lon XLONG --verbose

   optional arguments:
     -h, --help        show this help message and exit
     --input INPUT     Netcdf input file for creating ESMF mesh.
     --output OUTPUT   Name of the ESMF mesh created.
     --outdir OUT_DIR  Output directory (only if name of output mesh is not
                       defined)
     --lat LAT_NAME    Name of latitude varibale on netCDF input file. If none
                       given, looks to find variables that include 'lat'.
     --lon LON_NAME    Name of latitude varibale on netCDF input file. If none
                       given, looks to find variables that include 'lon'.
     --mask MASK_NAME  Name of mask varibale on netCDF input file. If none given,
                       create a fake mask with values of 1.
     --area AREA_NAME  Name of area variable on netCDF input file. If none given,
                       ESMF calculates element areas automatically.
     --overwrite       If meshfile exists, overwrite the meshfile.
     -v, --verbose     Increase output verbosity

Example: Making and visualizing a mesh file
-------------------------------------------

In this example, we will use ``mesh_maker`` to create a mesh file from a netCDF file with 2D latitudes and longitudes. On the sample input provided, those coordinates are saved on the ``LATIXY`` and ``LONGXY`` variables, respectively.

::

   input_file="python/ctsm/test/testinputs/surfdata_5x5_amazon_hist_78pfts_CMIP6_2000_c230517.nc"
   output_file="meshfile_5x5_amazon.nc"
   
   # Create the file. (Add --verbose for additional debugging information.)
   tools/site_and_regional/mesh_maker --input "${input_file}" --output "${output_file}" --lon LONGXY --lat LATIXY
   
   # Visualize the meshes
   tools/site_and_regional/mesh_plotter --input "${output_file}"

This produces two figures:

.. figure:: test_c240918_regional.png

.. figure:: test_c240918_global.png
