.. include:: ../substitutions.rst

.. _adding-resolution-names:

=========================
 Adding Resolution Names
=========================

If you are adding files for new resolutions which aren't covered in the namelist_definition file -- you'll need to add them in. The list of valid resolutions is in the ``id="res"`` entry in the ``$CTSMROOT/bld/namelist_files/namelist_definition_ctsm.xml`` file. You need to choose a name for your new resolution and simply add it to the comma delimited list of valid_values for the ``id="res"`` entry. The convention for global Gaussian grids is number_of_latitudes x number_of_longitudes. The convention for global finite volume grids is latitude_grid_size x longitude_grid_size where latitude and longitude is measured in degrees. The convention for unstructured HOMME grids is ne<size>np4, where <size> corresponds to the resolution. The higher <size> is the higher the resolution. So for example, ne60np4 is roughly half-degree while ne240np4 is roughly a eighth degree. For regional or single-point datasets the names have a grid size number_of_latitudes x number_of_longitudes followed by an underscore and then a descriptive name such as a City name followed by an abbreviation for the Country in caps. The only hard requirement is that names be unique for different grid files. Here's what the entry for resolutions looks like in the file:
::

   <entry id="res" type="char*30" category="default_settings"
          group="default_settings"
	  valid_values=
	  "512x1024,360x720cru,128x256,64x128,48x96,32x64,8x16,94x192,0.23x0.31,0.9x1.25,1.9x2.5,2.5x3.33,
	  4x5,10x15,5x5_amazon,1x1_tropicAtl,1x1_vacouverCAN,1x1_mexicocityMEX,1x1_asphaltjungleNJ,
	  1x1_brazil,1x1_urbanc_alpha,1x1_numaIA,1x1_smallvilleIA,0.5x0.5,3x3min,5x5min,10x10min,0.33x0.3,
	  ne4np4,ne16np4,ne30np4,ne60np4,ne120np4,ne240np4,wus12,us20,1km-merge-10min">
	  Horizontal resolutions
	  Note: 0.5x0.5, 5x5min, 10x10min, 3x3min, 1km-merge-10min, and 0.33x0.33 are only used for CLM tools
   </entry>

As you can see you just add your new resolution names to the end of the valid_values list.
