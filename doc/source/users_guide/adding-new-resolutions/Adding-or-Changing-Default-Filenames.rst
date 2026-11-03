.. include:: ../substitutions.rst

.. _changing-default-filenames:

============================
 Changing Default Filenames
============================

To add or change the default filenames you edit the ``$CTSMROOT/bld/namelist_files/namelist_definition_ctsm.xml`` and either change an existing filename or add a new one. Most entries in the default namelist files, include different attributes that describe the different properties that describe the differences in the datasets. Attributes include the resolution, year to simulation, range of years to simulate for transient datafiles, the land-mask, the representative concentration pathway (RCP) for future scenarios, and the type of biogeochemistry (BGC) model used. For example the surface dataset (``fsurdat``) for the 1.9x2.5 resolution is as follows:

::
  
  <fsurdat hgrid="0.9x1.25"  sim_year="1850" use_crop=".true." >
  lnd/clm2/surfdata_map/surfdata_0.9x1.25_78pfts_CMIP6_simyr1850_c170824.nc
  </fsurdat>


Other ``fsurdat`` files are distinguished from this one by their resolution (``hgrid``), simulation year (``sim_year``) and prognostic crop (``use_crop``) attributes.

