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

----------------------------
What are the required files?
----------------------------

Different types of simulations and different types of configurations for CLM require different lists of files. The |version|-BGC or Carbon Nitrogen (cn) Biogeochemistry model for example requires ``stream_fldfilename_ndep`` files, which are NOT required by CLMSP. Transient simulations also require transient datasets, and the names of these datasets are sometimes different from the static versions (sometimes both are required as in the dynamic PFT cases).

In the following table we list the different files used by CLM, they are listed in order of importance, dependencies, and customizing. So the required files are all near the top, and the files used only under different conditions are listed later, and files with the fewest dependencies are near the top, as are the files that are least likely to be customized.

.. _reqd-files-table:

Table 3-1. Required Files for Different Configurations and Simulation Types
---------------------------------------------------------------------------
.. todo::
    Insert table 3-1
