.. include:: ../substitutions.rst

.. _running-with-historical-co2-forcing:

=====================================
 Running with historical CO2 forcing
=====================================

In this case you want to run a simulation with stand-alone CLM responding to changes in CO2 for a historical period. For this example, we will start with the "I_1850-2000_CN" compset that has transient: land-use, Nitrogen and Aerosol deposition already. You could also use another compset if you didn't want these other features to be transient. In order to get CO2 to be transient we need to add a new streams file and add it to the list of streams in the user_nl_datm file. You also need a NetCDF datafile that datm can read that gives the variation. You could supply your own file, but we have a standard file that is used by CAM for this and our example will make use of this file.

.. note:: Most everything here has to do with changing datm rather than CLM to allow this to happen. As such the user that wishes to do this should first become more familiar with datm and read the `CESM Data Model User's Guide <https://esmci.github.io/cime/versions/cesm2.2/html/data_models/index.html>`_ especially as it pertains to the datm.

.. warning:: This section documents the process for doing something that is non-standard. There may be errors with the documentation and process, and you may have to do some work before all of this works for you. If that is the case, we recommend that you do further research into understanding the process and the files, as well as understanding the datm and how it works. You may have to read documentation found in the code for datm as well as "csm_share".

The datm has "streams" files that have rough XML-like syntax and specify the location and file to get data from, as well as information on the variable names and the data locations of the grid points. The datm expects specific variable names and the datm "maps" the expected variable names from the file to the names expected by datm. The file we are working with here is a file with a single-point, that covers the entire globe (so the vertices go from -90 to 90 degrees in latitude and 0 to 360 degrees in longitude). Since it's a single point it's a little easier to work with than datasets that may be at a given horizontal resolution. The datm also expects that variables will be in certain units, and only expects a limited number of variables so arbitrary fields can NOT be exchanged this way. However, the process would be similar for datasets that do contain more than one point.

The three things that are needed: a domain file, a data file, and a streams text file. The domain file is a CF-compliant NetCDF file that has information on the grid points (latitudes and longitudes for cell-centers and vertices, mask, fraction, and areas). The datafile is a CF-compliant NetCDF file with the data that will be mapped. The streams text file is the XML-like file that tells datm how to find the files and how to map the variables datm knows about to the variable names on the NetCDF files. Note, that in our case the domain file and the data file are the same file. In other cases, the domain file may be separate from the data file.

First we are going to create a case, and we will edit the ``user_nl_datm`` so that we add a CO2 data stream in. There is a streams text file available in ``$CTSMROOT/doc/UsersGuide/co2_streams.txt``, that includes file with a CO2 time-series from 1765 to 2007.

Example: Transient Simulation with Historical CO2
--------------------------------------------------------------
::

   > cd scripts
   > ./create_newcase -case DATM_CO2_TSERIES -res f19_g17_gl4 -compset IHistClm50BgcCrop
   > cd DATM_CO2_TSERIES

   # Historical CO2 will already be setup correctly for this compset
   # to check that look at the variables: CCSM_BGC,CLM_CO2_TYPE, and DATM_CO2_TSERIES
   > ./xmlquery CCSM_BGC,CLM_CO2_TYPE,DATM_CO2_TSERIES
   # Expect: CCSM_BGC=CO2A,CLM_CO2_TYPE=diagnostic,DATM_CO2_TSERIES=20tr
   > ./case.setup

   # Run preview namelist so we have the namelist in CaseDocs
   > ./preview_namelists

Once, you've done that you can build and run your case normally.

