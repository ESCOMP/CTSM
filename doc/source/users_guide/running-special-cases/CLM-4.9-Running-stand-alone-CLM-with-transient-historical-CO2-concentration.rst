.. _running-with-historical-co2-forcing:

=====================================
 Running with historical CO2 forcing
=====================================

In this case you want to run a simulation with stand-alone CLM responding to changes in CO2 for a historical period. 
For this example, we will start with the "I_1850-2000_CN" compset that has transient: land-use, Nitrogen and Aerosol deposition already. 
You could also use another compset if you didn't want these other features to be transient. 
In order to get CO2 to be transient we need to add a new streams file and add it to the list of streams in the user_nl_datm file. 
You also need a NetCDF datafile that datm can read that gives the variation. You could supply your own file, but we have a standard file that is used by CAM for this and our example will make use of this file.

.. note:: Most everything here has to do with changing datm rather than CLM to allow this to happen. As such the user that wishes to do this should first become more familiar with datm and read the `CESM Data Model User's Guide <CLM-URL>`_ especially as it pertains to the datm.

.. warning:: This section documents the process for doing something that is non-standard. There may be errors with the documentation and process, and you may have to do some work before all of this works for you. If that is the case, we recommend that you do further research into understanding the process and the files, as well as understanding the datm and how it works. You may have to read documentation found in the code for datm as well as "csm_share".

The datm has "streams" files that have rough XML-like syntax and specify the location and file to get data from, as well as information on the variable names and the data locations of the grid points. 
The datm expects specific variable names and the datm "maps" the expected variable names from the file to the names expected by datm. 
The file we are working with here is a file with a single-point, that covers the entire globe (so the vertices go from -90 to 90 degrees in latitude and 0 to 360 degrees in longitude). 
Since it's a single point it's a little easier to work with than datasets that may be at a given horizontal resolution. 
The datm also expects that variables will be in certain units, and only expects a limited number of variables so arbitrary fields can NOT be exchanged this way. 
However, the process would be similar for datasets that do contain more than one point.

The three things that are needed: a domain file, a data file, and a streams text file. 
The domain file is a CF-compliant NetCDF file that has information on the grid points (latitudes and longitudes for cell-centers and vertices, mask , fraction, and areas). 
The datafile is a CF-compliant NetCDF file with the data that will be mapped. 
The streams text file is the XML-like file that tells datm how to find the files and how to map the variables datm knows about to the variable names on the NetCDF files. Note, that in our case the domain file and the data file are the same file. In other cases, the domain file may be separate from the data file.

First we are going to create a case, and we will edit the ``user_nl_datm`` so that we add a CO2 data stream in. 
There is a streams text file available in ``models/lnd/clm/doc/UsersGuide/co2_streams.txt``, that includes file with a CO2 time-series from 1765 to 2007.


Example: Transient Simulation with Historical CO2
--------------------------------------------------------------
::

   > cd scripts
   > ./create_newcase -case DATM_CO2_TSERIES -res f19_g16 -compset I20TRCRUCLM45BGC 
   > cd DATM_CO2_TSERIES

   # Set CCSM_BGC to CO2A so that CO2 will be passed from atmosphere to land
   # Set CLM_CO2_TYPE to diagnostic so that the land will use the value sent from the atmosphere
   > ./xmlchange CCSM_BGC=CO2A,CLM_CO2_TYPE=diagnostic
   > ./case.setup

   # Create the streams file for CO2
   > cat << EOF >> datm.streams.txt.co2tseries

   <streamstemplate>
      <general_comment>
         This is a streams file to pass historical CO2 from datm8 to the other
	 surface models. It reads in a historical dataset derived from data used
	 by CAM. The getco2_historical.ncl script in models/lnd/clm2/tools/ncl_scripts
	 was used to convert the CAM file to a streams compatible format (adding domain
	 information and making CO2 have latitude/longitude even if only for a single 
	 point.
      </general_comment>
   <stream>
      <comment>
        Input stream description file for historical CO2 reconstruction data

        04 March 2010: Converted to form that can be used by datm8 by Erik Kluzek
        18 December 2009: Prepared by B. Eaton using data provided by 
        Jean-Francois Lamarque. All variables except f11 are directly from
        PRE2005_MIDYR_CONC.DAT. Data from 1765 to 2007 with 2006/2007 just
        a repeat of 2005.
      </comment>
      <dataSource>
         CLMNCEP
      </dataSource>
      <domainInfo>
         <variableNames>
            time    time
            lonc    lon
            latc    lat
            area    area
            mask    mask
         </variableNames>
         <filePath>
            $CSMDATA/atm/datm7/CO2
         </filePath>
         <fileNames>
            fco2_datm_1765-2007_c100614.nc
         </fileNames>
      </domainInfo>
      <fieldInfo>
         <variableNames>
            CO2        co2diag
         </variableNames>
         <filePath>
            $CSMDATA/atm/datm7/CO2
         </filePath>
         <fileNames>
            fco2_datm_1765-2007_c100614.nc
         </fileNames>
      </fieldInfo>
   </stream>
   </streamstemplate>

   EOF

   # And copy it to the run directory
   > cp datm.streams.txt.co2tseries $RUNDIR

   # Run preview namelist so we have the namelist in CaseDocs
   > ./preview_namelists

The first thing we will do is to edit the ``user_nl_datm`` file to add a CO2 file stream in. 
To do this we will copy a ``user_nl_datm`` in with the changes needed. The file ``addco2_user_nl_datm.user_nl`` is in ``models/lnd/clm/doc/UsersGuide`` and looks like this...
::

   dtlimit = 1.5,1.5,1.5,1.5,1.5
   fillalgo = 'nn','nn','nn','nn','nn'
   fillmask = 'nomask','nomask','nomask','nomask','nomask'
   mapalgo = 'bilinear','bilinear','bilinear','bilinear','nn'
   mapmask = 'nomask','nomask','nomask','nomask',nomask'
   streams = "datm.streams.txt.CLM_QIAN.Solar 1895 1948 1972  ", "datm.streams.txt.CLM_QIAN.Precip 1895 1948 1972  ",
             "datm.streams.txt.CLM_QIAN.TPQW 1895 1948 1972  ", "datm.streams.txt.presaero.trans_1850-2000 1849 1849 2006",
	     "datm.streams.txt.co2tseries 1766 1766 2005 "
   taxmode = 'cycle','cycle','cycle','cycle','extend'
   tintalgo = 'coszen','nearest','linear','linear','linear'

You just copy this into your case directory. But, also compare it to the version in ``CaseDocs`` to make sure the changes are just to add in the new CO2 stream. Check to see that filenames, and start, end and align years are correct.
::

   > cp ../../models/lnd/clm/doc/UsersGuide/addco2_user_nl_datm.user_nl user_nl_datm
   > diff user_nl_datm CaseDocs/datm_atm_in

Once, you've done that you can build and run your case normally.

.. warning:: This procedure assumes you are using a ``I20TRCRUCLM45BGC`` compset out of the box, with ``DATM_PRESAERO`` equal to trans_1850-2000. So it assumes standard CLM4.5 CRUNCEP atmosphere forcing, and transient prescribed aerosols from streams files. If your case changes anything here your ``user_nl_datm`` file will need to be adjusted to work with it.

.. note:: The intent of the ``user_nl_datm`` is to add an extra streams file for CO2 to the end of the streams variable, and other arrays associated with streams (adding mapalgo as a new array with bilinear for everything, but the CO2 file which should be "nn" for nearest neighbor). Other variables should be the same as the other stream values.

.. warning:: The streams file above is hard-coded for the path of the file on NCAR computers. To use it on an outside machine you'll need to edit the filepath in the streams file to point to the location where you have the file.

After going through these steps, you will have a case where you have datm reading in an extra streams text file that points to a data file with CO2 data on it that will send that data to the CLM.

