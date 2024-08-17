.. include:: ../substitutions.rst

.. _adding-ptclm-site-data:

============================
Adding PTCLMmkdata Site Data
============================

The "sitegroupname" option to PTCLMmkdata looks for groups of sites in the files in the ``PTCLM_sitedata`` directory under the PTCLMmkdata directory. You can add new names available for this option including your own lists of sites, by adding more files in this directory. There are three files for each "sitegroupname": ``$SITEGROUP_sitedata.txt``, ``$SITEGROUP_soildata.txt`` and ``$SITEGROUP_pftdata.txt`` (where ``$SITEGROUP`` is the name that would be entered as "sitegroupname" to PTCLMmkdata). Each file needs to have the same list of sites, but gives different information: site data, PFT data, and soil data respectively. Although the site codes need to be the same between the three files, the files do NOT have to be in the same order. Each file has a one-line header that lists the contents of each column which are separated by commas. The first column for each of the files is the "site_code" which must be consistent between the three files. The site code can be any unique character string, but in general we use the AmeriFlux site code.

Site data file:`` $SITEGROUP_sitedata.txt``): The header for this file is:
::

   site_code,name,state,lon,lat,elev,startyear,endyear,alignyear

The columns: name, state, and elevation are informational only. Name is a longer descriptive name of the site, and state is the state for U.S. sites or country for non U.S. sites. The columns: lon and lat are the longitude and latitude of the location in decimal degrees. The last three columns are the start and ending year for the data and the align year for an 1850 case for the data. The align year is currently unused.

Soil data file: ``$SITEGROUP_soildata.txt``): The header for this file is:
::

   site_code,soil_depth,n_layers,layer_depth,layer_sand%,layer_clay%

The first three fields after "site_code" are currently unused. The only two that are used are the percent sand and clay columns to set the soil texture.

PFT data file: ``$SITEGROUP_pftdata.txt```): The header for this file is:
::

   site_code,pft_f1,pft_c1,pft_f2,pft_c2,pft_f3,pft_c3,pft_f4,pft_c4,pft_f5,pft_c5

This file gives the vegetation coverage for the different vegetation types for the site. The file only supports up to five PFT's at the same time. The columns with "pft_f" are the fractions for each PFT, and the columns with "pft_c" is the integer index of the given PFT. Look at the pft-physiology file to see what the PFT index for each PFT type is.

----------------------------------------------------
Dynamic Land-Use Change Files for use by PTCLMmkdata
----------------------------------------------------

There is a mechanism for giving site-specific land-use change in PTCLMmkdata. Adding site specific files to the ``PTCLM_sitedata`` directory under PTCLMmkdata allows you to specify the change in vegetation and change in harvesting (for the CN model) for that site. Files are named: ``$SITE_dynpftdata.txt``. There is a sample file for the US-Ha1 site called: ``US-Ha1_dynpftdata.txt``. The file has a one-line header with the information that the file has, and then one-line for each year with a transition. The header line is as follows:
::

   trans_year,pft_f1,pft_c1,pft_f2,pft_c2,pft_f3,pft_c3,pft_f4,pft_c4,pft_f5,pft_c5,har_vh1,har_vh2,har_sh1,har_sh2,har_sh3,graze,hold_harv,hold_graze

This file only requires a line for each year where a transition or harvest happens. As in the "pftdata" file above "pft_f" refers to the fraction and "pft_c" refers to the PFT index, and only up to five vegetation types are allowed to co-exist. The last eight columns have to do with harvesting and grazing. The last two columns are whether to hold harvesting and/or grazing constant until the next transition year and will just be either 1 or 0. This file will be converted by the **PTCLM_sitedata/cnvrt_trnsyrs2_pftdyntxtfile.pl** script in the PTCLMmkdata directory to a format that **mksurfdata_esmf** can read that has an entry for each year for the range of years valid for the compset in question.

.. _converting-ameriflux-for-ptclmmkdata:

------------------------------------------------
Converting AmeriFlux Data for use by PTCLMmkdata
------------------------------------------------

AmeriFlux data comes in comma separated format and is available from: `http://public.ornl.gov/ameriflux/dataproducts.shtml <http://public.ornl.gov/ameriflux/dataproducts.shtml>`_. Before you download the data you need to agree to the usage terms.

Here is a copy of the usage terms from the web-site on June/13/2011.

"The AmeriFlux data provided on this site are freely available and were furnished by individual AmeriFlux scientists who encourage their use. Please kindly inform the appropriate AmeriFlux scientist(s) of how you are using the data and of any publication plans. Please acknowledge the data source as a citation or in the acknowledgments if the data are not yet published. If the AmeriFlux Principal Investigators (PIs) feel that they should be acknowledged or offered participation as authors, they will let you know and we assume that an agreement on such matters will be reached before publishing and/or use of the data for publication. If your work directly competes with the PI's analysis they may ask that they have the opportunity to submit a manuscript before you submit one that uses unpublished data. In addition, when publishing, please acknowledge the agency that supported the research. Lastly, we kindly request that those publishing papers using AmeriFlux data provide preprints to the PIs providing the data and to the data archive at the Carbon Dioxide Information Analysis Center (CDIAC)."

The above agreement applies to the "US-UMB" dataset imported into our repository as well, and Gil Bohrer is the PI on record for that dataset.

The CESM can NOT handle missing data, so we recommend using the "Level 4" Gap filled datasets. The fields will also need to be renamed. The "WS" column becomes "WIND", "PREC" becomes "PRECmms", "RH" stays as "RH", "TA" becomes "TBOT", "Rg" becomes "FSDS", "Rgl" becomes "FLDS", "PRESS" becomes "PSRF". "ZBOT" can just be set to the constant of "30" (m). The units of Temperature need to be converted from "Celsius" to "Kelvin" (use the value in ``SHR_CONST_TKFRZ`` in the file ``models/csm_share/shr/shr_const.F90`` of ``273.15``. The units of Pressure also need to be converted from "kPa" to "Pa". LATIXY, and LONGXY should also be set to the latitude and longitude of the site.

-----------------------------------------------------------------
Example: PTCLMmkdata transient example over a shorter time period
-----------------------------------------------------------------

This is an example of using PTCLMmkdata for Harvard Forest (AmeriFlux site code US-Ha1) for transient land use 1991-2006. In order to do this we would've needed to have converted the AmeriFlux data into NetCDF format as shown in :ref:`converting-ameriflux-for-ptclmmkdata` section above. Also note that this site has a site-specific dynamic land-use change file for it ``PTCLM_sitedata/US-Ha1_dynpftdata.txt`` in the PTCLMmkdata directory and this file will be used for land-use change and harvesting rather than the global dataset.

::

     > cd $CTSMROOT/tools/PTCLM
     # We are going to use forcing data over 1991 to 2006, but we need to start with
     # a transient compset to do so, so we use the 20th Century transient: 1850-2000
     # Note: When creating the fpftdyn dataset for this site it will use the
     #     PTCLM_sitedata/US-Ha1_dynpftdata.txt
     # file for land-use change and harvesting
     > ./PTCLMmkdata -s US-Ha1 -d $MYCSMDATA --sitegroupname AmeriFlux
     > mkdir $MYCSMDATA/atm/datm7/CLM1PT_data/1x1pt_US-Ha1
     > cd $MYCSMDATA/atm/datm7/CLM1PT_data/1x1pt_US-Ha1
     # Copy data in NetCDF format to this directory, filenames should be YYYY-MM.nc
     # The fieldnames on the file should be:
     #    FLDS,FSDS,LATIXY,   LONGXY,   PRECTmms,PSRF,RH,TBOT,WIND,ZBOT
     # With units
     #    W/m2,W/m2,degrees_N,degrees_E,mm/s,    Pa,  %, K,   m/s, m
     # The time coordinate units should be: days since YYYY-MM-DD 00:00:00
     > cd ../../../../../US-Ha1_I20TRCRUCLM45BGC
     # Now we need to set the start date to 1991, and make sure the align year is for 1991
     > ./xmlchange RUN_STARTDATE=1991-01-01,DATM_CLMNCEP_YR_ALIGN=1991
     # Similarly for Nitrogen deposition data we cycle over: 1991 to 2006
     > cat << EOF >> user_nl_clm
     model_year_align_ndep=1991,stream_year_first_ndep=1991,stream_year_last_ndep=2006
     EOF
