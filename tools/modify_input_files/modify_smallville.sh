#!/bin/bash

# Script that prepares landuse files for smallvilleIA.
# Load the nco module and run the script in the ctsm_pylib environment:
module load nco

# This script runs from the mksurfdata_esmf/Makefile.
# When running standalone, it may need "subset_data_single_point/" in front
# of each landuse.timeseries file name.
 file_to_2100="landuse.timeseries_1x1_smallvilleIA_SSP2-4.5_1850-2100_78pfts_c$(date +%y%m%d).nc"
 file_to_1855="landuse.timeseries_1x1_smallvilleIA_SSP2-4.5_1850-1855_78pfts_c$(date +%y%m%d).nc"
 file_lake="landuse.timeseries_1x1_smallvilleIA_SSP2-4.5_1850-1855_78pfts_dynLakes_c$(date +%y%m%d).nc"
 file_urban="landuse.timeseries_1x1_smallvilleIA_SSP2-4.5_1850-1855_78pfts_dynUrban_c$(date +%y%m%d).nc"
 file_pft="landuse.timeseries_1x1_smallvilleIA_SSP2-4.5_1850-1855_78pfts_dynPft_c$(date +%y%m%d).nc"

# Trim the file to just the years 1850-1855
ncks -d time,0,5 $file_to_2100 $file_to_1855

# Replace all values in the LAKE and CROP variables
ncap2 -s "PCT_LAKE=array(0.,0.,PCT_CROP); PCT_LAKE={0.,50.,25.,25.,25.,25.} ; PCT_LAKE_MAX=array(50.,50.,PCT_CROP_MAX); PCT_CROP=array(0.,0.,PCT_LAKE); PCT_CROP={0.,25.,12.,12.,12.,12.}; PCT_CROP_MAX=array(25.,25.,PCT_LAKE_MAX)" $file_to_1855 $file_lake

# Replace all values in the URBAN and CROP variables
ncap2 -s "PCT_URBAN=array(0.,0.,PCT_URBAN); PCT_URBAN={0.,0.,0.,20.,15.,0.,10.,8.,0.,10.,8.,0.,10.,8.,0.,10.,8.,0.} ; PCT_URBAN_MAX=array(0.,0.,PCT_URBAN_MAX); PCT_URBAN_MAX={20.,15.,0.}; PCT_CROP=array(0.,0.,PCT_LAKE); PCT_CROP={0.,25.,12.,12.,12.,12.}; PCT_CROP_MAX=array(25.,25.,PCT_LAKE_MAX)" $file_to_1855 $file_urban

# Update values in the pft, cft, harvest, and grazing variables as posted here:
# https://github.com/ESCOMP/CTSM/issues/1673#issuecomment-1879156989
ncap2 -s "PCT_NAT_PFT=array(0.,0.,PCT_NAT_PFT); PCT_NAT_PFT={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,100.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,100.,0.,100.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,100.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,50.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,50.,0.,25.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,75.,0.} ; PCT_NAT_PFT_MAX=array(0.,0.,PCT_NAT_PFT_MAX); PCT_NAT_PFT_MAX={100.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,100.,0.}; PCT_CFT=array(0.,0.,PCT_CFT); PCT_CFT={100.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,100.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,91.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,91.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,2.,4.,4.,6.,6.,8.,8.,10.,10.,42.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,4.,4.,4.,4.,4.,4.,4.,4.,4.,64.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; PCT_CFT_MAX=array(0.,0.,PCT_CFT_MAX); PCT_CFT_MAX={100.,2.,2.,3.,3.,4.,4.,5.,5.,91.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; PCT_CROP=array(0.,0.,PCT_CROP); PCT_CROP={0.,0.,100.,100.,50.,25.}; PCT_CROP_MAX=array(100.,100.,PCT_CROP_MAX); HARVEST_SH1=array(0.,0.,HARVEST_SH1); HARVEST_SH1={0.,0.,0.,0.,0.,0.}; HARVEST_SH2=array(0.,0.,HARVEST_SH2); HARVEST_SH2={0.,0.,0.,0.,0.,0.}; HARVEST_SH3=array(0.,0.,HARVEST_SH3); HARVEST_SH3={0.,0.,0.,0.,0.,0.}; HARVEST_VH1=array(0.,0.,HARVEST_VH1); HARVEST_VH1={0.,0.,0.,0.,0.,0.}; HARVEST_VH2=array(0.,0.,HARVEST_VH2); HARVEST_VH2={0.,0.,0.,0.,0.,0.}; GRAZING=array(0.,0.,GRAZING); GRAZING={0.,0.,0.,0.,0.,0.}" $file_to_1855 $file_pft

exit
