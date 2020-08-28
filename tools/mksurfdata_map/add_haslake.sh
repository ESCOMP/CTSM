#!/bin/bash

# Script to add ISLAKE variable to landuse.timeseries .nc file used for input. 
# Right now, this happens manually. Ideally should be incorporated in mksurfdat script. 

#dir=/glade/scratch/ivanderk/inputdata/clm_input/transRES/

#cd $dir


infile=landuse.timeseries_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr1850-2015_c190715.nc
#infile=landuse.timeseries_0.9x1.25_hist_16pfts_Irrig_CMIP6_simyr1850-2015_c190722.nc

#define number of years
nyears=166 

cdo -setname,HASLAKE -gtc,0 -selname,PCT_LAKE -seltimestep,$nyears $infile outfile.nc 
cdo merge $infile outfile.nc infile_copy.nc
mv infile_copy.nc $infile
rm outfile.nc 
ncatted -O -a long_name,HASLAKE,o,c,"Boolean all lake grid cells (including dynamical lake cells)" $infile
