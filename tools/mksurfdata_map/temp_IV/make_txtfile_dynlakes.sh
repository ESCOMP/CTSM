#!/bin/bash

# Script to create the .txt file for every year

file='/glade/u/home/ivanderk/clm5.0/tools/mksurfdata_map/landuse_timeseries_hist_dynlake_simyr1850-2017.txt'

for year in {1850..2017}; do

    # when writing echo message: make sure formatting is right! (201 in total). 
    
	echo "/glade/scratch/ivanderk/inputdata/rawdata/timeseries/mksurf_lake_0.05x0.05_hist_clm5_hydrolakes_${year}_c20200305.nc                                                                                   ${year}" >> $file 

done




