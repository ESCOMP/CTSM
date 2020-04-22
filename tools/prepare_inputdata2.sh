#!/bin/bash

#SBATCH --account=nn2806k 
#SBATCH --job-name=mkmapdata
#SBATCH --mem-per-cpu=256G --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --time=20:00:00

date="200422"
plotlat=(61.0243 60.8231 60.8328 60.9335 60.8203 60.8760 61.0866 60.5445 61.0355 60.8803 60.6652 60.6901)
plotlon=(8.12343 7.27596 7.17561 6.41504 8.70466 7.17666 6.63028 6.51468 9.07876 7.16982 6.33738 5.96487)
plotname=(ALP1 ALP2 ALP3 ALP4 SUB1 SUB2 SUB3 SUB4 BOR1 BOR2 BOR3 BOR4)


for i in 7 8 9 10 11
do

######## Make script grids:
module purge
module load NCL/6.6.2-intel-2019b
cd ~/ctsm/tools/mkmapdata
mkdir -p /cluster/shared/noresm/inputdata/share/scripgrids/fates_platform/${plotname[i]}

echo $i
echo ${plotlat[$i]} ${plotlon[$i]}
./mknoocnmap.pl -centerpoint ${plotlat[i]},${plotlon[i]} -name ${plotname[i]} -dx 0.01 -dy 0.01
mv ~/ctsm/tools/mkmapgrids/*.nc /cluster/shared/noresm/inputdata/share/scripgrids/fates_platform/${plotname[i]}/

#ALP1_noocean_c200417.nc SCRIPgrid_ALP1_nomask_c200417.nc map_ALP1_noocean_to_ALP1_nomask_aave_da_200417.nc  

######## Make domain file:
module purge
cd ~/ctsm/cime/tools/mapping/gen_domain_files
mkdir -p /cluster/shared/noresm/inputdata/share/domains/fates_platform/${plotname[i]}

#### Compile 
#cd src/
#../../../configure --macros-format Makefile --mpilib mpi-serial --machine saga
#. ./.env_mach_specific.sh ; gmake

#### Run
. ./src/.env_mach_specific.sh
./gen_domain -m /cluster/shared/noresm/inputdata/share/scripgrids/fates_platform/${plotname[i]}/map_${plotname[i]}_noocean_to_${plotname[i]}_nomask_aave_da_${date}.nc -o ${plotname[i]} -l ${plotname[i]}

mv domain.lnd.${plotname[i]}_${plotname[i]}.${date}.nc domain.lnd.${plotname[i]}.${date}.nc
mv domain*.nc /cluster/shared/noresm/inputdata/share/domains/fates_platform/${plotname[i]}/

#domain.ocn.1x1_APL1_ocn.200417.nc domain.lnd.1x1_APL1_1x1_APL1_ocn.200417.nc domain.ocn.1x1_APL1_1x1_APL1_ocn.200417.nc

######## Make mapping file:
module purge
cd ~/ctsm/tools/mkmapdata
mkdir -p /cluster/shared/noresm/inputdata/lnd/clm2/mappingdata/maps/fates_platform/${plotname[i]}/
#module load ESMF/8.0.0-intel-2019b
#module load NCO/4.9.1-intel-2019b 
#module load NCL/6.6.2-intel-2019b
#### Modify regridbatch.sh
#### Run regridbatch.sh
./regridbatch.sh 1x1_${plotname[i]} /cluster/shared/noresm/inputdata/share/scripgrids/fates_platform/${plotname[i]}/SCRIPgrid_${plotname[i]}_nomask_c${date}.nc 
mv map_* /cluster/shared/noresm/inputdata/lnd/clm2/mappingdata/maps/fates_platform/${plotname[i]}/

######## Make surface data file:
cd ~/ctsm/tools/mksurfdata_map
mkdir -p /cluster/shared/noresm/inputdata/lnd/clm2/surfdata_map/fates_platform/${plotname[i]}

#### Compile
##module load netCDF-Fortran/4.5.2-iimpi-2019b
##Modify Makefile.common
##gmake clean
##gmake

#### Run
# ./mksurfdata.pl -no-crop -res usrspec -usr_gname 1x1_km -usr_gdate 200417 -usr_mapdir /cluster/shared/noresm/inputdata/lnd/clm2/mappingdata/maps/fates_platform/ALP1 -dinlc /cluster/shared/noresm/inputdata -hirespft -years "2000" -allownofile
##** trouble getting vegtyp file with:
##Need to modify ~/ctsm/bld/namelist_files/namelist_defaults_clm4_5_tools.xml to have vegtyp file for 2000
## Or change the -years to 2005

./mksurfdata.pl -no-crop -res usrspec -usr_gname 1x1_${plotname[i]} -usr_gdate ${date} -usr_mapdir /cluster/shared/noresm/inputdata/lnd/clm2/mappingdata/maps/fates_platform/${plotname[i]} -dinlc /cluster/shared/noresm/inputdata -hirespft -years "2005" -allownofile

mv surfdata_1x1_${plotname[i]}_hist_16pfts_Irrig_CMIP6_simyr2005_c${date}.nc surfdata_${plotname[i]}_simyr2000.nc

mv surfdata*${plotname[i]}* /cluster/shared/noresm/inputdata/lnd/clm2/surfdata_map/fates_platform/${plotname[i]}


######## Create a case
#cd ~/ctsm/cime/scripts

#./create_newcase --case ../../../ctsm_cases/2000FATES_seedclim_${plotname[i]} --compset 2000_DATM%1PT_CLM50%FATES_SICE_SOCN_MOSART_SGLC_SWAV --res CLM_USRDAT --machine saga --run-unsupported --project nn2806k

#./create_clone --case ../../../ctsm_cases/2000FATES_seedclim_APL1 --clone ../../../ctsm_cases/2000FATES_seedclim1 --project nn2806k

######## Build the case
#### Need to modify  ~/ctsm/bld/namelist_files/namelist_defaults_usr_files.xml to point to the correct surface data 
#./case.build

######## Run the case
#### Need to modify  ~/ctsm/cime/src/components/data_comps/datm/cime_config/namelist_definition_datm.xml to use specific atm forcing. 

done

