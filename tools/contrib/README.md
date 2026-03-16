# $CTSMROOT/tools/contrib/README                                                    Jan/24/2019

The purpose of this directory is for users of CTSM to contribute scripts for pre or post processing or
case management of CTSM that others might find useful.  The script should have some documentation made 
available before adding it. These scripts may not be as well tested or supported as other CTSM
tools. They are also ONLY assumed to work on the NCAR supercomputer. So paths will be hardwired to
assume NCAR directory structures.

The python scripts require the following settings before running on Derecho:

``` shell
module load conda
../../py_env_create
conda activate ctsm_pylib
```

Brief description of scripts:

abm_raw.ncl
        Handle the peak crop-fire month for CMIP7
add_tillage_to_paramsfile.py
        Add tillage data to the parameter file
        Can this be deleted since we have modify paramfile scripts? EBK 3/12/2026
create_scrip_file.ncl
        Create a SCRIP grid file needed for running with WRF
CRUJRA_antarctica.ipynb
        Jupyter notebook to add forcing data over Antarctica for the CRUJRA forcing
CRUJRA_greenland.ipynb
        Jupyter notebook to add forcing data over Greenland for the CRUJRA forcing
run_clm_historical.v11.csh
	does all the setup and submission required to do a 1850-2023 CLM 
	historical simulation in five separate submissions
        v11 - Oleson, 07/2025

modify_singlept_site
        Modify some data on a surface dataset created by site_and_regional/subset_data

popden.ncl
        Script to modify the CMIP7 population density file for use by clm6_0

SpinupStability_SP_v10.ncl
        This script assesses the equilibrium state of a Satellite Phenology (SP) 
        spinup run, works on either monthly or annual mean history files - Keith
        Oleson 07/2025

SpinupStability_BGC_v11.ncl
        This script assesses the equilibrium state of a Biogeochemistry (BGC) 
        spinup run, works on either monthly or annual mean history files - Keith
        Oleson 07/2025

SpinupStability_BGC_v12_SE.ncl
        This script assesses the equilibrium state of a ne30pg3 (spectral element)
        Biogeochemistry (BGC) spinup run, works on either monthly or annual mean 
        history files - Oleson 07/2025

remove_duplicate_tests.py
        Script to rewrite the testlist removing duplicates
        Should this be removed? EBK 3/12/2026

run_clmtowers
        This script will run any number of flux tower sites.
        It's based on having created surface datasets with PTCLM.
        v1 - Keith Oleson, 8/2015

ssp_anomaly_forcing_smooth
        This script creates anomaly forcing for CMIP6 SSP scenarios that 
        can be used to run CTSM in CESM with datm.
        v0 -- Sean Swenson
        v1 - Peter Lawrence 3/2020
        v2 - Sean Swenson/Erik Kluzek 6/2022

test_rxcropmaturity_python.sh
        Test the prescribed crop maturity python tools

tweak_latlons.py
        Tweak the latitude/longitudes on a file so that a nearest neighbor mapping
        with ESMF before ESMF8.9.0 won't differ on processor count


