# CTSM Tools Site and Regional Directory

#### $CTSMROOT/tools/site_and_regional/README

The purpose of this directory is to contain all of the scripts that involve creating CTSM input data files
for single site as well as regional cases.

The python scripts require the following settings before running:

``` shell
# Do what's needed to make conda available on your system
../../py_env_create
conda activate ctsm_pylib
```

## Brief description of scripts:

subset_data
    This script extracts domain files, surface dataset, and DATM files
    at either a single point or a region using the global dataset.
    For extracting domain files, surface dataset, and DATM files at a single point, use:
        `./subset_data point`

    For extracting domain files, surface dataset, and DATM files at a region, use:
        `./subset_data region`

mesh_maker
    This script creates a mesh file from a netcdf file with valid lats and lons.

mesh_plotter
    This script plots a mesh file

modify_singlept_site_neon
    After running `subset_data` overwrite some fields with site-specific
    data for neon sites.

run_tower
    Wrapper script for running CTSM simulations for one or more
    NEON or PLUMBER2 sites for spin-up or transient run types.

neon_surf_wrapper
    Wrapper script that run `subset_data` to extract data for all neon points and then
    use modify_singlept_site_neon to update site-specific fields.
    This code uses neon_sites_dompft.csv to determine --dompft (dominant pft types) values.

plumber2_surf_wrapper
    Wrapper script that run `subset_data` to extract data for all PLUMBER2 points.
    This code uses PLUMBER2_sites.csv to determine --dompft (dominant pft types) values.

neon_gcs_upload
    Script to rename and upload NEON site finidat files to NEON server
    for use in transient startup cases

## DEPRECATED SCRIPTS:

mknoocnmap.pl 
        Script to create unity mapping dataset for single-point
        or regional studies over land-only (no ocean).

NCL Scripts available:

mkunitymap.ncl
        NCL script to create a unity map -- ran by above script (mknoocnmap.pl)

