# Creating an ESMF mesh file from a netCDF file

This gist includes instructions for creating and visualizing a mesh file from a netcdf file with valid 1D or 2D lats and lons coordinates.

* **ESMF Mesh file** aka **Unstructured Grid File Format** is a netcdf file (format) that includes the information about the grids coordinates and their connectivity to each other. 

Additional information about ESMF mesh files are available [here](https://earthsystemmodeling.org/docs/release/ESMF_8_0_1/ESMF_refdoc/node3.html#SECTION03028200000000000000).

------

In this example, we will use `./mesh_maker.py` which uses `mesh_type.py` to create a mesh file and visualize it.

1- First clone my fork and branch that includes these capabilities:
``` Shell
git clone https://github.com/negin513/ctsm.git ctsm_mesh
cd ctsm_mesh

git checkout subset_mesh_dask
```

2- Next run mesh_maker.py for a netcdf file:

```
cd tools/site_and_regional
```
Check all the avaialble options:

```
./mesh_maker.py --help
```

The output shows all available options for this script:
```
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This script creates ESMF unstructured GRID (mesh file) from a netcdf
file with valid lats and lons. Provided lats and lons can be 1D or 2D.

For example for running WRF-CTSM cases, the user can create a mesh
file for their domain :
    ./mesh_maker.py --input wrfinput_d01 --output my_region
        --lat XLAT --lon XLONG --verbose

optional arguments:
  -h, --help        show this help message and exit
  --input INPUT     Netcdf input file for creating ESMF mesh.
  --output OUTPUT   Name of the ESMF mesh created.
  --outdir OUT_DIR  Output directory (only if name of output mesh is not
                    defined)
  --lat LAT_NAME    Name of latitude varibale on netcdf input file. If none
                    given, looks to find variables that include 'lat'.
  --lon LON_NAME    Name of latitude varibale on netcdf input file. If none
                    given, looks to find variables that include 'lon'.
  --mask MASK_NAME  Name of mask varibale on netcdf input file. If none given,
                    create a fake mask with values of 1.
  --area AREA_NAME  Name of area variable on netcdf input file. If none given,
                    ESMF calculates element areas automatically.
  --overwrite       If meshfile exists, overwrite the meshfile.
  -v, --verbose     Increase output verbosity
  
  ```
  
Let's create a mesh file from a netcdf file with 1D lats and lons. On the sample files provided 1D lat and long coordinates are saved on `lsmlat` and `lsmlon` variables. 

```
./mesh_maker.py --input /glade/scratch/negins/example_files/surfdata_4x5_hist_78pfts_CMIP6_simyr1850_275.0-330.0_-40-15_c220705.nc --output test_mesh_1d.nc --lat lsmlat --lon lsmlon --overwrite
```
`--verbose` option also provide additional information for debugging. 

This script will create regional and global mesh plots. For example for the above files, the plos are:
test_mesh_1d_regional.png
![image](https://user-images.githubusercontent.com/17344536/200441736-972a8136-5c05-4bc9-9bca-b498d972914a.png)


test_mesh_1d_global.png

![image](https://user-images.githubusercontent.com/17344536/200441753-d06e95d1-d85b-4216-9c23-d11ba89a31e4.png)



------
  ## Creating Mesh files for a WRF domain:
For running WRF-CTSM cases, we need to create ESMF mesh files for the WRF domain. We can create mesh file from wrfinput (wrf initial condition files). wrfinput has 2D coordinate information on `XLAT` and `XLONG` variable.

For example, let's create a mesh file from a WRF input file for WRF-CTSM run.  
 ```
./mesh_maker.py --input /glade/scratch/negins/example_files/wrfinput_d01  --output test_mesh_wrf.nc --lat XLAT --lon XLONG --overwrite
```

This produce mesh files for running for our WRF domain. 

Here is how the regional plot looks like for this mesh file:

 ![image](https://user-images.githubusercontent.com/17344536/200442002-1ee5595c-9252-4934-a07c-2f6ad86aff1b.png)


 