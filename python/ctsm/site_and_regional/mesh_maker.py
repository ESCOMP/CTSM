from mesh_type import MeshType
import xarray as xr
import os
oned = True
if oned:
    ifile = "/glade/scratch/negins/this_region_4x5_new_2/surfdata_4x5_hist_78pfts_CMIP6_simyr1850_275.0-330.0_-40-15_c220705.nc"
else:
    ifile = "/glade/scratch/negins/wrf-ctsm_production/WRF/test/em_real_sim1_noahmp/wrfinput_d01"

import os, sys, getopt

if os.path.isfile(ifile):
    ds = xr.open_dataset(ifile, mask_and_scale=False, decode_times=False).transpose()
else:
    print("Input file could not find!")
    sys.exit(2)

    import os, sys, getopt

    if os.path.isfile(ifile):
        ds = xr.open_dataset(
            ifile, mask_and_scale=False, decode_times=False
        ).transpose()
    else:
        print("Input file could not find!")
        sys.exit(2)

if oned:
    lat_name = "lsmlat"
    lon_name = "lsmlon"
else:
    lat_name = "XLAT"
    lon_name = "XLONG"

lats = ds[lat_name]
lons = ds[lon_name]

hasTime = "Time" in ds[lat_name].dims
print(hasTime)

if hasTime:
    lats = lats[:, :, 0]
    lons = lons[:, :, 0]


this_mesh = MeshType(lats, lons)

this_mesh.calculate_corners()
filename = "/glade/scratch/negins/this_region_4x5_new_2/khar2.nc"

this_mesh.create_esmf(filename, area=None)
print("Done")
