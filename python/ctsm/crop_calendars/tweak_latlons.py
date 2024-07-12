# %%
import numpy as np
import xarray as xr
import os
import sys
from netCDF4 import Dataset
import contextlib

# -- add python/ctsm  to path (needed if we want to run this stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.mesh_maker import main as mesh_maker

topdir = "/glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/cropdata/calendars/processed/"
file_list_in = [
    "swindow_starts_ggcmi_crop_calendar_phase3_v1.01.2000-2000.20231005_145103.nc",
    "swindow_ends_ggcmi_crop_calendar_phase3_v1.01.2000-2000.20231005_145103.nc",
    "gdds_20230829_161011.nc",
    "gdd20bl.copied_from.gdds_20230829_161011.v2.nc",
    "sdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-hcru_hcru_mt13.2000-2000.20230728_165845.nc",
    "hdates_ggcmi_crop_calendar_phase3_v1.01_nninterp-hcru_hcru_mt13.2000-2000.20230728_165845.nc",
    "/glade/work/samrabin/cropCals_testing_20240626/gdds_20240712_114642_10x15_interpd_halfdeg.nc",
    "/glade/work/samrabin/gdd20_baselines/gswp3.10x15_interpd_halfdeg.1980-2009.nc",
]
file_mesh_in = (
    "/glade/campaign/cesm/cesmdata/inputdata/share/meshes/360x720_120830_ESMFmesh_c20210507_cdf5.nc"
)

file_list_out = []
coord_list = ["lat", "lon"]
COORD_DATATYPE = np.float64

# %% Define functions

def get_ds(topdir, file_in):
    if not os.path.exists(file_in):
        file_in = os.path.join(topdir, file_in)
    ds = xr.open_dataset(file_in)
    return file_in, ds

def get_tweak(ds_in, coord_str, init_tweak):
    """
    Get the tweak that will be applied to all datasets' lat/lon coordinates
    """
    da = ds_in[coord_str]
    coord2_orig = da.values.astype(COORD_DATATYPE)
    coord2 = coord2_orig
    tweak = init_tweak
    coord2 += tweak

    # This is necessary if precision is lower than float64
    max_tweak = 1e-2
    while np.any(coord2 == da.values):
        tweak *= 10
        if tweak > max_tweak:
            raise RuntimeError(f"Tweaking by +{max_tweak} failed to 'take'")
        coord2 = coord2_orig
        coord2 += tweak
    return tweak

def apply_tweak(ds_in, coord_str, tweak):
    # Apply tweak
    da = ds_in[coord_str]
    coord2 = da.values.astype(COORD_DATATYPE)
    coord2 += tweak
    if np.any(coord2 == da.values):
        raise RuntimeError('Tweak didn''t "take"')
    coord_tweak = np.full_like(coord2, tweak)

    # Ensure that no value is above maximum in input data. This is needed for mesh_maker to work.
    max_coord = np.max(da.values)
    where_toohigh = np.where(coord2 > max_coord)
    Ntoohigh = len(where_toohigh[0])
    if Ntoohigh != 1:
        raise RuntimeError(f"Expected 1 coordinate value too high; got {Ntoohigh}")
    coord2[where_toohigh] = max_coord
    coord_tweak[where_toohigh] = max_coord

    # Convert to DataArray
    new_coords_dict = {coord_str: coord2}
    da2 = xr.DataArray(
        data=coord2,
        coords=new_coords_dict,
        dims=da.dims,
        attrs=da.attrs,
    )

    # Replace coordinate in dataset
    ds_in[coord_str] = da2

    # Add a variable with the amount of the tweak
    tweak_attrs = {}
    if "standard_name" in da.attrs:
        coord_name = da.attrs["standard_name"]
    elif "long_name" in da.attrs:
        coord_name = da.attrs["long_name"].replace("coordinate_", "")
    else:
        coord_name = coord_str
    tweak_attrs["standard_name"] = coord_name + "_tweak"
    tweak_attrs[
        "long_name"
    ] = f"Amount {coord_name} was shifted to avoid ambiguous nearest neighbors"
    if "units" in da.attrs:
        tweak_attrs["units"] = da.attrs["units"]
    da_tweak = xr.DataArray(
        data=coord_tweak,
        coords=new_coords_dict,
        dims=da.dims,
        attrs=tweak_attrs,
    )
    tweak_name = coord_str + "_tweak"
    ds_in[tweak_name] = da_tweak

    return ds_in

def check(ds, f0_base, ds2, f_base, var):
    if not np.array_equal(ds[var].values, ds2[var].values):
        if not np.array_equal(ds[var].shape, ds2[var].shape):
            msg = f"{var} shapes differ b/w {f0_base} ({ds[var].shape}) and {f_base} ({ds2[var].shape})"
            raise RuntimeError(msg)
        max_diff = np.max(np.abs(ds[var].values - ds2[var].values))
        msg = f"{var}s differ between {f0_base} and {f_base}; max = {max_diff}"
        type0 = type(ds[var].values[0])
        type2 = type(ds2[var].values[0])
        if type0 != type2:
            msg += f"\nTypes also differ: {type0} vs. {type2}"
        raise RuntimeError(msg)

# %% Apply tweak to all files

# Set up empty dicts
tweak_dict = {}
coord_type_dict = {}
for coord in coord_list:
    tweak_dict[coord] = -np.inf

# Get tweaks
for file_in in file_list_in:
    file_in, ds = get_ds(topdir, file_in)
    for coord in coord_list:
        this_tweak = get_tweak(ds, coord, init_tweak=1e-6)
        if this_tweak > tweak_dict[coord]:
            tweak_dict[coord] = this_tweak
for coord in coord_list:
    print(f"Tweaking {coord} by {tweak_dict[coord]}")
print(" ")

# Apply tweaks
for file_in in file_list_in:
    file_in, ds = get_ds(topdir, file_in)

    for coord in coord_list:
        ds = apply_tweak(ds, coord, tweak_dict[coord])

    # Set up for save
    file_out = file_in.replace(".nc", ".tweaked_latlons.nc")
    with Dataset(file_in, "r") as netcdf_file:
        netcdf_format = netcdf_file.data_model

    # Save
    print(f"Saving {file_out}")
    ds.to_netcdf(file_out, format=netcdf_format)
    file_list_out.append(file_out)
print("Done")


# %% Ensure all files got the same tweaks

ds = xr.open_dataset(file_list_out[0])
f0_base = os.path.basename(file_list_out[0])

for filename in file_list_out[1:]:
    ds2 = xr.open_dataset(filename)
    f_base = os.path.basename(filename)
    for coord in coord_list:
        check(ds, f0_base, ds2, f_base, coord)
        check(ds, f0_base, ds2, f_base, coord + "_tweak")
print("All good!")


# %% Save new mesh file

outfile_name = os.path.basename(file_mesh_in)
outfile_name = outfile_name.replace(".nc", ".tweaked_latlons.nc")
outdir = os.path.dirname(file_list_out[0])
file_mesh_out = os.path.join(outdir, outfile_name)

@contextlib.contextmanager
def redirect_argv(arglist):
    argv_tmp = sys.argv[:]
    sys.argv=arglist
    yield
    sys.argv = argv_tmp


mesh_maker_args = [
    "mesh_maker",
    "--input",
    file_list_out[0],
    "--output",
    file_mesh_out,
    "--lat",
    "lat",
    "--lon",
    "lon",
    "--overwrite",
]
print(f"Saving {file_mesh_out}...")
with redirect_argv(mesh_maker_args):
    mesh_maker()

# Change format, if needed
with Dataset(file_mesh_in, "r") as netcdf_file:
    netcdf_format_in = netcdf_file.data_model
with Dataset(file_mesh_out, "r") as netcdf_file:
    netcdf_format_out = netcdf_file.data_model
if netcdf_format_in != netcdf_format_out:
    file_mesh_out_tmp = file_mesh_out + ".tmp"
    os.rename(file_mesh_out, file_mesh_out_tmp)
    ds = xr.open_dataset(file_mesh_out_tmp)
    ds.to_netcdf(file_mesh_out, format=netcdf_format_in)
    os.remove(file_mesh_out_tmp)


print("Done")
