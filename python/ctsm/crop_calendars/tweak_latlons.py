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
]
file_mesh_in = (
    "/glade/campaign/cesm/cesmdata/inputdata/share/meshes/360x720_120830_ESMFmesh_c20210507_cdf5.nc"
)

file_list_out = []
coord_list = ["lat", "lon"]

# %%


def apply_tweak(ds, coord_str, tweak=1e-6):
    # Apply tweak
    da = ds[coord_str]
    coord2 = da.values
    coord2 += tweak
    while np.any(coord2 == da.values):
        tweak *= 10
        coord2 = da.values
        coord2 += tweak
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
    ds[coord_str] = da2

    # Add a variable with the amount of the tweak
    tweak_attrs = {}
    if "standard_name" in da:
        coord_name = da.attrs["standard_name"]
    else:
        coord_name = da.attrs["long_name"].replace("coordinate_", "")
    tweak_attrs["standard_name"] = coord_name + "_tweak"
    tweak_attrs[
        "long_name"
    ] = f"Amount {coord_name} was shifted to avoid ambiguous nearest neighbors"
    tweak_attrs["units"] = da.attrs["units"]
    da_tweak = xr.DataArray(
        data=coord_tweak,
        coords=new_coords_dict,
        dims=da.dims,
        attrs=tweak_attrs,
    )
    tweak_name = coord_str + "_tweak"
    ds[tweak_name] = da_tweak

    return ds


# %% Apply tweak to all files

for filename in file_list_in:
    file_in = os.path.join(topdir, filename)
    ds = xr.open_dataset(file_in)

    for coord in coord_list:
        ds = apply_tweak(ds, coord, tweak=1e-6)

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

for filename in file_list_out[1:]:
    ds2 = xr.open_dataset(filename)
    for coord in coord_list:
        # Ensure that coordinates are the same
        var = coord
        assert np.array_equal(ds[var].values, ds2[var].values)
        # Ensure that tweaks were the same
        var = coord + "_tweak"
        assert np.array_equal(ds[var].values, ds2[var].values)
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
# %%
