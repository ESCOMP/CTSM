import numpy as np

def import_coord_1d(ds, coordName):
    """Import 1-d coordinate variable

    Args:
        ds (xarray Dataset): Dataset whose coordinate you want to import.
        coordName (str): Name of coordinate to import

    Returns:
        xarray DataArray: DataArray corresponding to the requested coordinate.
    """
    da = ds[coordName]
    if len(da.dims) != 1:
        abort(f"Expected 1 dimension for {coordName}; found {len(da.dims)}: {da.dims}")
    return da, len(da)


def import_coord_2d(ds, coordName, varName):
    """Import 2-d latitude or longitude variable from a CESM history file (e.g., name LATIXY or LONGXY) and return it as a 1-d DataArray that can be used as a coordinate for writing CESM input files

    Args:
        ds (xarray Dataset): Dataset whose coordinate you want to import.
        coordName (str): Name of coordinate to import
        varName (str): Name of variable with dimension coordName

    Returns:
        xarray DataArray: 1-d variable that can be used as a coordinate for writing CESM input files
        int: Length of that variable
    """
    da = ds[varName]
    thisDim = [x for x in da.dims if coordName in x]
    if len(thisDim) != 1:
        abort(f"Expected 1 dimension name containing {coordName}; found {len(thisDim)}: {thisDim}")
    thisDim = thisDim[0]
    otherDim = [x for x in da.dims if coordName not in x]
    if len(otherDim) != 1:
        abort(
            f"Expected 1 dimension name not containing {coordName}; found {len(otherDim)}: {otherDim}"
        )
    otherDim = otherDim[0]
    da = da.astype(np.float32)
    da = da.isel({otherDim: [0]}).squeeze().rename({thisDim: coordName}).rename(coordName)
    da = da.assign_coords({coordName: da.values})
    da.attrs["long_name"] = "coordinate " + da.attrs["long_name"]
    da.attrs["units"] = da.attrs["units"].replace(" ", "_")
    return da, len(da)

