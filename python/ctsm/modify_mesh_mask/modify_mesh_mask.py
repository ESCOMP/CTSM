"""
Run this code by using the following wrapper script:
/tools/modify_mesh_mask/mesh_mask_modifier

The wrapper script includes a full description and instructions.
"""

import logging

from math import isclose
import numpy as np
import xarray as xr

from ctsm.config_utils import lon_range_0_to_360

logger = logging.getLogger(__name__)


class ModifyMeshMask:
    """
    Description
    -----------
    """

    # The mesh_mask_modifier tool reads landmask, while the modify_fsurdat tool
    # reads landmask_diff from the landmask file. Sample landmask file:
    # fill_indianocean_slevis.nc located here as of 2022/06/30:
    # /glade/work/slevis/git/mksurfdata_toolchain/tools/modify_fsurdat/fill_indian_ocean/
    # Read landmask_diff here only for consistency checks
    def __init__(self, my_data, landmask_file, lat_dimname, lon_dimname, lat_varname, lon_varname):

        self.file = my_data

        # landmask from user-specified .nc file in the .cfg file
        self._landmask_file = xr.open_dataset(landmask_file)

        assert lat_varname in self._landmask_file.variables
        assert lon_varname in self._landmask_file.variables
        assert lat_dimname in self._landmask_file.dims
        assert lon_dimname in self._landmask_file.dims

        self.latvar = self._landmask_file[lat_varname]
        self.lonvar = self._landmask_file[lon_varname]
        self.lsmlat = self._landmask_file.dims[lat_dimname]
        self.lsmlon = self._landmask_file.dims[lon_dimname]

    @classmethod
    def init_from_file(
        cls, file_in, landmask_file, lat_dimname, lon_dimname, lat_varname, lon_varname
    ):
        """Initialize a ModifyMeshMask object from file_in"""
        logger.info("Opening file to be modified: %s", file_in)
        my_file = xr.open_dataset(file_in)
        return cls(my_file, landmask_file, lat_dimname, lon_dimname, lat_varname, lon_varname)

    def set_mesh_mask(self, var):
        """
        Sets 1d mask variable var = ocnmask = not(landmask).
        Assumes the 1d vector is in the same south-to-north west-to-east
        order as the 2d array.
        """

        # Initialize
        ncount = 0
        landmask_diff = self._landmask_file.landmask_diff
        landmask = self._landmask_file.landmask

        for row in range(self.lsmlat):  # rows from landmask file
            logger.info("row = %d", row + 1)
            for col in range(self.lsmlon):  # cols from landmask file
                errmsg = (
                    "landmask not 0 or 1 at row, col, value = "
                    + f"{row} {col} {landmask[row, col]}"
                )
                assert isclose(landmask[row, col], 0, abs_tol=1e-9) or isclose(
                    landmask[row, col], 1, abs_tol=1e-9
                ), errmsg
                errmsg = (
                    "landmask_diff not 0 or 1 at row, col, value = "
                    + f"{row} {col} {landmask_diff[row, col]}"
                )
                assert isclose(landmask_diff[row, col], 0, abs_tol=1e-9) or isclose(
                    landmask_diff[row, col], 1, abs_tol=1e-9
                ), errmsg
                if int(landmask_diff[row, col]) == 1:
                    errmsg = (
                        "landmask should = landmask_diff where the "
                        + f"latter equals 1, but here landmask = 0 at row, col = {row} {col}"
                    )
                    assert int(landmask[row, col]) == 1, errmsg
                # Reshape landmask into the
                # elementCount dimension of the mesh file.
                # In the process overwrite self.file[var].
                ocnmask = np.logical_not(int(landmask[row, col]))
                self.file[var][ncount] = ocnmask

                # All else in this function supports error checking

                # lon and lat from the landmask file
                if len(self.latvar.sizes) == 2:
                    lat_new = float(self.latvar[row, col])
                    lon_new = float(self.lonvar[row, col])
                else:
                    lat_new = float(self.latvar[row])
                    lon_new = float(self.lonvar[col])
                # ensure lon range of 0-360 rather than -180 to 180
                lon_new = lon_range_0_to_360(lon_new)
                # lon and lat from the mesh file
                lat_mesh = float(self.file["centerCoords"][ncount, 1])
                lon_mesh = float(self.file["centerCoords"][ncount, 0])
                # ensure lon range of 0-360 rather than -180 to 180
                lon_mesh = lon_range_0_to_360(lon_mesh)

                errmsg = (
                    "Must be equal: "
                    " lat_new = "
                    + str(lat_new)
                    + " lat_mesh = "
                    + str(lat_mesh)
                    + " (at ncount = "
                    + str(ncount)
                    + ")"
                )
                assert isclose(lat_new, lat_mesh, abs_tol=1e-5), errmsg
                errmsg = (
                    "Must be equal: "
                    " lon_new = "
                    + str(lon_new)
                    + " lon_mesh = "
                    + str(lon_mesh)
                    + " (at ncount = "
                    + str(ncount)
                    + ")"
                )
                assert isclose(lon_new, lon_mesh, abs_tol=1e-5), errmsg

                # increment counter
                ncount = ncount + 1

        # Error check
        element_count = int(max((self.file["elementCount"])) + 1)
        errmsg = "element_count =" + str(element_count) + "ncount =" + str(ncount) + "must be equal"
        assert ncount == element_count, errmsg
