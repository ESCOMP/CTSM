"""
Run this code by using the following wrapper script:
/tools/modify_input_files/mesh_mask_modifier

The wrapper script includes a full description and instructions.
"""

import logging

from math import isclose
import numpy as np
import xarray as xr

from ctsm.utils import abort
from ctsm.longitude import Longitude

logger = logging.getLogger(__name__)


class ModifyMeshMask:
    """
    Description
    -----------
    """

    # The mesh_mask_modifier tool reads landmask, while the modify_fsurdat tool
    # reads mod_lnd_props from the landmask file. Sample landmask
    # file: fill_indianocean_slevis.nc located here as of 2022/8/1:
    # /glade/work/slevis/git/mksurfdata_toolchain/tools/modify_input_files ...
    # ... /islas_examples/modify_fsurdat/fill_indian_ocean/
    # Read mod_lnd_props here only for consistency checks
    def __init__(
        self,
        my_data,
        *,
        landmask_file,
        lat_dimname,
        lon_dimname,
        lat_varname,
        lon_varname,
        lon_type,
    ):

        self.file = my_data

        # landmask from user-specified .nc file in the .cfg file
        self._landmask_file = xr.open_dataset(landmask_file)

        assert lat_varname in self._landmask_file.variables
        assert lon_varname in self._landmask_file.variables
        assert lat_dimname in self._landmask_file.dims
        assert lon_dimname in self._landmask_file.dims

        # latvar, lonvar same as right-hand-sides, which may be 1D or 2D
        self.latvar = self._landmask_file[lat_varname][:, ...]
        self.lonvar = self._landmask_file[lon_varname][..., :]
        self.lsmlat = self._landmask_file.dims[lat_dimname]
        self.lsmlon = self._landmask_file.dims[lon_dimname]
        self.lon_type = lon_type

        lonvar_first = self.lonvar[..., 0].data.max()
        lonvar_last = self.lonvar[..., -1].data.max()

        # If self.file["centerCoords"][0, 0] < 0, this suggests a
        # -180 to 180 grid, so make self.lonvar -180 to 180 if not already.
        # If self.file["centerCoords"][0, 0] >= 0, this suggests a
        # 0 to 360 grid, so make self.lonvar 0 to 360 if not already.
        if (self.file["centerCoords"][0, 0] < 0 <= lonvar_first) or (
            self.file["centerCoords"][0, 0] >= 0 > lonvar_first
        ):
            logger.info(
                "first lon_mesh = %s, last lon_mesh = %s",
                self.file["centerCoords"][0, 0].data,
                self.file["centerCoords"][-1, 0].data,
            )
            logger.info("first lonvar = %s, last lonvar = %s.", lonvar_first, lonvar_last)
            explanation = (
                "For consistency in their order, changing lonvar to "
                + "lon_mesh's convention and later (in set_mesh_mask) "
                + "changing any negative longitude values to their "
                + "corresponding positive values."
            )
            logger.info(explanation)
            self._landmask_file = self._landmask_file.roll(lsmlon=self.lsmlon // 2)
        self.lonvar = self._landmask_file[lon_varname][..., :]  # update lonvar

    @classmethod
    def init_from_file(
        cls, *, file_in, landmask_file, lat_dimname, lon_dimname, lat_varname, lon_varname, lon_type
    ):
        """Initialize a ModifyMeshMask object from file_in"""
        logger.info("Opening file to be modified: %s", file_in)
        my_file = xr.open_dataset(file_in)
        return cls(
            my_file,
            landmask_file=landmask_file,
            lat_dimname=lat_dimname,
            lon_dimname=lon_dimname,
            lat_varname=lat_varname,
            lon_varname=lon_varname,
            lon_type=lon_type,
        )

    def set_mesh_mask(self, var):
        """
        Sets 1d mask variable var = ocnmask = not(landmask).
        Assumes the 1d vector is in the same south-to-north west-to-east
        order as the 2d array.
        """

        # Initialize
        ncount = 0
        mod_lnd_props = self._landmask_file.mod_lnd_props
        landmask = self._landmask_file.landmask

        for row in range(self.lsmlat):  # rows from landmask file
            logger.info("row = %d", row + 1)
            for col in range(self.lsmlon):  # cols from landmask file
                errmsg = (
                    "landmask not 0 or 1 at row, col, value = "
                    + f"{row} {col} {landmask[row, col]}"
                )
                assert landmask[row, col] == 0 or landmask[row, col] == 1, errmsg
                errmsg = (
                    "mod_lnd_props not 0 or 1 at row, col, value = "
                    + f"{row} {col} {mod_lnd_props[row, col]}"
                )
                assert mod_lnd_props[row, col] == 0 or mod_lnd_props[row, col] == 1, errmsg
                if int(mod_lnd_props[row, col]) == 1:
                    errmsg = (
                        "landmask should = mod_lnd_props where the "
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
                latvar_scalar = None
                if len(self.latvar.sizes) == 2:
                    latvar_scalar = float(self.latvar[row, col])
                    lonvar_scalar = float(self.lonvar[row, col])
                elif len(self.latvar.sizes) == 1:
                    latvar_scalar = float(self.latvar[row])
                    lonvar_scalar = float(self.lonvar[col])
                else:
                    errmsg = (
                        "ERROR: Expecting latvar.sizes == 1 or 2, not "
                        + f"{len(self.latvar.sizes)}"
                    )
                    abort(errmsg)

                # lon and lat from the mesh file
                lat_mesh = float(self.file["centerCoords"][ncount, 1])
                lon_mesh = float(self.file["centerCoords"][ncount, 0])
                # ensure lon range of 0-360 rather than -180 to 180
                lonvar_scalar = Longitude(lonvar_scalar, self.lon_type).get(360)

                errmsg = (
                    "Must be equal: "
                    " latvar_scalar = "
                    + str(latvar_scalar)
                    + " lat_mesh = "
                    + str(lat_mesh)
                    + " (at ncount = "
                    + str(ncount)
                    + ")"
                )
                assert isclose(latvar_scalar, lat_mesh, abs_tol=1e-5), errmsg
                errmsg = (
                    "Must be equal: "
                    " lonvar_scalar = "
                    + str(lonvar_scalar)
                    + " lon_mesh = "
                    + str(lon_mesh)
                    + " (at ncount = "
                    + str(ncount)
                    + ")"
                )
                assert isclose(lonvar_scalar, lon_mesh, abs_tol=1e-5), errmsg

                # increment counter
                ncount = ncount + 1

        # Error check
        element_count = int(max((self.file["elementCount"])) + 1)
        errmsg = "element_count =" + str(element_count) + "ncount =" + str(ncount) + "must be equal"
        assert ncount == element_count, errmsg
