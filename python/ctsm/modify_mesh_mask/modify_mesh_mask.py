"""
Run this code by using the following wrapper script:
/tools/modify_mesh_mask/mesh_mask_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging

from math import isclose
import numpy as np
import xarray as xr

from ctsm.git_utils import get_ctsm_git_short_hash
from ctsm.utils import update_metadata
from ctsm.config_utils import lon_range_0_to_360

logger = logging.getLogger(__name__)

class ModifyMeshMask:
    """
    Description
    -----------
    """
    # TODO landmask here needs all land/ocn, not just the section
    # being changed for modify_fsurdat, so let's include both in the
    # file: landmask_all for modify_meshes and landmask_change for
    # modify_fsurdat.
    def __init__(self, my_data, landmask_file):

        self.file = my_data

        # landmask from user-specified .nc file in the .cfg file
        self._landmask_file = xr.open_dataset(landmask_file)
        self.rectangle = self._landmask_file.landmask  # (lsmlat, lsmlon)
        self.lat_2d = self._landmask_file.lat  # (lsmlat)
        self.lon_2d = self._landmask_file.lon  # (lsmlon)
        self.lsmlat = self._landmask_file.lsmlat
        self.lsmlon = self._landmask_file.lsmlon

        self.not_rectangle = np.logical_not(self.rectangle)


    @classmethod
    def init_from_file(cls, file_in, landmask_file):
        """Initialize a ModifyMeshMask object from file_in"""
        logger.info('Opening file to be modified: %s', file_in)
        my_file = xr.open_dataset(file_in)
        return cls(my_file, landmask_file)


    def write_output(self, file_in, file_out):
        """
        Description
        -----------
        Write output file

        Arguments
        ---------
        file_in:
            (str) Command line entry of input file
        file_out:
            (str) Command line entry of output file
        """

        # update attributes
        # TODO Better as dictionary?
        title = 'Modified mesh file'
        summary = 'Modified mesh file'
        contact = 'N/A'
        data_script = os.path.abspath(__file__) + " -- " + get_ctsm_git_short_hash()
        description = 'Modified this file: ' + file_in
        update_metadata(self.file, title=title, summary=summary,
                        contact=contact, data_script=data_script,
                        description=description)

        # mode 'w' overwrites file if it exists
        self.file.to_netcdf(path=file_out, mode='w',
                            format="NETCDF3_64BIT")
        logger.info('Successfully created: %s', file_out)
        self.file.close()


    def set_mesh_mask(self, var):
        """
        Sets 1d mask variable "var" = 2d mask variable "not_rectangle".
        Assumes that 1d vector is in same south-to-north west-to-east
        order as the 2d array.
        """

        ncount = 0  # initialize
        for row in self.lsmlat:  # rows from landmask file
            logger.info('row = %d', row)
            for col in self.lsmlon:  # cols from landmask file
                # Reshape landmask file's mask (not_rectangle) into the
                # elementCount dimension of the mesh file.
                # In the process overwrite self.file[var].
                self.file[var][ncount] = self.not_rectangle[row, col]

                # All else in this function supports error checking

                # lon and lat from the landmask file
                lat_new = float(self.lat_2d[row])
                lon_new = float(self.lon_2d[col])
                # ensure lon range of 0-360 rather than -180 to 180
                lon_new = lon_range_0_to_360(lon_new)
                # lon and lat from the mesh file
                lat_mesh = float(self.file['centerCoords'][ncount, 1])
                lon_mesh = float(self.file['centerCoords'][ncount, 0])
                # ensure lon range of 0-360 rather than -180 to 180
                lon_mesh = lon_range_0_to_360(lon_mesh)

                errmsg = 'Must be equal: ' \
                         ' lat_new = ' + str(lat_new) + \
                         ' lat_mesh = ' + str(lat_mesh) + \
                         ' (at ncount = ' + str(ncount) + ')'
                assert isclose(lat_new, lat_mesh, abs_tol=1e-5), errmsg
                errmsg = 'Must be equal: ' \
                         ' lon_new = ' + str(lon_new) + \
                         ' lon_mesh = ' + str(lon_mesh) + \
                         ' (at ncount = ' + str(ncount) + ')'
                assert isclose(lon_new, lon_mesh, abs_tol=1e-5), errmsg

                # increment counter
                ncount = ncount + 1

        # Error check
        element_count = int(max((self.file['elementCount'])) + 1)
        errmsg = 'element_count =' + str(element_count) + \
                 'ncount =' + str(ncount) + 'must be equal'
        assert ncount == element_count, errmsg
