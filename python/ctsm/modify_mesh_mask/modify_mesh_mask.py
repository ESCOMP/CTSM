"""
Run this code by using the following wrapper script:
/tools/modify_mesh_mask/mesh_mask_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging

import numpy as np
import xarray as xr

from ctsm.git_utils import get_ctsm_git_short_hash
from ctsm.utils import abort, update_metadata
from ctsm.config_utils import lon_range_0_to_360

logger = logging.getLogger(__name__)

class ModifyMeshMask:
    """
    Description
    -----------
    Notes
    - Started from a copy of python/ctsm/modify_fsurdat/modify_fsurdat.py
    - Functions here haven't changed, so...
      TODO Change fsurdat to mesh_mask as needed OR make these functions
           generic by moving to utils?
    """

    def __init__(self, my_data, lon_1, lon_2, lat_1, lat_2, landmask_file):

        self.file = my_data

        self.rectangle = self._get_rectangle(
            lon_1=lon_1, lon_2=lon_2,
            lat_1=lat_1, lat_2=lat_2,
            longxy=self.file.LONGXY, latixy=self.file.LATIXY)

        if landmask_file is not None:
            # overwrite self.not_rectangle with data from
            # user-specified .nc file in the .cfg file
            self._landmask_file = xr.open_dataset(landmask_file)
            self.rectangle = self._landmask_file.landmask

        self.not_rectangle = np.logical_not(self.rectangle)


    @classmethod
    def init_from_file(cls, fsurdat_in, lon_1, lon_2, lat_1, lat_2, landmask_file):
        """Initialize a ModifyFsurdat object from file fsurdat_in"""
        logger.info('Opening fsurdat_in file to be modified: %s', fsurdat_in)
        my_file = xr.open_dataset(fsurdat_in)
        return cls(my_file, lon_1, lon_2, lat_1, lat_2, landmask_file)


    @staticmethod
    def _get_rectangle(lon_1, lon_2, lat_1, lat_2, longxy, latixy):
        """
        Description
        -----------
        """

        # ensure that lon ranges 0-360 in case user entered -180 to 180
        lon_1 = lon_range_0_to_360(lon_1)
        lon_2 = lon_range_0_to_360(lon_2)

        # determine the rectangle(s)
        # TODO This is not really "nearest" for the edges but isel didn't work
        rectangle_1 = (longxy >= lon_1)
        rectangle_2 = (longxy <= lon_2)
        eps = np.finfo(np.float32).eps  # to avoid roundoff issue
        rectangle_3 = (latixy >= (lat_1 - eps))
        rectangle_4 = (latixy <= (lat_2 + eps))

        if lon_1 <= lon_2:
            # rectangles overlap
            union_1 = np.logical_and(rectangle_1, rectangle_2)
        else:
            # rectangles don't overlap: stradling the 0-degree meridian
            union_1 = np.logical_or(rectangle_1, rectangle_2)

        if lat_1 < -90 or lat_1 > 90 or lat_2 < -90 or lat_2 > 90:
            errmsg = 'lat_1 and lat_2 need to be in the range -90 to 90'
            abort(errmsg)
        elif lat_1 <= lat_2:
            # rectangles overlap
            union_2 = np.logical_and(rectangle_3, rectangle_4)
        else:
            # rectangles don't overlap: one in the north, one in the south
            union_2 = np.logical_or(rectangle_3, rectangle_4)

        # union rectangles overlap
        rectangle = np.logical_and(union_1, union_2)

        return rectangle


    def write_output(self, fsurdat_in, fsurdat_out):
        """
        Description
        -----------
        Write output file

        Arguments
        ---------
        fsurdat_in:
            (str) Command line entry of input surface dataset
        fsurdat_out:
            (str) Command line entry of output surface dataset
        """

        # update attributes
        # TODO Better as dictionary?
        title = 'Modified fsurdat file'
        summary = 'Modified fsurdat file'
        contact = 'N/A'
        data_script = os.path.abspath(__file__) + " -- " + get_ctsm_git_short_hash()
        description = 'Modified this file: ' + fsurdat_in
        update_metadata(self.file, title=title, summary=summary,
                        contact=contact, data_script=data_script,
                        description=description)

        # abort if output file already exists
        file_exists = os.path.exists(fsurdat_out)
        if file_exists:
            errmsg = 'Output file already exists: ' + fsurdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        self.file.to_netcdf(path=fsurdat_out, mode='w',
                            format="NETCDF3_64BIT")
        logger.info('Successfully created fsurdat_out: %s', fsurdat_out)
        self.file.close()


    def setvar_lev0(self, var, val):
        """
        Sets 2d variable var to value val in user-defined rectangle,
        defined as "other" in the function
        """
        self.file[var] = self.file[var].where(
            self.not_rectangle, other=val)
