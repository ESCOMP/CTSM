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
    Started from a copy of python/ctsm/modify_fsurdat/modify_fsurdat.py
    Deleted unnecessary functions.

    Modified __init__ and added new function set_mesh_mask.

    Other functions remain identical; point to them or repeat code here?
    """

    # TODO If landmask_file ends up being required, lon_1,2 and lat_1,2 will be
    #      removed here and elsewhere
    def __init__(self, my_data, lon_1, lon_2, lat_1, lat_2, landmask_file):

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
    def init_from_file(cls, fsurdat_in, lon_1, lon_2, lat_1, lat_2, landmask_file):
        """Initialize a ModifyFsurdat object from file fsurdat_in"""
        logger.info('Opening fsurdat_in file to be modified: %s', fsurdat_in)
        my_file = xr.open_dataset(fsurdat_in)
        return cls(my_file, lon_1, lon_2, lat_1, lat_2, landmask_file)


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


    def set_mesh_mask(self, var):
        """
        Sets 2d variable var = (1 - rectangle)
        Find the common lat/lon values between var and rectangle and set
        var(centerCoords) = not_rectangle(lat, lon)
        """
        # TODO Research faster way of finding common lat/lon values

        for element_count in self.file['elementCount']:
            for row in self.lsmlat:
                if abs(self.file['centerCoords'][element_count, 1] - \
                       self.lat_2d[row]) < 0.001:
                    for col in self.lsmlon:
                        if abs(self.file['centerCoords'][element_count, 0] - \
                               self.lon_2d[col]) < 0.001:
                            print(element_count)  # Keep for now
                            print(row)  # TODO Remove
                            print(col)  # TODO Remove
                            print(self.file[var][element_count])
                            print(self.not_rectangle[row, col])
                            self.file[var][element_count] = \
                                self.not_rectangle[row, col]
                            break
