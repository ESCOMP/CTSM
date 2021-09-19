"""
Run this code by using the following wrapper script:
/tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import os

import xarray as xr

from ctsm.utils import abort, get_git_sha, update_metadata

class ModifyFsurdat:
    """
    A case to encapsulate function that modifies fsurdat file.

    ...

    Attributes
    ----------

    Methods
    -------
    modify:
        Modify input surface dataset
    """

    def __init__(self, fsurdat_in, fsurdat_out, overwrite_single_pft,
                 dom_pft, zero_nonveg_lu, uni_snowpack,
                 no_saturation_excess):
        self._fsurdat_in = fsurdat_in
        self._fsurdat_out = fsurdat_out
        self._overwrite_single_pft = overwrite_single_pft
        self._dom_pft = dom_pft
        self._zero_nonveg_lu = zero_nonveg_lu
        self._uni_snowpack = uni_snowpack
        self._no_saturation_excess = no_saturation_excess


    def modify(self):
        """
        Description
        -----------
        Modify variable of output file
        """

        print ("Creating surface dataset file")
        file_in = self._fsurdat_in
        print("Open file: " + file_in)
        file = xr.open_dataset(file_in)

        # modify surface data properties
        if self._overwrite_single_pft:
            file['PCT_NAT_PFT'][:,:,:] = 0
            file['PCT_NAT_PFT'][:,:,self._dom_pft] = 100
        if self._zero_nonveg_lu:
            file['PCT_NATVEG'][:,:]  = 100
            file['PCT_CROP'][:,:]    = 0
            file['PCT_LAKE'][:,:]    = 0.
            file['PCT_WETLAND'][:,:] = 0.
            file['PCT_URBAN'][:,:,]   = 0.
            file['PCT_GLACIER'][:,:] = 0.
        if self._uni_snowpack:
            file['STD_ELEV'][:,:] = 20.
        if self._no_saturation_excess:
            file['FMAX'][:,:] = 0.

        # specify dimension order
        file = file.transpose(u'time', u'cft', u'lsmpft', u'natpft', u'nglcec',
                              u'nglcecp1', u'nlevsoi', u'nlevurb', u'numrad',
                              u'numurbl', 'lsmlat', 'lsmlon')

        # update attributes
        title = 'Modified fsurdat file'
        summary = 'Modified fsurdat file'
        contact = 'N/A'
        data_script = os.path.abspath(__file__) + " -- " + get_git_sha()
        description = 'Modified this file: ' + file_in
        update_metadata(file, title, summary, contact, data_script, description)

        # abort if output file already exists
        file_exists = os.path.exists(self._fsurdat_out)
        if file_exists:
            errmsg = 'Output file already exists: ' + self._fsurdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        file.to_netcdf(path=self._fsurdat_out, mode='w')
        print('Successfully created file (fsurdat_out) :' + self._fsurdat_out)
        file.close()
