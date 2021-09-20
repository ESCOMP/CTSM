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

    def __init__(self, fsurdat_in):

        print("Open file: " + fsurdat_in)
        self.file = xr.open_dataset(fsurdat_in)


    def write_output(self, fsurdat_in, fsurdat_out):
        """
        Description
        -----------
        Write output file
        """

        # specify dimension order
        self.file = self.file.transpose(u'time', u'cft', u'lsmpft', u'natpft',
                                        u'nglcec', u'nglcecp1', u'nlevsoi',
                                        u'nlevurb', u'numrad', u'numurbl',
                                        'lsmlat', 'lsmlon')

        # update attributes
        title = 'Modified fsurdat file'
        summary = 'Modified fsurdat file'
        contact = 'N/A'
        data_script = os.path.abspath(__file__) + " -- " + get_git_sha()
        description = 'Modified this file: ' + fsurdat_in
        update_metadata(self.file, title, summary, contact, data_script,
                        description)

        # abort if output file already exists
        file_exists = os.path.exists(fsurdat_out)
        if file_exists:
            errmsg = 'Output file already exists: ' + fsurdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        self.file.to_netcdf(path=fsurdat_out, mode='w')
        print('Successfully created fsurdat_out:' + fsurdat_out)
        self.file.close()


    def overwrite_single_pft(self, dom_pft):
        """
        Description
        -----------
        """

        self.file['PCT_NAT_PFT'][:,:,:] = 0
        self.file['PCT_NAT_PFT'][:,:,dom_pft] = 100


    def zero_nonveg_lu(self):
        """
        Description
        -----------
        """

        self.file['PCT_NATVEG'][:,:] = 100
        self.file['PCT_CROP'][:,:] = 0
        self.file['PCT_LAKE'][:,:] = 0
        self.file['PCT_WETLAND'][:,:] = 0
        self.file['PCT_URBAN'][:,:,] = 0
        self.file['PCT_GLACIER'][:,:] = 0


    def uni_snowpack(self):
        """
        Description
        -----------
        """

        self.file['STD_ELEV'][:,:] = 20


    def no_saturation_excess(self):
        """
        Description
        -----------
        """

        self.file['FMAX'][:,:] = 0
