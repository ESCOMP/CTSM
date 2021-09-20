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
    write_output:
        Write out the modified fsurdat (surface dataset)

    TODO:
        Complete the list of methods
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
        _title = 'Modified fsurdat file'
        _summary = 'Modified fsurdat file'
        _contact = 'N/A'
        _data_script = os.path.abspath(__file__) + " -- " + get_git_sha()
        _description = 'Modified this file: ' + fsurdat_in
        update_metadata(self.file, title=_title, summary=_summary,
                        contact=_contact, data_script=_data_script,
                        description=_description)

        # abort if output file already exists
        file_exists = os.path.exists(fsurdat_out)
        if file_exists:
            errmsg = 'Output file already exists: ' + fsurdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        self.file.to_netcdf(path=fsurdat_out, mode='w')
        print('Successfully created fsurdat_out: ' + fsurdat_out)
        self.file.close()


    def dom_pft(self, dom_pft):
        """
        Description
        -----------
        Replace fsurdat file's PCT_NAT_PFT values with:
        - 100 for dom_pft selected by user
        - 0 for all other PFTs
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


    def max_sat_area(self, max_sat_area):
        """
        Description
        -----------
        Replace fsurdat file's FMAX values with a constant selected by user
        """

        self.file['FMAX'][:,:] = max_sat_area
