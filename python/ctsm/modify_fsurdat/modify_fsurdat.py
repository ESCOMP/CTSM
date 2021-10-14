"""
Run this code by using the following wrapper script:
/tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import os

import numpy as np
import xarray as xr

from ctsm.utils import abort, get_git_sha, update_metadata, lon_range_0_to_360

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

    TODO: List of methods here?
    """

    def __init__(self, fsurdat_in):

        print("Open file: " + fsurdat_in)
        self.file = xr.open_dataset(fsurdat_in)


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
        self.file.to_netcdf(path=fsurdat_out, mode='w',
                            format="NETCDF3_64BIT")
        print('Successfully created fsurdat_out: ' + fsurdat_out)
        self.file.close()


    def dom_nat_pft(self, dom_nat_pft):
        """
        Description
        -----------
        Replace fsurdat file's PCT_NAT_PFT values with:
        - 100 for dom_nat_pft selected by user
        - 0 for all other non-crop PFTs

        Arguments
        ---------
        dom_nat_pft:
            (int) Command line entry of PFT to be set to 100% everywhere
        """

        self.file['PCT_NAT_PFT'][:,:,:] = 0
        self.file['PCT_NAT_PFT'][dom_nat_pft,:,:] = 100


    def zero_nonveg(self):
        """
        Description
        -----------
        Set all landunit weights to 0 except the natural vegetation landunit.
        Set that one to 100%.

        Arguments
        ---------
        """

        self.file['PCT_NATVEG'][:,:] = 100
        self.file['PCT_CROP'][:,:] = 0
        self.file['PCT_LAKE'][:,:] = 0
        self.file['PCT_WETLAND'][:,:] = 0
        self.file['PCT_URBAN'][:,:] = 0
        self.file['PCT_GLACIER'][:,:] = 0


    def std_elev(self, std_elev):
        """
        Description
        -----------
        Set STD_ELEV to specify uniform snowpack.

        Arguments
        ---------
        std_elev:
            (float) command line entry of STD_ELEV value
        """

        self.file['STD_ELEV'][:,:] = std_elev


    def max_sat_area(self, max_sat_area):
        """
        Description
        -----------
        Replace fsurdat file's FMAX values with a constant selected by user

        Arguments
        ---------
        max_sat_area:
            (float) command line entry of FMAX value
        """

        min_max_sat_area = 0
        max_max_sat_area = 1
        if min_max_sat_area <= max_sat_area <= max_max_sat_area:
            self.file['FMAX'][:,:] = max_sat_area
        else:
            errmsg = 'Argument --max_sat_area requires values from ' + \
                     str(min_max_sat_area) + ' to ' + str(max_max_sat_area)
            abort(errmsg)


    def set_in_rectangle(self, idealized, lon_in_1, lon_in_2, lat_in_1,
                        lat_in_2, dom_nat_pft, lai, sai, hgt_top, hgt_bot,
                        zero_nonveg, std_elev, max_sat_area):
        """
        Description
        -----------
        Set fsurdat variables in a rectangle defined by lon/lat limits

        Arguments
        ---------
        lon_in_1:
            (int) westernmost edge of rectangle
        lon_in_2:
            (int) easternmost edge of rectangle
        lat_in_1:
            (int) southernmost edge of rectangle
        lat_in_2:
            (int) northernmost edge of rectangle
        dom_nat_pft:
            (int) user-selected PFT to be set to 100% in rectangle
        lai:
            (list of floats) user-defined leaf area index
        sai:
            (list of floats) user-defined stem area index
        hgt_top:
            (list of floats) user-defined height at top of vegetation canopy (m)
        hgt_bot:
            (list of floats) user-defined height at bottom of vegetation canopy (m)
        """

        # default to bare soil
        if dom_nat_pft == -999:
            dom_nat_pft = 0
        if dom_nat_pft == 0:
            lai = [0,0,0,0,0,0,0,0,0,0,0,0]
            sai = [0,0,0,0,0,0,0,0,0,0,0,0]
            hgt_top = [0,0,0,0,0,0,0,0,0,0,0,0]
            hgt_bot = [0,0,0,0,0,0,0,0,0,0,0,0]

        # Currently type(lon_in_*) = int with required range 0-360.
        # If instead of requiring integer values, we decide to allow floats,
        # then ensure that lon ranges 0-360 in case user entered -180 to 180.
        lon_1 = lon_range_0_to_360(lon_in_1)
        lon_2 = lon_range_0_to_360(lon_in_2)

        # determine the rectangle(s)
        # TODO This is not really "nearest" for the edges but isel didn't work
        rectangle_1 = (self.file.LONGXY >= lon_1)
        rectangle_2 = (self.file.LONGXY <= lon_2)
        eps = np.finfo(np.float32).eps  # to avoid roundoff issue
        rectangle_3 = (self.file.LATIXY >= (lat_in_1 - eps))
        rectangle_4 = (self.file.LATIXY <= (lat_in_2 + eps))

        if lon_1 <= lon_2:
            # rectangles overlap
            union_1 = np.logical_and(rectangle_1, rectangle_2)
        else:
            # rectangles don't overlap (stradling the 0-degree meridian)
            union_1 = np.logical_or(rectangle_1, rectangle_2)

        if lat_in_1 <= lat_in_2:
            # rectangles overlap
            union_2 = np.logical_and(rectangle_3, rectangle_4)
        else:
            # rectangles don't overlap (one in the north, on in the south)
            union_2 = np.logical_or(rectangle_3, rectangle_4)

        # rectangles overlap
        rectangle = np.logical_and(union_1, union_2)
        not_rectangle = np.logical_not(rectangle)
        print(rectangle)

        # Overwrite in rectangle(s)
        # ------------------------
        # If idealized, the user makes changes to variables as follows.
        # "other" assigns the corresponding value in the rectangle.
        # Values outside the rectangle are preserved.
        # ------------------------
        # 2D variables
        # max inundated fraction
        self.file['F0'] = \
         self.file['F0'].where(not_rectangle, other=0)
        # max saturated area
        self.file['FMAX'] = \
         self.file['FMAX'].where(not_rectangle, other=max_sat_area)
        # standard deviation of elevation
        self.file['STD_ELEV'] = \
         self.file['STD_ELEV'].where(not_rectangle, other=std_elev)
        # mean topographic slope
        self.file['SLOPE'] = \
         self.file['SLOPE'].where(not_rectangle, other=0)
        # is 10 seem reasonable?
        self.file['zbedrock'] = \
         self.file['zbedrock'].where(not_rectangle, other=10)
        # value representing loam
        self.file['SOIL_COLOR'] = \
         self.file['SOIL_COLOR'].where(not_rectangle, other=15)

        self.file['PFTDATA_MASK'] = \
         self.file['PFTDATA_MASK'].where(not_rectangle, other=1)
        self.file['LANDFRAC_PFT'] = \
         self.file['LANDFRAC_PFT'].where(not_rectangle, other=1)
        self.file['PCT_WETLAND'] = \
         self.file['PCT_WETLAND'].where(not_rectangle, other=0)
        self.file['PCT_CROP'] = \
         self.file['PCT_CROP'].where(not_rectangle, other=0)
        self.file['PCT_LAKE'] = \
         self.file['PCT_LAKE'].where(not_rectangle, other=0)
        self.file['PCT_URBAN'] = \
         self.file['PCT_URBAN'].where(not_rectangle, other=0)
        self.file['PCT_GLACIER'] = \
         self.file['PCT_GLACIER'].where(not_rectangle, other=0)
        self.file['PCT_NATVEG'] = \
         self.file['PCT_NATVEG'].where(not_rectangle, other=100)

        for lev in self.file.nlevsoi:
            # set next three 3D variables to values representing loam
            self.file['PCT_SAND'][lev,:,:] = \
             self.file['PCT_SAND'][lev,:,:].where(not_rectangle, other=43)
            self.file['PCT_CLAY'][lev,:,:] = \
             self.file['PCT_CLAY'][lev,:,:].where(not_rectangle, other=18)
            self.file['ORGANIC'][lev,:,:] = \
             self.file['ORGANIC'][lev,:,:].where(not_rectangle, other=0)

        for pft in self.file.natpft:
            # initialize 3D variable
            self.file['PCT_NAT_PFT'][pft,:,:] = \
             self.file['PCT_NAT_PFT'][pft,:,:].where(not_rectangle, other=0)
            for mon in (self.file.time - 1):
                # initialize 4D variables
                self.file['MONTHLY_LAI'][mon,pft,:,:] = \
                 self.file['MONTHLY_LAI'][mon,pft,:,:].where(not_rectangle, other=0)
                self.file['MONTHLY_SAI'][mon,pft,:,:] = \
                 self.file['MONTHLY_SAI'][mon,pft,:,:].where(not_rectangle, other=0)
                self.file['MONTHLY_HEIGHT_TOP'][mon,pft,:,:] = \
                 self.file['MONTHLY_HEIGHT_TOP'][mon,pft,:,:].where(not_rectangle, other=0)
                self.file['MONTHLY_HEIGHT_BOT'][mon,pft,:,:] = \
                 self.file['MONTHLY_HEIGHT_BOT'][mon,pft,:,:].where(not_rectangle, other=0)
                if pft == dom_nat_pft:
                    # set 4D variables to values for dom_nat_pft
                    self.file['MONTHLY_LAI'][mon,pft,:,:] = \
                     self.file['MONTHLY_LAI'][mon,pft,:,:].where(not_rectangle, other=lai[int(mon)])
                    self.file['MONTHLY_SAI'][mon,pft,:,:] = \
                     self.file['MONTHLY_SAI'][mon,pft,:,:].where(not_rectangle, other=sai[int(mon)])
                    self.file['MONTHLY_HEIGHT_TOP'][mon,pft,:,:] = \
                     self.file['MONTHLY_HEIGHT_TOP'][mon,pft,:,:].where(not_rectangle, other=hgt_top[int(mon)])
                    self.file['MONTHLY_HEIGHT_BOT'][mon,pft,:,:] = \
                     self.file['MONTHLY_HEIGHT_BOT'][mon,pft,:,:].where(not_rectangle, other=hgt_bot[int(mon)])

        for crop in self.file.cft:
            cft_local = crop - (pft + 1)
            # initialize 3D variable
            self.file['PCT_CFT'][cft_local,:,:] = \
             self.file['PCT_CFT'][cft_local,:,:].where(not_rectangle, other=0)
        # set 3D variables
        self.file['PCT_NAT_PFT'][dom_nat_pft,:,:] = \
         self.file['PCT_NAT_PFT'][dom_nat_pft,:,:].where(not_rectangle, other=100)
        # required even though PCT_CROP = 0; the sum of all crops must = 100
        self.file['PCT_CFT'][0,:,:] = \
         self.file['PCT_CFT'][0,:,:].where(not_rectangle, other=100)



#       -------------------------------------------------------------------

#       # If not idealized, the user wants changes in specific vars in the
#       # rectangle and everything else unchanged. In such cases call
#       # - dom_nat_pft
#       # - zero_nonveg
#       # - std_elev
#       # - max_sat_area
#       # - soil_color (function not written, yet)
