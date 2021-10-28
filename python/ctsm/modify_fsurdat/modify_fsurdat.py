"""
Run this code by using the following wrapper script:
/tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

import os

import numpy as np
import xarray as xr

from ctsm.utils import abort, get_git_sha, update_metadata, lon_range_0_to_360

class ModifyFsurdat:
    """

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


    def set_dom_nat_pft(self, dom_nat_pft, lai, sai, hgt_top, hgt_bot,
                        not_rectangle):
        """
        Description
        -----------
        In rectangle selected by user (or default -90 to 90 and 0 to 360),
        replace fsurdat file's PCT_NAT_PFT with:
        - 100 for dom_nat_pft selected by user
        - 0 for all other non-crop PFTs
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_nat_pft

        Arguments
        ---------
        dom_nat_pft:
            (int) User's entry of PFT to be set to 100% everywhere
        lai:
            (float) User's entry of MONTHLY_LAI for their dom_nat_pft
        sai:
            (float) User's entry of MONTHLY_SAI for their dom_nat_pft
        hgt_top:
            (float) User's entry of MONTHLY_HEIGHT_TOP for their dom_nat_pft
        hgt_bot:
            (float) User's entry of MONTHLY_HEIGHT_BOT for their dom_nat_pft
        not_rectangle:
            (xarray dataarray) Inverse of land rectangle selected by user
        """

        for pft in self.file.natpft:
            # initialize 3D variable; set outside the loop below
            self.file['PCT_NAT_PFT'][pft,:,:] = \
             self.file['PCT_NAT_PFT'][pft,:,:]. \
                  where(not_rectangle, other=0)

            self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                  var_name='MONTHLY_LAI', var=lai,
                                  not_rectangle=not_rectangle)
            self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                  var_name='MONTHLY_SAI', var=sai,
                                  not_rectangle=not_rectangle)
            self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                  var_name='MONTHLY_HEIGHT_TOP', var=hgt_top,
                                  not_rectangle=not_rectangle)
            self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                  var_name='MONTHLY_HEIGHT_BOT', var=hgt_bot,
                                  not_rectangle=not_rectangle)

        # set 3D variable
        self.file['PCT_NAT_PFT'][dom_nat_pft,:,:] = \
         self.file['PCT_NAT_PFT'][dom_nat_pft,:,:]. \
              where(not_rectangle, other=100)


    def set_lai_sai_hgts(self, pft, dom_nat_pft, var_name, var, not_rectangle):
        """
        Description
        -----------
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_nat_pft. Else do nothing.
        """
        if len(var) == 12:
            for mon in self.file.time - 1:  # loop over 12 months
                # initialize 4D variable; set in "if pft" below
                self.file[var_name][mon,pft,:,:] = \
                 self.file[var_name][mon,pft,:,:]. \
                      where(not_rectangle, other=0)
                if pft == dom_nat_pft:
                    # set 4D variable to value for dom_nat_pft
                    self.file[var_name][mon,pft,:,:] = \
                     self.file[var_name][mon,pft,:,:]. \
                          where(not_rectangle, other=var[int(mon)])
        elif len(var) != 0:
            message = 'Error: This variable should have 12 entries in the ' \
                      'configure file: ' + var_name
            print(message)  # TODO do this via logging


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


    def set_in_rectangle(self, _idealized, _lon_in_1, _lon_in_2, _lat_in_1,
                        _lat_in_2, _dom_nat_pft, _lai, _sai, _hgt_top,
                        _hgt_bot, _zero_nonveg, _std_elev, _soil_color,
                        _max_sat_area):
        """
        Description
        -----------
        Set fsurdat variables in a rectangle defined by lon/lat limits

        Arguments
        ---------
        _lon_in_1:
            (int) westernmost edge of rectangle
        _lon_in_2:
            (int) easternmost edge of rectangle
        _lat_in_1:
            (int) southernmost edge of rectangle
        _lat_in_2:
            (int) northernmost edge of rectangle
        _dom_nat_pft:
            (int) user-selected PFT to be set to 100% in rectangle
        _lai:
            (list of floats) user-defined leaf area index
        _sai:
            (list of floats) user-defined stem area index
        _hgt_top:
            (list of floats) user-defined height at top of vegetation
            canopy (m)
        _hgt_bot:
            (list of floats) user-defined height at bottom of vegetation
            canopy (m)
        _zero_nonveg:
        _std_elev:
        _soil_color:
        _max_sat_area:
        """

        # Currently type(lon_in_*) = int with required range 0-360.
        # If instead of requiring integer values, we decide to allow floats,
        # then ensure that lon ranges 0-360 in case user entered -180 to 180.
        lon_1 = lon_range_0_to_360(_lon_in_1)
        lon_2 = lon_range_0_to_360(_lon_in_2)

        # determine the rectangle(s)
        # TODO This is not really "nearest" for the edges but isel didn't work
        rectangle_1 = (self.file.LONGXY >= lon_1)
        rectangle_2 = (self.file.LONGXY <= lon_2)
        eps = np.finfo(np.float32).eps  # to avoid roundoff issue
        rectangle_3 = (self.file.LATIXY >= (_lat_in_1 - eps))
        rectangle_4 = (self.file.LATIXY <= (_lat_in_2 + eps))

        if lon_1 <= lon_2:
            # rectangles overlap
            union_1 = np.logical_and(rectangle_1, rectangle_2)
        else:
            # rectangles don't overlap (stradling the 0-degree meridian)
            union_1 = np.logical_or(rectangle_1, rectangle_2)

        if _lat_in_1 <= _lat_in_2:
            # rectangles overlap
            union_2 = np.logical_and(rectangle_3, rectangle_4)
        else:
            # rectangles don't overlap (one in the north, on in the south)
            union_2 = np.logical_or(rectangle_3, rectangle_4)

        # rectangles overlap
        rectangle = np.logical_and(union_1, union_2)
        not_rectangle = np.logical_not(rectangle)

        # Overwrite in rectangle(s)
        # ------------------------
        # If idealized, the user makes changes to variables as follows.
        # "other" assigns the corresponding value in the rectangle.
        # Values outside the rectangle are preserved.
        # ------------------------
        if _idealized:
            if _dom_nat_pft is None:  # default to bare soil when not user-set
                _dom_nat_pft = 0
            if _dom_nat_pft == 0:  # values corresponding to bare soil
                _lai = [0] * 12
                _sai = [0] * 12
                _hgt_top = [0] * 12
                _hgt_bot = [0] * 12
            if _soil_color is None:  # default to loam when not user-set
                _soil_color = 15
            if _std_elev is None:  # other default values when not user-set
                _std_elev = 0
            if _max_sat_area is None:
                _max_sat_area = 0

            # 2D variables
            # max inundated fraction
            self.file['F0'] = \
             self.file['F0'].where(not_rectangle, other=0)
            # max saturated area
            self.file['FMAX'] = \
             self.file['FMAX'].where(not_rectangle, other=_max_sat_area)
            # standard deviation of elevation
            self.file['STD_ELEV'] = \
             self.file['STD_ELEV'].where(not_rectangle, other=_std_elev)
            # mean topographic slope
            self.file['SLOPE'] = \
             self.file['SLOPE'].where(not_rectangle, other=0)
            # zbedrock
            self.file['zbedrock'] = \
             self.file['zbedrock'].where(not_rectangle, other=10)
            # value representing loam
            self.file['SOIL_COLOR'] = \
             self.file['SOIL_COLOR'].where(not_rectangle, other=_soil_color)

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

            for crop in self.file.cft:
                cft_local = crop - (max(self.file.natpft) + 1)
                # initialize 3D variable; set outside the loop below
                self.file['PCT_CFT'][cft_local,:,:] = \
                 self.file['PCT_CFT'][cft_local,:,:]. \
                      where(not_rectangle, other=0)

            # set 3D variable
            # required although PCT_CROP = 0: the sum of all crops must = 100
            self.file['PCT_CFT'][0,:,:] = \
             self.file['PCT_CFT'][0,:,:].where(not_rectangle, other=100)

            # set 3D and 4D variables PCT_NAT_PFT, MONLTHLY_LAI, MONTHLY_SAI,
            # MONTHLY_HEIGHT_TOP, MONTHLY_HEIGHT_TOP
            self.set_dom_nat_pft(dom_nat_pft=_dom_nat_pft, lai=_lai, sai=_sai,
                                 hgt_top=_hgt_top, hgt_bot=_hgt_bot,
                                 not_rectangle=not_rectangle)

        else:  # not idealized
            # Changing specific vars in the rectangle and leaving everything
            # else unchanged
            if _dom_nat_pft is not None:
                self.set_dom_nat_pft(dom_nat_pft=_dom_nat_pft,
                                     lai=_lai, sai=_sai,
                                     hgt_top=_hgt_top, hgt_bot=_hgt_bot,
                                     not_rectangle=not_rectangle)
            if _max_sat_area is not None:
                # max saturated area
                self.file['FMAX'] = \
                 self.file['FMAX'].where(not_rectangle, other=_max_sat_area)
            if _std_elev is not None:
                # standard deviation of elevation
                self.file['STD_ELEV'] = \
                 self.file['STD_ELEV'].where(not_rectangle, other=_std_elev)
            if _soil_color is not None:
                # value representing loam
                self.file['SOIL_COLOR'] = \
                 self.file['SOIL_COLOR'].where(not_rectangle,
                                               other=_soil_color)
            # TODO Next also zero_nonveg or is idealized option sufficient?
