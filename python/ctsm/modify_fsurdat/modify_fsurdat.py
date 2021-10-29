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
        self._file = xr.open_dataset(fsurdat_in)


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
        data_script = os.path.abspath(__file__) + " -- " + get_git_sha()
        description = 'Modified this file: ' + fsurdat_in
        update_metadata(self._file, title=title, summary=summary,
                        contact=contact, data_script=data_script,
                        description=description)

        # abort if output file already exists
        file_exists = os.path.exists(fsurdat_out)
        if file_exists:
            errmsg = 'Output file already exists: ' + fsurdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        self._file.to_netcdf(path=fsurdat_out, mode='w',
                            format="NETCDF3_64BIT")
        print('Successfully created fsurdat_out: ' + fsurdat_out)
        self._file.close()


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

        for pft in self._file.natpft:
            # initialize 3D variable; set outside the loop below
            self._file['PCT_NAT_PFT'][pft,:,:] = \
             self._file['PCT_NAT_PFT'][pft,:,:]. \
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
        self._file['PCT_NAT_PFT'][dom_nat_pft,:,:] = \
         self._file['PCT_NAT_PFT'][dom_nat_pft,:,:]. \
              where(not_rectangle, other=100)


    def set_lai_sai_hgts(self, pft, dom_nat_pft, var_name, var, not_rectangle):
        """
        Description
        -----------
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_nat_pft. Else do nothing.
        """
        if len(var) == 12:
            if dom_nat_pft == 0:  # bare soil
                var = [0] * 12
            for mon in self._file.time - 1:  # loop over 12 months
                if pft == dom_nat_pft:
                    # set 4D variable to value for dom_nat_pft
                    self._file[var_name][mon,pft,:,:] = \
                     self._file[var_name][mon,pft,:,:]. \
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

        self._file['PCT_NATVEG'][:,:] = 100
        self._file['PCT_CROP'][:,:] = 0
        self._file['PCT_LAKE'][:,:] = 0
        self._file['PCT_WETLAND'][:,:] = 0
        self._file['PCT_URBAN'][:,:] = 0
        self._file['PCT_GLACIER'][:,:] = 0


    def set_in_rectangle(self, idealized, lon_in_1, lon_in_2, lat_in_1,
                        lat_in_2, dom_nat_pft, lai, sai, hgt_top,
                        hgt_bot, zero_nonveg, std_elev, soil_color,
                        max_sat_area):
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
            (list of floats) user-defined height at top of vegetation
            canopy (m)
        hgt_bot:
            (list of floats) user-defined height at bottom of vegetation
            canopy (m)
        zero_nonveg:
        std_elev:
        soil_color:
        max_sat_area:
        """

        # Currently type(lon_in_*) = int with required range 0-360.
        # If instead of requiring integer values, we decide to allow floats,
        # then ensure that lon ranges 0-360 in case user entered -180 to 180.
        lon_1 = lon_range_0_to_360(lon_in_1)
        lon_2 = lon_range_0_to_360(lon_in_2)

        # determine the rectangle(s)
        # TODO This is not really "nearest" for the edges but isel didn't work
        rectangle_1 = (self._file.LONGXY >= lon_1)
        rectangle_2 = (self._file.LONGXY <= lon_2)
        eps = np.finfo(np.float32).eps  # to avoid roundoff issue
        rectangle_3 = (self._file.LATIXY >= (lat_in_1 - eps))
        rectangle_4 = (self._file.LATIXY <= (lat_in_2 + eps))

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

        # Overwrite in rectangle(s)
        # ------------------------
        # If idealized, the user makes changes to variables as follows.
        # "other" assigns the corresponding value in the rectangle.
        # Values outside the rectangle are preserved.
        # ------------------------
        if idealized:
            if dom_nat_pft is None:  # not user-set
                dom_nat_pft = 0  # default to bare soil
            if soil_color is None:  # not user-set
                soil_color = 15  # default to loam
            if std_elev is None:  # not user-set
                std_elev = 0
            if max_sat_area is None:  # not user-set
                max_sat_area = 0

            # 2D variables
            # max inundated fraction
            self._file['F0'] = \
             self._file['F0'].where(not_rectangle, other=0)
            # max saturated area
            self._file['FMAX'] = \
             self._file['FMAX'].where(not_rectangle, other=max_sat_area)
            # standard deviation of elevation
            self._file['STD_ELEV'] = \
             self._file['STD_ELEV'].where(not_rectangle, other=std_elev)
            # mean topographic slope
            self._file['SLOPE'] = \
             self._file['SLOPE'].where(not_rectangle, other=0)
            # zbedrock
            self._file['zbedrock'] = \
             self._file['zbedrock'].where(not_rectangle, other=10)
            # value representing loam
            self._file['SOIL_COLOR'] = \
             self._file['SOIL_COLOR'].where(not_rectangle, other=soil_color)

            self._file['PFTDATA_MASK'] = \
             self._file['PFTDATA_MASK'].where(not_rectangle, other=1)
            self._file['LANDFRAC_PFT'] = \
             self._file['LANDFRAC_PFT'].where(not_rectangle, other=1)
            self._file['PCT_WETLAND'] = \
             self._file['PCT_WETLAND'].where(not_rectangle, other=0)
            self._file['PCT_CROP'] = \
             self._file['PCT_CROP'].where(not_rectangle, other=0)
            self._file['PCT_LAKE'] = \
             self._file['PCT_LAKE'].where(not_rectangle, other=0)
            self._file['PCT_URBAN'] = \
             self._file['PCT_URBAN'].where(not_rectangle, other=0)
            self._file['PCT_GLACIER'] = \
             self._file['PCT_GLACIER'].where(not_rectangle, other=0)
            self._file['PCT_NATVEG'] = \
             self._file['PCT_NATVEG'].where(not_rectangle, other=100)

            for lev in self._file.nlevsoi:
                # set next three 3D variables to values representing loam
                self._file['PCT_SAND'][lev,:,:] = \
                 self._file['PCT_SAND'][lev,:,:].where(not_rectangle, other=43)
                self._file['PCT_CLAY'][lev,:,:] = \
                 self._file['PCT_CLAY'][lev,:,:].where(not_rectangle, other=18)
                self._file['ORGANIC'][lev,:,:] = \
                 self._file['ORGANIC'][lev,:,:].where(not_rectangle, other=0)

            for crop in self._file.cft:
                cft_local = crop - (max(self._file.natpft) + 1)
                # initialize 3D variable; set outside the loop below
                self._file['PCT_CFT'][cft_local,:,:] = \
                 self._file['PCT_CFT'][cft_local,:,:]. \
                      where(not_rectangle, other=0)

            # set 3D variable
            # required although PCT_CROP = 0: the sum of all crops must = 100
            self._file['PCT_CFT'][0,:,:] = \
             self._file['PCT_CFT'][0,:,:].where(not_rectangle, other=100)

            # set 3D and 4D variables PCT_NAT_PFT, MONLTHLY_LAI, MONTHLY_SAI,
            # MONTHLY_HEIGHT_TOP, MONTHLY_HEIGHT_TOP
            self.set_dom_nat_pft(dom_nat_pft=dom_nat_pft,
                                 lai=lai, sai=sai,
                                 hgt_top=hgt_top, hgt_bot=hgt_bot,
                                 not_rectangle=not_rectangle)

        else:  # not idealized
            # Changing specific vars in the rectangle and leaving everything
            # else unchanged
            if dom_nat_pft is not None:
                self.set_dom_nat_pft(dom_nat_pft=dom_nat_pft,
                                     lai=lai, sai=sai,
                                     hgt_top=hgt_top, hgt_bot=hgt_bot,
                                     not_rectangle=not_rectangle)
            if max_sat_area is not None:
                # max saturated area
                self._file['FMAX'] = \
                 self._file['FMAX'].where(not_rectangle, other=max_sat_area)
            if std_elev is not None:
                # standard deviation of elevation
                self._file['STD_ELEV'] = \
                 self._file['STD_ELEV'].where(not_rectangle, other=std_elev)
            if soil_color is not None:
                # value representing loam
                self._file['SOIL_COLOR'] = \
                 self._file['SOIL_COLOR'].where(not_rectangle,
                                                other=soil_color)
            # TODO Next also zero_nonveg or is idealized option sufficient?
