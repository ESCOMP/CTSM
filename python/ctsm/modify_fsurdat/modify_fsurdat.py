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

    def __init__(self, fsurdat_in, lon_1, lon_2, lat_1, lat_2):

        print("Open file: " + fsurdat_in)
        self._file = xr.open_dataset(fsurdat_in)

        self._not_rectangle = self._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2,
            lat_1=lat_1, lat_2=lat_2,
            longxy=self._file.LONGXY, latixy=self._file.LATIXY)


    def _get_not_rectangle(self, lon_1, lon_2, lat_1, lat_2, longxy, latixy):
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

        if lat_1 <= lat_2:
            # rectangles overlap
            union_2 = np.logical_and(rectangle_3, rectangle_4)
        else:
            # rectangles don't overlap: one in the north, one in the south
            union_2 = np.logical_or(rectangle_3, rectangle_4)

        # union rectangles overlap
        rectangle = np.logical_and(union_1, union_2)
        _not_rectangle = np.logical_not(rectangle)


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


    def set_dom_nat_pft(self, dom_nat_pft, lai, sai, hgt_top, hgt_bot):
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
        """

        for pft in self._file.natpft:
            # initialize 3D variable; set outside the loop below
            self._file['PCT_NAT_PFT'][pft,:,:] = \
             self._file['PCT_NAT_PFT'][pft,:,:]. \
                  where(self._not_rectangle, other=0)

            var_name='MONTHLY_LAI'
            if len(lai) == 12:
                self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                      var_name=var_name, var=lai)
            elif len(lai) != 0:
                message = 'Error: The variable lai should have 12 ' \
                          'entries in the configure file: ' + var_name
                print(message)  # TODO do this via logging

            var_name='MONTHLY_SAI'
            if len(sai) == 12:
                self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                      var_name=var_name, var=sai)
            elif len(sai) != 0:
                message = 'Error: The variable sai should have 12 ' \
                          'entries in the configure file: ' + var_name
                print(message)  # TODO do this via logging

            var_name='MONTHLY_HEIGHT_TOP'
            if len(hgt_top) == 12:
                self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                      var_name=var_name, var=hgt_top)
            elif len(hgt_top) != 0:
                message = 'Error: Variable hgt_top should have 12 ' \
                          'entries in the configure file: ' + var_name
                print(message)  # TODO do this via logging

            var_name='MONTHLY_HEIGHT_BOT'
            if len(hgt_bot) == 12:
                self.set_lai_sai_hgts(pft=pft, dom_nat_pft=dom_nat_pft,
                                      var_name=var_name, var=hgt_bot)
            elif len(hgt_bot) != 0:
                message = 'Error: Variable hgt_bot should have 12 ' \
                          'entries in the configure file: ' + var_name
                print(message)  # TODO do this via logging

        # set 3D variable
        self._file['PCT_NAT_PFT'][dom_nat_pft,:,:] = \
         self._file['PCT_NAT_PFT'][dom_nat_pft,:,:]. \
              where(self._not_rectangle, other=100)


    def set_lai_sai_hgts(self, pft, dom_nat_pft, var_name, var):
        """
        Description
        -----------
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_nat_pft. Else do nothing.
        """
        if dom_nat_pft == 0:  # bare soil
            var = [0] * 12
        for mon in self._file.time - 1:  # loop over 12 months
            if pft == dom_nat_pft:
                # set 4D variable to value for dom_nat_pft
                self._file[var_name][mon,pft,:,:] = \
                 self._file[var_name][mon,pft,:,:]. \
                      where(self._not_rectangle, other=var[int(mon)])


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


    def set_idealized(self):
        """
        Description
        -----------
        Set fsurdat variables in a rectangle defined by lon/lat limits
        """

        # Overwrite in rectangle(s)
        # ------------------------
        # If idealized, the user makes changes to variables as follows.
        # "other" assigns the corresponding value in the rectangle.
        # Values outside the rectangle are preserved.
        # ------------------------

        # Default values
        zbedrock = 10
        max_sat_area = 0  # max saturated area
        max_inundated = 0  # max inundated fraction
        std_elev = 0  # standard deviation of elevation
        slope = 0  # mean topographic slope
        pftdata_mask = 1
        landfrac_pft = 1
        pct_nat_veg = 100
        pct_not_nat_veg = 0
        pct_cft = 100  # PCT_CROP = 0 but sum(PCT_CFT) must = 100
        pct_sand = 43  # loam
        pct_clay = 18  # loam
        soil_color = 15  # loam
        organic = 0

        # 2D variables
        self._file['F0'] = \
         self._file['F0'].where(self._not_rectangle, other=max_inundated)
        self._file['FMAX'] = \
         self._file['FMAX'].where(self._not_rectangle, other=max_sat_area)
        self._file['STD_ELEV'] = \
         self._file['STD_ELEV'].where(self._not_rectangle, other=std_elev)
        self._file['SLOPE'] = \
         self._file['SLOPE'].where(self._not_rectangle, other=slope)
        self._file['zbedrock'] = \
         self._file['zbedrock'].where(self._not_rectangle, other=zbedrock)
        self._file['SOIL_COLOR'] = \
         self._file['SOIL_COLOR'].where(self._not_rectangle, other=soil_color)
        self._file['PFTDATA_MASK'] = \
         self._file['PFTDATA_MASK'].where(self._not_rectangle, other=pftdata_mask)
        self._file['LANDFRAC_PFT'] = \
         self._file['LANDFRAC_PFT'].where(self._not_rectangle, other=landfrac_pft)
        self._file['PCT_WETLAND'] = \
         self._file['PCT_WETLAND'].where(self._not_rectangle, other=pct_not_nat_veg)
        self._file['PCT_CROP'] = \
         self._file['PCT_CROP'].where(self._not_rectangle, other=pct_not_nat_veg)
        self._file['PCT_LAKE'] = \
         self._file['PCT_LAKE'].where(self._not_rectangle, other=pct_not_nat_veg)
        self._file['PCT_URBAN'] = \
         self._file['PCT_URBAN'].where(self._not_rectangle, other=pct_not_nat_veg)
        self._file['PCT_GLACIER'] = \
         self._file['PCT_GLACIER'].where(self._not_rectangle, other=pct_not_nat_veg)
        self._file['PCT_NATVEG'] = \
         self._file['PCT_NATVEG'].where(self._not_rectangle, other=pct_nat_veg)

        for lev in self._file.nlevsoi:
            # set next three 3D variables to values representing loam
            self._file['PCT_SAND'][lev,:,:] = \
             self._file['PCT_SAND'][lev,:,:].where(self._not_rectangle, other=pct_sand)
            self._file['PCT_CLAY'][lev,:,:] = \
             self._file['PCT_CLAY'][lev,:,:].where(self._not_rectangle, other=pct_clay)
            self._file['ORGANIC'][lev,:,:] = \
             self._file['ORGANIC'][lev,:,:].where(self._not_rectangle, other=organic)

        for crop in self._file.cft:
            cft_local = crop - (max(self._file.natpft) + 1)
            # initialize 3D variable; set outside the loop below
            self._file['PCT_CFT'][cft_local,:,:] = \
             self._file['PCT_CFT'][cft_local,:,:]. \
                  where(self._not_rectangle, other=pct_not_nat_veg)

        # set 3D variable
        self._file['PCT_CFT'][0,:,:] = \
         self._file['PCT_CFT'][0,:,:].where(self._not_rectangle, other=pct_cft)
