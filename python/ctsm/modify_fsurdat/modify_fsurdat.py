"""
Run this code by using the following wrapper script:
/tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging

import numpy as np
import xarray as xr

from ctsm.utils import abort, get_git_sha, update_metadata, lon_range_0_to_360

logger = logging.getLogger(__name__)

class ModifyFsurdat:
    """
    Description
    -----------
    """

    def __init__(self, fsurdat_in, lon_1, lon_2, lat_1, lat_2, landmask_file):

        logger.info(
            'Opening fsurdat_in file to be modified: %s', fsurdat_in)
        self.file = xr.open_dataset(fsurdat_in)

        self.not_rectangle = self._get_not_rectangle(
            lon_1=lon_1, lon_2=lon_2,
            lat_1=lat_1, lat_2=lat_2,
            longxy=self.file.LONGXY, latixy=self.file.LATIXY)

        if landmask_file is not None:
            # overwrite self.not_rectangle with data from
            # user-specified .nc file in the .cfg file
            self._landmask_file = xr.open_dataset(landmask_file)
            rectangle = self._landmask_file.landmask
            self.not_rectangle = np.logical_not(rectangle)


    @staticmethod
    def _get_not_rectangle(lon_1, lon_2, lat_1, lat_2, longxy, latixy):
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
        not_rectangle = np.logical_not(rectangle)

        return not_rectangle


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

        # initialize 3D variable
        self.setvar('PCT_NAT_PFT', val=0, third_dim=-1)
        # set 3D variable value for dom_nat_pft
        self.setvar('PCT_NAT_PFT', val=100, third_dim=dom_nat_pft)

        # dictionary of 4d variables to loop over
        vars_4d = {'MONTHLY_LAI': lai,
                   'MONTHLY_SAI': sai,
                   'MONTHLY_HEIGHT_TOP': hgt_top,
                   'MONTHLY_HEIGHT_BOT': hgt_bot}
        for var, val in vars_4d.items():
            if val is not None:
                self.set_lai_sai_hgts(dom_nat_pft=dom_nat_pft,
                                      var=var, val=val)


    def set_lai_sai_hgts(self, dom_nat_pft, var, val):
        """
        Description
        -----------
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_nat_pft. Else do nothing.
        """
        if dom_nat_pft == 0:  # bare soil: var must equal 0
            val = [0] * 12
        if len(val) != 12:
            errmsg = 'Error: Variable should have exactly 12 ' \
                     'entries in the configure file: ' + var
            abort(errmsg)
        for mon in self.file.time - 1:  # loop over 12 months
            # set 4D variable to value for dom_nat_pft
            self.setvar(var, val[int(mon)], third_dim=dom_nat_pft,
                         fourth_dim=mon)


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


    def setvar(self, var, val, third_dim=None, fourth_dim=None):
        """
        Set variable var to value val in the *rectangle* domain, which
        is denoted as *other* in the *not_rectangle* domain.
        """

        if fourth_dim is None and third_dim is None:
            self.file[var] = self.file[var].where(
                self.not_rectangle, other=val)
        elif fourth_dim is None and third_dim == -1:
            self.file[var][:,:,:] = self.file[var][:,:,:].where(
                self.not_rectangle, other=val)
        elif fourth_dim is None and third_dim >= 0:
            self.file[var][third_dim,:,:] = \
                self.file[var][third_dim,:,:].where(
                self.not_rectangle, other=val)
        elif fourth_dim is not None and third_dim is not None:
            self.file[var][fourth_dim,third_dim,:,:] = \
                self.file[var][fourth_dim,third_dim,:,:].where(
                self.not_rectangle, other=val)
        else:
            errmsg = "Unexpected dimension setting for variable " + var + ". It's likely that the problem is in the third_dim and/or fourth_dim settings when calling setvar"
            abort(errmsg)


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
        # if pct_nat_veg had to be set to less than 100, then each special
        # landunit would have to receive a unique pct value rather than the
        # common value used here in pct_not_nat_veg = 0
        pct_nat_veg = 100  # do not change; works with pct_not_nat_veg = 0
        pct_not_nat_veg = 0  # do not change; works with pct_nat_veg = 100
        pct_sand = 43  # loam
        pct_clay = 18  # loam
        soil_color = 15  # loam
        organic = 0

        # 2D variables
        self.setvar('F0', max_inundated)
        self.setvar('FMAX', max_sat_area)
        self.setvar('STD_ELEV', std_elev)
        self.setvar('SLOPE', slope)
        self.setvar('zbedrock', zbedrock)
        self.setvar('SOIL_COLOR', soil_color)
        self.setvar('PFTDATA_MASK', pftdata_mask)
        self.setvar('LANDFRAC_PFT', landfrac_pft)
        self.setvar('PCT_WETLAND', pct_not_nat_veg)
        self.setvar('PCT_CROP', pct_not_nat_veg)
        self.setvar('PCT_LAKE', pct_not_nat_veg)
        self.setvar('PCT_URBAN', pct_not_nat_veg)
        self.setvar('PCT_GLACIER', pct_not_nat_veg)
        self.setvar('PCT_NATVEG', pct_nat_veg)

        for lev in self.file.nlevsoi:
            # set next three 3D variables to values representing loam
            self.setvar('PCT_SAND', val=pct_sand, third_dim=lev)
            self.setvar('PCT_CLAY', val=pct_clay, third_dim=lev)
            self.setvar('ORGANIC', val=organic, third_dim=lev)

        for crop in self.file.cft:
            cft_local = crop - (max(self.file.natpft) + 1)
            # initialize 3D variable; set outside the loop below
            self.setvar('PCT_CFT', val=0, third_dim=cft_local)

        # set 3D variable
        # NB. sum(PCT_CFT) must = 100 even though PCT_CROP = 0
        self.setvar('PCT_CFT', val=100, third_dim=0)
