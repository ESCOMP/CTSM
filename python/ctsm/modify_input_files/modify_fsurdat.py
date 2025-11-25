"""
Run this code by using the following wrapper script:
/tools/modify_input_files/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

import os
import logging

from math import isclose
import numpy as np
import xarray as xr

from ctsm.utils import abort, update_metadata
from ctsm.git_utils import get_ctsm_git_short_hash
from ctsm.longitude import Longitude

logger = logging.getLogger(__name__)


class ModifyFsurdat:
    """
    Description
    -----------
    """

    def __init__(
        self,
        my_data,
        *,
        lon_1,
        lon_2,
        lat_1,
        lat_2,
        landmask_file,
        lat_dimname,
        lon_dimname,
        lon_type=None,
    ):

        if lon_type is not None:
            lon_1 = Longitude(lon_1, lon_type)
            lon_2 = Longitude(lon_2, lon_type)

        self.numurbl = 3  # Number of urban density types
        self.file = my_data
        if "numurbl" in self.file.dims:
            self.numurbl = self.file.dims["numurbl"]
        else:
            abort("numurbl is not a dimension on the input surface dataset file and needs to be")

        self.rectangle = self._get_rectangle(
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            longxy=self.file.LONGXY,
            latixy=self.file.LATIXY,
        )

        if landmask_file is not None:
            # overwrite self.not_rectangle with data from
            # user-specified .nc file in the .cfg file
            landmask_ds = xr.open_dataset(landmask_file)
            self.rectangle = landmask_ds.mod_lnd_props.data
            # CF convention has dimension and coordinate variable names the same
            if lat_dimname is None:  # set to default
                lat_dimname = "lsmlat"
            if lon_dimname is None:  # set to default
                lon_dimname = "lsmlon"
            lsmlat = landmask_ds.dims[lat_dimname]
            lsmlon = landmask_ds.dims[lon_dimname]

            for row in range(lsmlat):  # rows from landmask file
                for col in range(lsmlon):  # cols from landmask file
                    errmsg = (
                        "landmask_ds.mod_lnd_props not 0 or 1 at "
                        + f"row, col, value = {row} {col} {self.rectangle[row, col]}"
                    )
                    assert isclose(self.rectangle[row, col], 0, abs_tol=1e-9) or isclose(
                        self.rectangle[row, col], 1, abs_tol=1e-9
                    ), errmsg

        self.not_rectangle = np.logical_not(self.rectangle)

    @classmethod
    def init_from_file(
        cls,
        *,
        fsurdat_in,
        lon_1,
        lon_2,
        lat_1,
        lat_2,
        landmask_file,
        lat_dimname,
        lon_dimname,
        lon_type,
    ):
        """Initialize a ModifyFsurdat object from file fsurdat_in"""
        logger.info("Opening fsurdat_in file to be modified: %s", fsurdat_in)
        my_file = xr.open_dataset(fsurdat_in)
        return cls(
            my_file,
            lon_1=lon_1,
            lon_2=lon_2,
            lat_1=lat_1,
            lat_2=lat_2,
            landmask_file=landmask_file,
            lat_dimname=lat_dimname,
            lon_dimname=lon_dimname,
            lon_type=lon_type,
        )

    @staticmethod
    def _get_rectangle(*, lon_1: Longitude, lon_2: Longitude, lat_1, lat_2, longxy, latixy):
        """
        Description
        -----------
        """

        # Get the longitudes in the correct format
        if not all(isinstance(x, Longitude) for x in [lon_1, lon_2]):
            raise TypeError("lon_1 and lon_2 must be of type Longitude")
        lon_type = 360
        lon_1 = lon_1.get(lon_type)
        lon_2 = lon_2.get(lon_type)

        # determine the rectangle(s)
        # TODO This is not really "nearest" for the edges but isel didn't work
        rectangle_1 = longxy >= lon_1
        rectangle_2 = longxy <= lon_2
        eps = np.finfo(np.float32).eps  # to avoid roundoff issue
        rectangle_3 = latixy >= (lat_1 - eps)
        rectangle_4 = latixy <= (lat_2 + eps)

        if lon_1 <= lon_2:
            # rectangles overlap
            union_1 = np.logical_and(rectangle_1, rectangle_2)
        else:
            # rectangles don't overlap: stradling the 0-degree meridian
            union_1 = np.logical_or(rectangle_1, rectangle_2)

        union_2 = None
        if lat_1 < -90 or lat_1 > 90 or lat_2 < -90 or lat_2 > 90:
            errmsg = "lat_1 and lat_2 need to be in the range -90 to 90"
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

    def get_urb_dens(self):
        """Get the number of urban density classes"""
        return self.numurbl

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
        title = "Modified fsurdat file"
        summary = "Modified fsurdat file"
        contact = "N/A"
        data_script = os.path.abspath(__file__) + " -- " + get_ctsm_git_short_hash()
        description = "Modified this file: " + fsurdat_in
        update_metadata(
            self.file,
            title=title,
            summary=summary,
            contact=contact,
            data_script=data_script,
            description=description,
        )

        # abort if output file already exists
        file_exists = os.path.exists(fsurdat_out)
        if file_exists:
            errmsg = "Output file already exists: " + fsurdat_out
            abort(errmsg)

        # mode 'w' overwrites file if it exists
        self.file.to_netcdf(path=fsurdat_out, mode="w", format="NETCDF3_64BIT")
        logger.info("Successfully created fsurdat_out: %s", fsurdat_out)
        self.file.close()

    def evenly_split_cropland(self):
        """
        Description
        -----------
        In rectangle selected by user (or default -90 to 90 and 0 to 360),
        replace fsurdat file's PCT_CFT with equal values for all crop types.
        """
        pct_cft = np.full_like(self.file["PCT_CFT"].values, 100 / self.file.dims["cft"])
        self.file["PCT_CFT"] = xr.DataArray(
            data=pct_cft, attrs=self.file["PCT_CFT"].attrs, dims=self.file["PCT_CFT"].dims
        )

    def set_dom_pft(self, *, dom_pft, lai, sai, hgt_top, hgt_bot):
        """
        Description
        -----------
        In rectangle selected by user (or default -90 to 90 and 0 to 360),
        replace fsurdat file's PCT_NAT_PFT or PCT_CFT with:
        - 100 for dom_pft selected by user
        - 0 for all other PFTs/CFTs
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_pft

        Arguments
        ---------
        dom_pft:
            (int) User's entry of PFT/CFT to be set to 100% everywhere
                  If user left this UNSET in the configure file, then it
                  will default to 0 (bare ground).
        lai:
            (float) User's entry of MONTHLY_LAI for their dom_pft
        sai:
            (float) User's entry of MONTHLY_SAI for their dom_pft
        hgt_top:
            (float) User's entry of MONTHLY_HEIGHT_TOP for their dom_pft
        hgt_bot:
            (float) User's entry of MONTHLY_HEIGHT_BOT for their dom_pft
        """

        # If dom_pft is a cft, add PCT_NATVEG to PCT_CROP in the rectangle
        # and remove same from PCT_NATVEG, i.e. set PCT_NATVEG = 0.
        if dom_pft > max(self.file.natpft):  # dom_pft is a cft (crop)
            self.file["PCT_CROP"] = self.file["PCT_CROP"] + self.file["PCT_NATVEG"].where(
                self.rectangle, other=0
            )
            self.setvar_lev0("PCT_NATVEG", 0)

            for cft in self.file.cft:
                cft_local = cft - (max(self.file.natpft) + 1)
                # initialize 3D variable; set outside the loop below
                self.setvar_lev1("PCT_CFT", val=0, lev1_dim=cft_local)

            # set 3D variable
            self.setvar_lev1("PCT_CFT", val=100, lev1_dim=dom_pft - (max(self.file.natpft) + 1))
        else:  # dom_pft is a pft (not a crop)
            for pft in self.file.natpft:
                # initialize 3D variable; set outside the loop below
                self.setvar_lev1("PCT_NAT_PFT", val=0, lev1_dim=pft)
            # set 3D variable value for dom_pft
            self.setvar_lev1("PCT_NAT_PFT", val=100, lev1_dim=dom_pft)

        # dictionary of 4d variables to loop over
        vars_4d = {
            "MONTHLY_LAI": lai,
            "MONTHLY_SAI": sai,
            "MONTHLY_HEIGHT_TOP": hgt_top,
            "MONTHLY_HEIGHT_BOT": hgt_bot,
        }
        for var, val in vars_4d.items():
            if val is not None:
                self.set_lai_sai_hgts(dom_pft=dom_pft, var=var, val=val)

    def check_varlist(
        self, settings, allow_uppercase_vars=False, source="input settings dictionary"
    ):
        """
        Check a list of variables from a dictionary of settings
        """
        settings_return = {}
        varlist = settings.keys()
        for var in varlist:
            varname = var
            val = settings[varname]
            if not var in self.file:
                if not allow_uppercase_vars:
                    errmsg = "Error: Variable " + varname + " is NOT in the " + source
                    abort(errmsg)
                if not varname.upper() in self.file:
                    errmsg = "Error: Variable " + varname.upper() + " is NOT in the " + source
                    abort(errmsg)
                varname = varname.upper()

            settings_return[varname] = val
            #
            # Check that dimensions are as expected
            #
            if len(self.file[varname].dims) == 2:
                if not isinstance(val, float):
                    abort(
                        "For 2D vars, there should only be a single value for variable = "
                        + varname
                        + " in "
                        + source
                    )
            elif len(self.file[varname].dims) >= 3:
                dim1 = int(self.file.sizes[self.file[varname].dims[0]])
                if not isinstance(val, list):
                    abort(
                        "For higher dimensional vars, the variable needs to be expressed "
                        + "as a list of values of the dimension size = "
                        + str(dim1)
                        + " for variable="
                        + varname
                        + " in "
                        + source
                    )
                if len(val) != dim1:
                    abort(
                        "Variable "
                        + varname
                        + " is "
                        + str(len(val))
                        + " is of the wrong size. It should be = "
                        + str(dim1)
                        + " in "
                        + source
                    )
        return settings_return

    def set_varlist(self, settings, cfg_path="unknown-config-file"):
        """
        Set a list of variables from a dictionary of settings
        """
        for var in settings.keys():
            if var in self.file:
                if len(self.file[var].dims) == 2:
                    if not isinstance(settings[var], float):
                        abort(
                            "For 2D vars, there should only be a single value for variable = " + var
                        )
                    self.setvar_lev0(var, settings[var])
                elif len(self.file[var].dims) == 3:
                    dim1 = int(self.file.sizes[self.file[var].dims[0]])
                    vallist = settings[var]
                    if not isinstance(vallist, list):
                        abort(
                            "For higher dimensional vars, there must be a list of values "
                            + "for variable= "
                            + var
                            + " from the config file = "
                            + cfg_path
                        )
                    if len(vallist) != dim1:
                        abort(
                            "Variable " + var + " is of the wrong size. It should be = " + str(dim1)
                        )
                    for lev1 in range(dim1):
                        self.setvar_lev1(var, vallist[lev1], lev1_dim=lev1)
                elif len(self.file[var].dims) == 4:
                    dim_lev1 = int(self.file.sizes[self.file[var].dims[1]])
                    dim_lev2 = int(self.file.sizes[self.file[var].dims[0]])
                    vallist = settings[var]
                    for lev1 in range(dim_lev1):
                        for lev2 in range(dim_lev2):
                            self.setvar_lev2(var, vallist[lev2], lev1_dim=lev1, lev2_dim=lev2)
                else:
                    abort(
                        "Error: Variable "
                        + var
                        + " is a higher dimension than currently allowed = "
                        + str(self.file[var].dims)
                    )
            else:
                errmsg = "Error: Variable " + var + " is NOT in the file"
                abort(errmsg)

    def set_lai_sai_hgts(self, dom_pft, var, val):
        """
        Description
        -----------
        If user has specified lai, sai, hgt_top, hgt_bot, replace these with
        values selected by the user for dom_pft. Else do nothing.
        """
        months = int(max(self.file.time))  # 12 months
        if dom_pft == 0:  # bare soil: var must equal 0
            val = [0] * months
        if len(val) != months:
            errmsg = (
                "Error: Variable should have exactly "
                + months
                + " entries in the configure file: "
                + var
            )
            abort(errmsg)
        for mon in self.file.time - 1:  # loop over 12 months
            # set 4D variable to value for dom_pft
            self.setvar_lev2(var, val[int(mon)], lev1_dim=dom_pft, lev2_dim=mon)

    def zero_nonveg(self):
        """
        Description
        -----------
        Set all landunit weights to 0 except the natural vegetation landunit.
        Set that one to 100%.
        """

        self.setvar_lev0("PCT_NATVEG", 100)
        self.setvar_lev0("PCT_CROP", 0)
        self.setvar_lev0("PCT_LAKE", 0)
        self.setvar_lev0("PCT_WETLAND", 0)
        self.setvar_lev0("PCT_URBAN", 0)
        self.setvar_lev0("PCT_GLACIER", 0)
        self.setvar_lev0("PCT_OCEAN", 0)

    def setvar_lev0(self, var, val):
        """
        Sets 2d variable var to value val in user-defined rectangle,
        defined as "other" in the function
        """
        self.file[var] = self.file[var].where(self.not_rectangle, other=val)

    def setvar_lev1(self, var, val, lev1_dim):
        """
        Sets 3d variable var to value val in user-defined rectangle,
        defined as "other" in the function
        """
        self.file[var][lev1_dim, ...] = self.file[var][lev1_dim, ...].where(
            self.not_rectangle, other=val
        )

    def setvar_lev2(self, var, val, lev1_dim, lev2_dim):
        """
        Sets 4d variable var to value val in user-defined rectangle,
        defined as "other" in the function
        """
        self.file[var][lev2_dim, lev1_dim, ...] = self.file[var][lev2_dim, lev1_dim, ...].where(
            self.not_rectangle, other=val
        )

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
        std_elev = 0  # standard deviation of elevation
        slope = 0  # mean topographic slope
        landfrac_pft = 1
        landfrac_mksurfdata = 1
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
        self.setvar_lev0("FMAX", max_sat_area)
        self.setvar_lev0("STD_ELEV", std_elev)
        self.setvar_lev0("SLOPE", slope)
        self.setvar_lev0("zbedrock", zbedrock)
        self.setvar_lev0("SOIL_COLOR", soil_color)
        self.setvar_lev0("LANDFRAC_PFT", landfrac_pft)
        self.setvar_lev0("LANDFRAC_MKSURFDATA", landfrac_mksurfdata)
        self.setvar_lev0("PCT_WETLAND", pct_not_nat_veg)
        self.setvar_lev0("PCT_CROP", pct_not_nat_veg)
        self.setvar_lev0("PCT_LAKE", pct_not_nat_veg)
        self.setvar_lev0("PCT_URBAN", pct_not_nat_veg)
        self.setvar_lev0("PCT_GLACIER", pct_not_nat_veg)
        self.setvar_lev0("PCT_OCEAN", pct_not_nat_veg)
        self.setvar_lev0("PCT_NATVEG", pct_nat_veg)

        for lev in self.file.nlevsoi:
            # set next three 3D variables to values representing loam
            self.setvar_lev1("PCT_SAND", val=pct_sand, lev1_dim=lev)
            self.setvar_lev1("PCT_CLAY", val=pct_clay, lev1_dim=lev)
            self.setvar_lev1("ORGANIC", val=organic, lev1_dim=lev)

        for crop in self.file.cft:
            cft_local = crop - (max(self.file.natpft) + 1)
            # initialize 3D variable; set outside the loop below
            self.setvar_lev1("PCT_CFT", val=0, lev1_dim=cft_local)

        # set 3D variable
        # NB. sum(PCT_CFT) must = 100 even though PCT_CROP = 0
        self.setvar_lev1("PCT_CFT", val=100, lev1_dim=0)
