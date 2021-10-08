"""
Run this code by using the following wrapper script:
/tools/modify_fsurdat/fsurdat_modifier

The wrapper script includes a full description and instructions.
"""

#  Import libraries
import os

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


    def land_swath(self, lon_in_1, lon_in_2, lat_in_1, lat_in_2,
                   dom_nat_pft, lai, sai, hgt_top, hgt_bot):
        """
        Description
        -----------
        Make a swath of land defined by lon/lat limits and make all else ocean.

        Arguments
        ---------
        lon_in_1:
            (int) westernmost edge of land swath
        lon_in_2:
            (int) easternmost edge of land swath
        lat_in_1:
            (int) southernmost edge of land swath
        lat_in_2:
            (int) northernmost edge of land swath
        dom_nat_pft:
            (int) user-selected PFT for all land area to be set to 100%
        lai:
            (float) user-defined leaf area index
        sai:
            (float) user-defined stem area index
        hgt_top:
            (float) user-defined height at top of vegetation canopy (m)
        hgt_bot:
            (float) user-defined height at bottom of vegetation canopy (m)
        """

        # Currently lon_in_1 & lon_in_2 are integers with required range 0-360.
        # If instead of requiring integer values, we decide to allow floats,
        # then ensure that lon ranges 0-360 in case user entered -180 to 180.
        lon_1 = lon_range_0_to_360(lon_in_1)
        lon_2 = lon_range_0_to_360(lon_in_2)

        # Find nearest longitude indices
        # TODO Placeholder. This is not really "nearest" but isel didn't work...
        temp = self.file.LONGXY[1,:]
        lon_idx_1 = int(min(temp.lsmlon.where(temp >= lon_1, drop=True)))
        lon_idx_2 = int(min(temp.lsmlon.where(temp >= lon_2, drop=True)))

        # Find nearest latitude indices
        # TODO Placeholder. This is not really "nearest" but isel didn't work...
        temp = self.file.LATIXY[:,1]
        lat_idx_1 = int(min(temp.lsmlat.where(temp >= lat_in_1, drop=True)))
        lat_idx_2 = int(min(temp.lsmlat.where(temp >= lat_in_2, drop=True)))

        # initialize to global values
        # set PFTDATA_MASK & LANDFRAC_PFT to 0 over ocean and to 1 over land
        #     PCT_NATVEG to 0 over ocean and to 100 over land
        #     PCT_WETLAND to 100 over ocean and to 0 over land
        #     the rest to 0 which means that there's no need to update
        #     other crop, lake, urban, glacier, and F0-related variables
        self.file['PFTDATA_MASK'][:,:] = 0  # partly overwrite below
        self.file['LANDFRAC_PFT'][:,:] = 0  # partly overwrite below
        self.file['PCT_NATVEG'][:,:] = 0  # partly overwrite below
        self.file['PCT_WETLAND'][:,:] = 100  # partly overwrite below
        self.file['PCT_CROP'][:,:] = 0
        self.file['PCT_LAKE'][:,:] = 0
        self.file['PCT_URBAN'][:,:] = 0
        self.file['PCT_GLACIER'][:,:] = 0
        self.file['F0'][:,:] = 0  # max inundated fraction
        self.file['FMAX'][:,:] = 0  # --max_sat_area can override this
        self.file['STD_ELEV'][:,:] = 0  # --std_elev can override this
        self.file['SLOPE'][:,:] = 0  # mean topographic slope
        # non-zero value seems appropriate; reasonable?
        self.file['zbedrock'][:,:] = 10
        # set next three to values representing loam
        self.file['SOIL_COLOR'][:,:] = 15
        self.file['PCT_SAND'][:,:,:] = 43
        self.file['PCT_CLAY'][:,:,:] = 18
        self.file['ORGANIC'][:,:,:] = 0
        # set remaining variables according to the dom_nat_pft option
        # but initialize to bare soil globally first
        self.file['PCT_CFT'][:,:,:] = 0
        self.file['PCT_CFT'][15,:,:] = 100  # required when PCT_CROP = 0
        self.file['PCT_NAT_PFT'][:,:,:] = 0
        self.file['PCT_NAT_PFT'][0,:,:] = 100
        self.file['MONTHLY_LAI'][:,:,:,:] = 0
        self.file['MONTHLY_SAI'][:,:,:,:] = 0
        self.file['MONTHLY_HEIGHT_TOP'][:,:,:,:] = 0
        self.file['MONTHLY_HEIGHT_BOT'][:,:,:,:] = 0
        if dom_nat_pft == -999:
            dom_nat_pft = 0  # default to bare soil
        elif dom_nat_pft > 0:
            # update lai, sai, and heights globally because they will be
            # used only if grid cell is vegetated
            self.file['MONTHLY_LAI'][:,dom_nat_pft,:,:] = lai
            self.file['MONTHLY_SAI'][:,dom_nat_pft,:,:] = sai
            self.file['MONTHLY_HEIGHT_TOP'][:,dom_nat_pft,:,:] = hgt_top
            self.file['MONTHLY_HEIGHT_BOT'][:,dom_nat_pft,:,:] = hgt_bot

        # overwrite with land swath
        if lon_idx_1 > lon_idx_2 and lat_idx_1 < lat_idx_2:
            # If only the lon indices are in descending order,
            # wrap around the 0-degree meridian
            self.file['PFTDATA_MASK'][lat_idx_1:lat_idx_2,:lon_idx_2] = 1
            self.file['LANDFRAC_PFT'][lat_idx_1:lat_idx_2,:lon_idx_2] = 1
            self.file['PCT_WETLAND'][lat_idx_1:lat_idx_2,:lon_idx_2] = 0
            self.file['PCT_NATVEG'][lat_idx_1:lat_idx_2,:lon_idx_2] = 100
            self.file['PCT_NAT_PFT'][0,lat_idx_1:lat_idx_2,:lon_idx_2] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,lat_idx_1:lat_idx_2,:lon_idx_2] = 100

            self.file['PFTDATA_MASK'][lat_idx_1:lat_idx_2,lon_idx_1:] = 1
            self.file['LANDFRAC_PFT'][lat_idx_1:lat_idx_2,lon_idx_1:] = 1
            self.file['PCT_WETLAND'][lat_idx_1:lat_idx_2,lon_idx_1:] = 0
            self.file['PCT_NATVEG'][lat_idx_1:lat_idx_2,lon_idx_1:] = 100
            self.file['PCT_NAT_PFT'][0,lat_idx_1:lat_idx_2,lon_idx_1:] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,lat_idx_1:lat_idx_2,lon_idx_1:] = 100
        elif lon_idx_1 < lon_idx_2 and lat_idx_1 > lat_idx_2:
            # If only the lat indices are in descending order, make two land
            # swaths, one in the north and one in the south rather than one
            # between the lat indices
            self.file['PFTDATA_MASK'][:lat_idx_2,lon_idx_1:lon_idx_2] = 1
            self.file['LANDFRAC_PFT'][:lat_idx_2,lon_idx_1:lon_idx_2] = 1
            self.file['PCT_WETLAND'][:lat_idx_2,lon_idx_1:lon_idx_2] = 0
            self.file['PCT_NATVEG'][:lat_idx_2,lon_idx_1:lon_idx_2] = 100
            self.file['PCT_NAT_PFT'][0,:lat_idx_2,lon_idx_1:lon_idx_2] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,:lat_idx_2,lon_idx_1:lon_idx_2] = 100

            self.file['PFTDATA_MASK'][lat_idx_1:,lon_idx_1:lon_idx_2] = 1
            self.file['LANDFRAC_PFT'][lat_idx_1:,lon_idx_1:lon_idx_2] = 1
            self.file['PCT_WETLAND'][lat_idx_1:,lon_idx_1:lon_idx_2] = 0
            self.file['PCT_NATVEG'][lat_idx_1:,lon_idx_1:lon_idx_2] = 100
            self.file['PCT_NAT_PFT'][0,lat_idx_1:,lon_idx_1:lon_idx_2] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,lat_idx_1:,lon_idx_1:lon_idx_2] = 100
        elif lon_idx_1 > lon_idx_2 and lat_idx_1 > lat_idx_2:
            # If both the lon and the lat indices are in descending order,
            # wrap around the 0-degree meridian AND make a land swath in the
            # north and another in the south rather than only one between the
            # lat indices
            self.file['PFTDATA_MASK'][:lat_idx_2,:lon_idx_2] = 1
            self.file['LANDFRAC_PFT'][:lat_idx_2,:lon_idx_2] = 1
            self.file['PCT_WETLAND'][:lat_idx_2,:lon_idx_2] = 0
            self.file['PCT_NATVEG'][:lat_idx_2,:lon_idx_2] = 100
            self.file['PCT_NAT_PFT'][0,:lat_idx_2,:lon_idx_2] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,:lat_idx_2,:lon_idx_2] = 100

            self.file['PFTDATA_MASK'][:lat_idx_2,lon_idx_1:] = 1
            self.file['LANDFRAC_PFT'][:lat_idx_2,lon_idx_1:] = 1
            self.file['PCT_WETLAND'][:lat_idx_2,lon_idx_1:] = 0
            self.file['PCT_NATVEG'][:lat_idx_2,lon_idx_1:] = 100
            self.file['PCT_NAT_PFT'][0,:lat_idx_2,lon_idx_1:] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,:lat_idx_2,lon_idx_1:] = 100

            self.file['PFTDATA_MASK'][lat_idx_1:,:lon_idx_2] = 1
            self.file['LANDFRAC_PFT'][lat_idx_1:,:lon_idx_2] = 1
            self.file['PCT_WETLAND'][lat_idx_1:,:lon_idx_2] = 0
            self.file['PCT_NATVEG'][lat_idx_1:,:lon_idx_2] = 100
            self.file['PCT_NAT_PFT'][0,lat_idx_1:,:lon_idx_2] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,lat_idx_1:,:lon_idx_2] = 100

            self.file['PFTDATA_MASK'][lat_idx_1:,lon_idx_1:] = 1
            self.file['LANDFRAC_PFT'][lat_idx_1:,lon_idx_1:] = 1
            self.file['PCT_WETLAND'][lat_idx_1:,lon_idx_1:] = 0
            self.file['PCT_NATVEG'][lat_idx_1:,lon_idx_1:] = 100
            self.file['PCT_NAT_PFT'][0,lat_idx_1,lon_idx_1:] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft,lat_idx_1,lon_idx_1:] = 100
        else:
            # Simple case: both lon and lat indices in ascending order
            self.file['PFTDATA_MASK'][lat_idx_1:lat_idx_2, \
                                      lon_idx_1:lon_idx_2] = 1
            self.file['LANDFRAC_PFT'][lat_idx_1:lat_idx_2, \
                                      lon_idx_1:lon_idx_2] = 1
            self.file['PCT_WETLAND'][lat_idx_1:lat_idx_2, \
                                     lon_idx_1:lon_idx_2] = 0
            self.file['PCT_NATVEG'][lat_idx_1:lat_idx_2, \
                                    lon_idx_1:lon_idx_2] = 100
            self.file['PCT_NAT_PFT'][0, \
                                     lat_idx_1:lat_idx_2, \
                                     lon_idx_1:lon_idx_2] = 0
            self.file['PCT_NAT_PFT'][dom_nat_pft, \
                                     lat_idx_1:lat_idx_2, \
                                     lon_idx_1:lon_idx_2] = 100
