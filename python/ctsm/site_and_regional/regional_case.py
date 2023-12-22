"""
This module includes the definition for a RegionalCase classs.
"""
# -- Import libraries
# -- Import Python Standard Libraries
import logging
import os

# -- 3rd party libraries
import numpy as np

# -- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase, USRDAT_DIR
from ctsm.utils import add_tag_to_filename

logger = logging.getLogger(__name__)


class RegionalCase(BaseCase):
    """
    A class to encapsulate regional cases.

    ...
    Attributes
    ----------
    lat1 : float
        first (left) latitude of a region.
    lat1 : float
        second (right) latitude of a region.
    lon1 : float
        first (bottom) longitude of a region.
    lon2 : float
        second (top) longitude of a region.
    reg_name: str -- default = None
        Region's name
    create_domain : bool
        flag for creating domain file
    create_surfdata : bool
        flag for creating surface dataset
    create_landuse : bool
        flag for creating landuse file
    create_datm : bool
        flag for creating DATM files
    create_user_mods : bool
        flag for creating user mods files and folders
    overwrite : bool
        flag for over-writing files if they already exist


    Methods
    -------
    create_tag
        Create a tag for this region which is either
        region's name or a combination of bounds of this
        region lat1-lat2_lon1-lon2

    create_domain_at_reg
        Create domain file at this region

    create_surfdata_at_reg
        Create surface dataset at this region

    create_landuse_at_reg
        Create landuse file at this region

    """

    def __init__(
        self,
        lat1,
        lat2,
        lon1,
        lon2,
        reg_name,
        create_domain,
        create_surfdata,
        create_landuse,
        create_datm,
        create_user_mods,
        out_dir,
        overwrite,
    ):
        """
        Initializes RegionalCase with the given arguments.
        """
        super().__init__(
            create_domain,
            create_surfdata,
            create_landuse,
            create_datm,
            create_user_mods,
            overwrite,
        )
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name
        self.out_dir = out_dir
        self.create_tag()

    def create_tag(self):
        """
        Create a tag for a region which is either the region name
        or
        the lat1-lat2_lon1-lon2 if the region name does not exist.
        """
        if self.reg_name:
            self.tag = self.reg_name
        else:
            self.tag = "{}-{}_{}-{}".format(
                str(self.lon1), str(self.lon2), str(self.lat1), str(self.lat2)
            )

    def create_domain_at_reg(self, indir, file):
        """
        Create domain file for this RegionalCase class.
        """

        # specify files
        fdomain_in = os.path.join(indir, file)
        fdomain_out = add_tag_to_filename(fdomain_in, self.tag)
        logger.info("fdomain_in:  %s", fdomain_in)
        logger.info("fdomain_out: %s", os.path.join(self.out_dir, fdomain_out))
        logger.info("Creating domain file at region: %s", self.tag)

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fdomain_in, "xc", "yc", "ni", "nj")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(nj=yind, ni=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fdomain_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.out_dir, fdomain_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fdomain_out) %s", wfile)
        f_in.close()
        f_out.close()

    def create_surfdata_at_reg(self, indir, file, user_mods_dir, specify_fsurf_out):
        """
        Create surface data file for this RegionalCase class.
        """

        logger.info("Creating surface dataset file at region: %s", self.tag)

        # specify files
        fsurf_in = os.path.join(indir, file)
        if specify_fsurf_out is None:
            fsurf_out = add_tag_to_filename(fsurf_in, self.tag, replace_res=True)
        else:
            fsurf_out = specify_fsurf_out

        logger.info("fsurf_in:  %s", fsurf_in)
        logger.info("fsurf_out: %s", os.path.join(self.out_dir, fsurf_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fsurf_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fsurf_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.out_dir, fsurf_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("created file (fsurf_out) %s", wfile)
        f_in.close()
        f_out.close()

        # write to user_nl_clm if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "fsurdat = '${}'".format(os.path.join(USRDAT_DIR, fsurf_out))
                self.write_to_file(line, nl_clm)

    def create_landuse_at_reg(self, indir, file, user_mods_dir):
        """
        Create land use data file for this RegionalCase class.
        """

        logger.info("Creating landuse file at region: %s", self.tag)

        # specify files
        fluse_in = os.path.join(indir, file)
        fluse_out = add_tag_to_filename(fluse_in, self.tag, replace_res=True)
        logger.info("fluse_in:  %s", fluse_in)
        logger.info("fluse_out: %s", os.path.join(self.out_dir, fluse_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fluse_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.out_dir, fluse_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fluse_out) %s", wfile)
        f_in.close()
        f_out.close()

        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                # line = "landuse = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                line = "flanduse_timeseries = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                self.write_to_file(line, nl_clm)
