"""
This module includes the definition for a RegionalCase classs.
"""
# -- Import libraries
# -- Import Python Standard Libraries
import logging

# -- 3rd party libraries
import numpy as np

# -- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase

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
    ):
        """
        Initializes RegionalCase with the given arguments.
        """
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm)
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name

    def create_tag(self):
        """
        Create a tag for a region which is either the region name
        or 
        the lat1-lat2_lon1-lon2 if the region name does not exist.
        """
        if self.reg_name:
            self.tag = self.reg_name
        else:
            self.tag = (
                str(self.lon1)
                + "-"
                + str(self.lon2)
                + "_"
                + str(self.lat1)
                + "-"
                + str(self.lat2)
            )

    def create_domain_at_reg(self):
        """
        Create domain file for a region.
        """
        logger.info(
            "----------------------------------------------------------------------"
        )
        logger.info("Creating domain file at region:" + self.tag)

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(self.fdomain_in, "xc", "yc", "ni", "nj")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(nj=yind, ni=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fdomain_in

        wfile = self.fdomain_out
        # mode 'w' overwrites file
        f_out.to_netcdf(path=wfile, mode="w", format='NETCDF3_64BIT')
        logger.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f_in.close()
        f_out.close()

    def create_surfdata_at_reg(self):
        """
        Create surface dataset for a region.
        """
        logger.info(
            "----------------------------------------------------------------------"
        )
        logger.info("Creating surface dataset file at region:" + self.tag)

        # create 1d coordinate variables to enable sel() method
        filename = self.fsurf_in
        f_in = self.create_1d_coord(filename, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fsurf_in

        # mode 'w' overwrites file
        f_out.to_netcdf(path=self.fsurf_out, mode="w", format='NETCDF3_64BIT')
        logger.info("Successfully created file (fsurf_out)" + self.fsurf_out)

        f_in.close()
        f_out.close()

    def create_landuse_at_reg(self):
        """
        Create landuse file for a region.
        """
        logger.info(
            "----------------------------------------------------------------------"
        )
        logger.info("Creating landuse file at region:" + self.tag)

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(
            self.fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat"
        )
        lat = f_in["lat"]
        lon = f_in["lon"]

        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f_out = f_in.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fluse_in

        wfile = self.fluse_out
        # mode 'w' overwrites file
        f_out.to_netcdf(path=wfile, mode="w", format='NETCDF3_64BIT')
        logger.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f_in.close()
        f_out.close()
