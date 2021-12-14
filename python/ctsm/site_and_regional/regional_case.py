"""
This module includes the definition for a RegionalCase classs.
"""
#-- Import libraries
#-- Import Python Standard Libraries
import logging

#-- 3rd party libraries
import numpy as np
import xarray as xr

#-- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase

logger = logging.getLogger(__name__)

class RegionalCase(BaseCase):
    """
    A case to encapsulate regional cases.

    ...
    Attributes
    ----------
    plat : float
        latitude
    plon : float
        longitude
    site_name: str -- default = None
        Site name



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
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm)
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name

    def create_tag(self):
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
        #logging.debug ("Creating domain file at region"+ self.lon1.__str__()+"-"+self.lat2.__str__()+" "+self.lat1.__str__()+"-"+self.lat2.__str__())
        logging.info("Creating domain file at region:"+ self.tag)
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
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f_in.close()
        f_out.close()

    def create_surfdata_at_reg(self):
        #logging.debug ("Creating surface dataset file at region"+ self.lon1.__str__()+"-"+self.lat2.__str__()+" "+self.lat1.__str__()+"-"+self.lat2.__str__())
        logging.info("Creating surface dataset file at region:"+ self.tag)
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
        f_out.to_netcdf(path=self.fsurf_out, mode="w")
        logging.info("created file (fsurf_out)" + self.fsurf_out)
        # f1.close();
        f_in.close()
        f_out.close()

    def create_landuse_at_reg(self):
        #logging.debug ("Creating landuse file at region"+ self.lon1.__str__()+"-"+self.lat2.__str__()+" "+self.lat1.__str__()+"-"+self.lat2.__str__())
        logging.info("Creating landuse file at region:"+ self.tag)
        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(self.fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
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
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f_in.close()
        f_out.close()


