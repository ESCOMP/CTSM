import logging

import numpy as np
import xarray as xr

from ctsm.site_and_regional.base_case import BaseCase

logger = logging.getLogger(__name__)

class RegionalCase(BaseCase):
    """
    A case to encapsulate regional cases.
    ...
    Attributes
    ----------
    lat1 : float
        start latitude
    lat2 : float
        end latitude
    lon1 : float
        start longitude
    lon2 : float
        end longitude
    reg_name: str -- default = None
        region name
    tag : str
        ending tag for output file naming
    fluse_out : str
        file name of output subset land use file
    fluse_in : str
        file name of input land use file to subset
    fsurf_out : str
        file name of output subset surface data file
    fsurf_in : str
        file name of input surface data to subset
    fdomain_out : str
        file name of output domain subset domain file
    fdomain_in : str
        file name of input domain file to subset

    Methods
    -------
    create_tag
        create a tag for a region which is the region name
        or the "lon1-lon2-lat1-lat2" format if the region name does not exist.
    create_domain_at_reg
        Create domain file for a region
    create_landuse_at_reg:
        Create landuse file file for a region
    create_surfdata_at_reg:
        Create surface dataset file for a region
    create_datmdomain_at_reg:
        Create DATM domain file for a region
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
            output_dir,
    ):
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm, create_user_mods)
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name
        self.output_dir = output_dir
        self.tag = None
        self.fluse_out = None
        self.fluse_in = None
        self.fsurf_out = None
        self.fsurf_in = None
        self.fdomain_out = None
        self.fdomain_in = None

    def create_tag(self):
        if self.reg_name:
            self.tag = self.reg_name
        else:
            self.tag = "{}-{}_{}-{}".format(str(self.lon1), str(self.lon2), str(self.lat1), str(self.lat2))

    def create_domain_at_reg(self):
        #logging.debug ("Creating domain file at region"+ self.lon1.__str__()+"-"+self.lat2.__str__()+" "+self.lat1.__str__()+"-"+self.lat2.__str__())
        logging.info("Creating domain file at region:"+ self.tag)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fdomain_in, "xc", "yc", "ni", "nj")
        lat = f2["lat"]
        lon = f2["lon"]
        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3 = f2.isel(nj=yind, ni=xind)

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fdomain_in

        wfile = self.fdomain_out
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f2.close()
        f3.close()

    def create_surfdata_at_reg(self):
        #logging.debug ("Creating surface dataset file at region"+ self.lon1.__str__()+"-"+self.lat2.__str__()+" "+self.lat1.__str__()+"-"+self.lat2.__str__())
        logging.info("Creating surface dataset file at region:"+ self.tag)
        # create 1d coordinate variables to enable sel() method
        filename = self.fsurf_in
        f2 = self.create_1d_coord(filename, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f2["lat"]
        lon = f2["lon"]
        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3 = f2.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fsurf_in

        # mode 'w' overwrites file
        f3.to_netcdf(path=self.fsurf_out, mode="w")
        logging.info("created file (fsurf_out)" + self.fsurf_out)
        # f1.close();
        f2.close()
        f3.close()

    def create_landuse_at_reg(self):
        #logging.debug ("Creating landuse file at region"+ self.lon1.__str__()+"-"+self.lat2.__str__()+" "+self.lat1.__str__()+"-"+self.lat2.__str__())
        logging.info("Creating landuse file at region:"+ self.tag)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f2["lat"]
        lon = f2["lon"]
        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3 = f2.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fluse_in

        wfile = self.fluse_out
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f2.close()
        f3.close()
