"""
Holds the class RegionalCase
"""
import logging
import os

import numpy as np

from ctsm.site_and_regional.base_case import BaseCase, USRDAT_DIR

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
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm,
                         create_user_mods)
        self.lat1 = lat1
        self.lat2 = lat2
        self.lon1 = lon1
        self.lon2 = lon2
        self.reg_name = reg_name
        self.output_dir = output_dir
        self.tag = None

    def create_tag(self):
        if self.reg_name:
            self.tag = self.reg_name
        else:
            self.tag = "{}-{}_{}-{}".format(str(self.lon1), str(self.lon2), str(self.lat1),
                                            str(self.lat2))

    def create_domain_at_reg(self, indir, file):

        # specify files
        fdomain_in = os.path.join(indir, file)
        fdomain_out = os.path.join(self.output_dir, "domain.lnd.fv1.9x2.5_gx1v7." +
                                          self.tag + "_170518.nc")
        logging.info("fdomain_in:  %s", fdomain_in)
        logging.info("fdomain_out: %s", fdomain_out)
        logging.info("Creating domain file at region: %s", self.tag)

        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(fdomain_in, "xc", "yc", "ni", "nj")
        lat = f2["lat"]
        lon = f2["lon"]
        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3 = f2.isel(nj=yind, ni=xind)

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = fdomain_in

        # mode 'w' overwrites file
        f3.to_netcdf(path=fdomain_out, mode="w")
        logging.info("Successfully created file (fdomain_out) at %s", fdomain_out)
        f2.close()
        f3.close()

    def create_surfdata_at_reg(self, indir, file, user_mods_dir):

        logging.info("Creating surface dataset file at region: %s", self.tag)

        # specify files
        fsurf_in = os.path.join(indir, file)
        fsurf_out = os.path.join(self.output_dir,
                                        "surfdata_1.9x2.5_78pfts_CMIP6_simyr1850_" + self.tag
                                        + "_c170824.nc")
        logging.info("fsurf_in:  %s", fsurf_in)
        logging.info("fsurf_out: %s", fsurf_out)

        # create 1d coordinate variables to enable sel() method
        filename = fsurf_in
        f2 = self.create_1d_coord(filename, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f2["lat"]
        lon = f2["lon"]
        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3 = f2.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = fsurf_in

        # mode 'w' overwrites file
        f3.to_netcdf(path=fsurf_out, mode="w")
        logging.info("created file (fsurf_out) %s", fsurf_out)
        # f1.close();
        f2.close()
        f3.close()

        # write to user_nl_clm if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "fsurdat = '${}'".format(os.path.join(USRDAT_DIR, fsurf_out))
                self.write_to_file(line, nl_clm)

    def create_landuse_at_reg(self, indir, file, user_mods_dir):
        logging.info("Creating landuse file at region: %s", self.tag)

        # specify files
        fluse_in = os.path.join(indir, file)
        fluse_out = os.path.join(self.output_dir,
                                        "landuse.timeseries_1.9x2"
                                        ".5_hist_78pfts_CMIP6_simyr1850-2015_" +
                                        self.tag + ".c170824.nc")
        logging.info("fluse_in:  %s", fluse_in)
        logging.info("fluse_out: %s", fluse_out)

        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        lat = f2["lat"]
        lon = f2["lon"]
        # subset longitude and latitude arrays
        xind = np.where((lon >= self.lon1) & (lon <= self.lon2))[0]
        yind = np.where((lat >= self.lat1) & (lat <= self.lat2))[0]
        f3 = f2.isel(lsmlat=yind, lsmlon=xind)

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = fluse_in

        # mode 'w' overwrites file
        f3.to_netcdf(path=fluse_out, mode="w")
        logging.info("Successfully created file (fluse_out) %s", fluse_out)
        f2.close()
        f3.close()

        # write to user_nl_clm data if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "landuse = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                self.write_to_file(line, nl_clm)
