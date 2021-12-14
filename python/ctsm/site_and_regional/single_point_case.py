"""
This module includes the definition for SinglePointCase class.
"""

#-- Import libraries
#-- Import Python Standard Libraries
import os
import logging

#-- 3rd party libraries
import numpy as np
import xarray as xr

#-- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase

logger = logging.getLogger(__name__)

class SinglePointCase(BaseCase):
    """
    A class to encapsulate everything for single point cases.

    ...

    Attributes
    ----------
    plat : float
        latitude of the single point
    plon : float
        longitude of the single point
    site_name: str -- default = None
        Site name
    create_domain : bool
        flag for creating domain file
    create_surfdata : bool
        flag for creating surface dataset
    create_landuse : bool
        flag for creating landuse file
    create_datm : bool
        flag for creating DATM files
    overwrite_single_pft : bool
        flag to overwrite the whole grid 100% single PFT.
    dominant_pft : int
        dominant pft type for this single point
    zero_nonveg_landunits : bool
        flag for setting all non-vegetation landunits to zero
    overwrite_single_pft : bool
        flag for creating datasets using uniform snowpack
    saturation_excess : bool
        flag for making dataset using saturation excess

    Methods
    -------
    create_tag
        create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.

    create_domain_at_point
        Create domain file at a single point.

    create_landuse_at_point:
        Create landuse file at a single point.

    create_surfdata_at_point:
        Create surface dataset at a single point.

    create_datmdomain_at_point:
        Create DATM domain file at a single point.

    extract_datm_at:
        Extract DATM for one file at a single point.

    create_datm_at_point:
        Extract all DATM data at a single point.
    """

    def __init__(
        self,
        plat,
        plon,
        site_name,
        create_domain,
        create_surfdata,
        create_landuse,
        create_datm,
        overwrite_single_pft,
        dominant_pft,
        zero_nonveg_landunits,
        uniform_snowpack,
        saturation_excess,
    ):
        """ 
        Initializes SinglePointCase with the given arguments.

        """

        super().__init__(create_domain, create_surfdata, create_landuse, create_datm)
        self.plat = plat
        self.plon = plon
        self.site_name = site_name
        self.overwrite_single_pft = overwrite_single_pft
        self.dominant_pft = dominant_pft
        self.zero_nonveg_landunits = zero_nonveg_landunits
        self.uniform_snowpack = uniform_snowpack
        self.saturation_excess = saturation_excess

    def create_tag(self):
        """
        Create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.
        """
        if self.site_name:
            self.tag = self.site_name
        else:
            self.tag = str(self.plon) + "_" + str(self.plat)


    def create_domain_at_point(self):
        """
        Create domain file for this SinglePointCase class.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating domain file at "+ self.plon.__str__()+" "+ self.plat.__str__()+".")

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(self.fdomain_in, "xc", "yc", "ni", "nj")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(ni=self.plon, nj=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["nj", "ni"])

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fdomain_in

        wfile = self.fdomain_out
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdomain_out)" + self.fdomain_out)
        f_in.close()
        f_out.close()

    def create_landuse_at_point(self):
        """
        Create landuse file at a single point.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating landuse file at "+ self.plon.__str__()+" "+ self.plat.__str__()+".")

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(self.fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["lsmlat", "lsmlon"])

        # specify dimension order
        # f_out = f_out.transpose('time','lat','lon')
        f_out = f_out.transpose(u"time", u"cft", u"natpft", u"lsmlat", u"lsmlon")
        # f_out['YEAR'] = f_out['YEAR'].squeeze()

        # revert expand dimensions of YEAR
        year = np.squeeze(np.asarray(f_out["YEAR"]))
        x = xr.DataArray(year, coords={"time": f_out["time"]}, dims="time", name="YEAR")
        x.attrs["units"] = "unitless"
        x.attrs["long_name"] = "Year of PFT data"
        f_out["YEAR"] = x

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fluse_in

        wfile = self.fluse_out
        # mode 'w' overwrites file
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (luse_out)" + self.fluse_out+ ".")
        f_in.close()
        f_out.close()

    def create_surfdata_at_point(self):
        """
        Create surface data file at a single point.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating surface dataset file at "+ self.plon.__str__()+" "+ self.plat.__str__()+".")

        # create 1d coordinate variables to enable sel() method
        filename = self.fsurf_in
        f_in = self.create_1d_coord(filename, "LONGXY", "LATIXY", "lsmlon", "lsmlat")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["lsmlat", "lsmlon"]).copy(deep=True)

        # modify surface data properties
        if self.overwrite_single_pft:
            f_out["PCT_NAT_PFT"][:, :, :] = 0
            f_out["PCT_NAT_PFT"][:, :, self.dominant_pft] = 100
        if self.zero_nonveg_landunits:
            f_out["PCT_NATVEG"][:, :] = 100
            f_out["PCT_CROP"][:, :] = 0
            f_out["PCT_LAKE"][:, :] = 0.0
            f_out["PCT_WETLAND"][:, :] = 0.0
            f_out["PCT_URBAN"][
                :,
                :,
            ] = 0.0
            f_out["PCT_GLACIER"][:, :] = 0.0
        if self.uniform_snowpack:
            f_out["STD_ELEV"][:, :] = 20.0
        if not self.saturation_excess:
            f_out["FMAX"][:, :] = 0.0

        # specify dimension order
        # f_out = f_out.transpose(u'time', u'cft', u'natpft', u'lsmlat', u'lsmlon')
        f_out = f_out.transpose(
            u"time",
            u"cft",
            u"lsmpft",
            u"natpft",
            u"nglcec",
            u"nglcecp1",
            u"nlevsoi",
            u"nlevurb",
            u"numrad",
            u"numurbl",
            "lsmlat",
            "lsmlon",
        )

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fsurf_in
        del f_out.attrs["History_Log"]
        # mode 'w' overwrites file
        f_out.to_netcdf(path=self.fsurf_out, mode="w")
        logging.info("Successfully created file (fsurf_out) :" + self.fsurf_out)
        f_in.close()
        f_out.close()

    def create_datmdomain_at_point(self):
        """
        Create DATM domain file at a single point
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating DATM domain file at "+ self.plon.__str__()+" "+ self.plat.__str__()+".")

        # create 1d coordinate variables to enable sel() method
        filename = self.fdatmdomain_in
        f_in = self.create_1d_coord(filename, "xc", "yc", "ni", "nj")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(ni=self.plon, nj=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["nj", "ni"])
        wfile = self.fdatmdomain_out

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = self.fdatmdomain_in

        # mode 'w' overwrites file
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdatmdomain_out) :" + self.fdatmdomain_out)
        f_in.close()
        f_out.close()

    def extract_datm_at(self, file_in, file_out):
        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(file_in, "LONGXY", "LATIXY", "lon", "lat")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(lon=self.plon, lat=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["lat", "lon"])

        # specify dimension order
        f_out = f_out.transpose(u"scalar", "time", "lat", "lon")

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = file_in

        # mode 'w' overwrites file
        f_out.to_netcdf(path=file_out, mode="w")
        logging.info("Successfully created file :" + file_out)
        f_in.close()
        f_out.close()

    def create_datm_at_point(self):
        """
        Create all DATM dataset at a point.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating DATM files at "+ self.plon.__str__()+" "+ self.plat.__str__()+".")
        # --  specify subdirectory names and filename prefixes
        solrdir = "Solar/"
        precdir = "Precip/"
        tpqwldir = "TPHWL/"
        prectag = "clmforc.GSWP3.c2011.0.5x0.5.Prec."
        solrtag = "clmforc.GSWP3.c2011.0.5x0.5.Solr."
        tpqwtag = "clmforc.GSWP3.c2011.0.5x0.5.TPQWL."

        # --  create data files
        infile = []
        outfile = []
        for y in range(self.datm_syr, self.datm_eyr + 1):
            ystr = str(y)
            for m in range(1, 13):
                mstr = str(m)
                if m < 10:
                    mstr = "0" + mstr

                dtag = ystr + "-" + mstr

                fsolar = self.dir_input_datm + solrdir + solrtag + dtag + ".nc"
                fsolar2 = self.dir_output_datm + solrtag + self.tag + "." + dtag + ".nc"
                fprecip = self.dir_input_datm + precdir + prectag + dtag + ".nc"
                fprecip2 = (
                    self.dir_output_datm + prectag + self.tag + "." + dtag + ".nc"
                )
                ftpqw = self.dir_input_datm + tpqwldir + tpqwtag + dtag + ".nc"
                ftpqw2 = self.dir_output_datm + tpqwtag + self.tag + "." + dtag + ".nc"

                infile += [fsolar, fprecip, ftpqw]
                outfile += [fsolar2, fprecip2, ftpqw2]

        nm = len(infile)
        for n in range(nm):
            logging.debug(outfile[n])
            file_in = infile[n]
            file_out = outfile[n]
            self.extract_datm_at(file_in, file_out)

        logging.info("All DATM files are created in: "+ self.dir_output_datm+".")
