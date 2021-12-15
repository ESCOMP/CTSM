"""
This module includes the definition for SinglePointCase class.
"""

# -- Import libraries
# -- Import Python Standard Libraries
import logging
import os

# -- 3rd party libraries
import numpy as np
import xarray as xr

# -- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase, USRDAT_DIR, DatmFiles

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
    create_user_mods : bool
        flag for creating user mods directories and files
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
            create_user_mods,
            overwrite_single_pft,
            dominant_pft,
            zero_nonveg_landunits,
            uniform_snowpack,
            saturation_excess,
            output_dir,
    ):
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm,
                         create_user_mods)
        self.plat = plat
        self.plon = plon
        self.site_name = site_name
        self.overwrite_single_pft = overwrite_single_pft
        self.dominant_pft = dominant_pft
        self.zero_nonveg_landunits = zero_nonveg_landunits
        self.uniform_snowpack = uniform_snowpack
        self.saturation_excess = saturation_excess
        self.output_dir = output_dir
        self.tag = None

    def create_tag(self):
        """
        Create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.
        """
        if self.site_name:
            self.tag = self.site_name
        else:
            self.tag = "{}_{}".format(str(self.plon), str(self.plat))

    def create_domain_at_point(self, indir, file):
        """
        Create domain file for this SinglePointCase class.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating domain file at %s, %s.", self.plon.__str__(), self.plat.__str__())

        # specify files
        fdomain_in = os.path.join(indir, file)
        fdomain_out = self.add_tag_to_filename(fdomain_in, self.tag)
        logging.info("fdomain_in:  %s", fdomain_in)
        logging.info("fdomain_out: %s", os.path.join(self.output_dir, fdomain_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fdomain_in, "xc", "yc", "ni", "nj")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(ni=self.plon, nj=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["nj", "ni"])

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fdomain_in

        wfile = os.path.join(self.output_dir, fdomain_out)
        f_out.to_netcdf(path=fdomain_out, mode="w")
        logging.info("Successfully created file (fdomain_out) %s", wfile)
        f_in.close()
        f_out.close()

    def create_landuse_at_point(self, indir, file, user_mods_dir):
        """
        Create landuse file at a single point.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating land use file at %s, %s.", self.plon.__str__(), self.plat.__str__())

        # specify files
        fluse_in = os.path.join(indir, file)
        fluse_out = self.add_tag_to_filename(fluse_in, self.tag)
        logging.info("fluse_in:  %s", fluse_in)
        logging.info("fluse_out: %s", os.path.join(self.output_dir, fluse_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(
            fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat"
        )

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
        temp_xr = xr.DataArray(year, coords={"time": f_out["time"]}, dims="time", name="YEAR")
        temp_xr.attrs["units"] = "unitless"
        temp_xr.attrs["long_name"] = "Year of PFT data"
        f_out["YEAR"] = temp_xr

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fluse_in

        wfile = os.path.join(self.output_dir, fluse_out)
        # mode 'w' overwrites file
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fluse_out), %s", wfile)
        f_in.close()
        f_out.close()

        # write to user_nl_clm data if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "landuse = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                self.write_to_file(line, nl_clm)

    def create_surfdata_at_point(self, indir, file, user_mods_dir):
        """
        Create surface data file at a single point.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info(
            "Creating surface dataset file at %s, %s", self.plon.__str__(), self.plat.__str__())

        # specify file
        fsurf_in = os.path.join(indir, file)
        fsurf_out = self.add_tag_to_filename(fsurf_in, self.tag)
        logging.info("fsurf_in:  %s", fsurf_in)
        logging.info("fsurf_out: %s", os.path.join(self.output_dir, fsurf_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fsurf_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["lsmlat", "lsmlon"]).copy(deep=True)

        # update the plon and plat to match the surface data
        # we do this so that if we create user_mods the PTS_LON and PTS_LAT in CIME match
        # the surface data coordinates - which is required
        self.plat = f_out.coords["lsmlat"].values[0]
        self.plon = f_out.coords["lsmlon"].values[0]

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
        f_out.attrs["Created_from"] = fsurf_in
        del f_out.attrs["History_Log"]
        # mode 'w' overwrites file
        wfile = os.path.join(self.output_dir, fsurf_out)
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fsurf_out) %s", wfile)
        f_in.close()
        f_out.close()

        # write to user_nl_clm if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "fsurdat = '${}'".format(os.path.join(USRDAT_DIR, fsurf_out))
                self.write_to_file(line, nl_clm)

    def create_datmdomain_at_point(self, datm_tuple: DatmFiles):
        """
        Create DATM domain file at a single point
        """
        logging.info("----------------------------------------------------------------------")
        logging.info(
            "Creating DATM domain file at %s, %s", self.plon.__str__(), self.plat.__str__())

        # specify files
        fdatmdomain_in = os.path.join(datm_tuple.indir, datm_tuple.fdomain_in)
        datm_file = self.add_tag_to_filename(fdatmdomain_in, self.tag)
        fdatmdomain_out = os.path.join(datm_tuple.outdir, datm_file)
        logging.info("fdatmdomain_in:  %s", fdatmdomain_in)
        logging.info("fdatmdomain out: %s", os.path.join(self.output_dir, fdatmdomain_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fdatmdomain_in, "xc", "yc", "ni", "nj")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(ni=self.plon, nj=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["nj", "ni"])

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fdatmdomain_in

        # mode 'w' overwrites file
        wfile = os.path.join(self.output_dir, fdatmdomain_out)
        f_out.to_netcdf(path=wfile, mode="w")
        logging.info("Successfully created file (fdatmdomain_out) : %s", wfile)
        f_in.close()
        f_out.close()

    def extract_datm_at(self, file_in, file_out):
        """
        Create a DATM dataset at a point.
        """
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
        logging.info("Successfully created file : %s", file_out)
        f_in.close()
        f_out.close()

    def write_shell_commands(self, file):
        """
        writes out xml commands commands to a file (i.e. shell_commands) for single-point runs
        """
        # write_to_file surrounds text with newlines
        with open(file, 'w') as nl_file:
            self.write_to_file("# Change below line if you move the subset data directory", nl_file)
            self.write_to_file("./xmlchange {}={}".format(USRDAT_DIR, self.output_dir), nl_file)
            self.write_to_file("./xmlchange PTS_LON={}".format(str(self.plon)), nl_file)
            self.write_to_file("./xmlchange PTS_LAT={}".format(str(self.plat)), nl_file)
            self.write_to_file("./xmlchange MPILIB=mpi-serial", nl_file)

    def write_datm_streams_lines(self, streamname, datmfiles, file):
        """
        writes out lines for the user_nl_datm_streams file for a specific DATM stream
        for using subset DATM data at a single point

        streamname - stream name (e.g. TPQW)
        datmfiles - comma-separated list (str) of DATM file names
        file - file connection to user_nl_datm_streams file
        """
        self.write_to_file("{}:datafiles={}".format(streamname, ','.join(datmfiles)), file)
        self.write_to_file("{}:mapalgo=none".format(streamname), file)
        self.write_to_file("{}:meshfile=none".format(streamname), file)

    def create_datm_at_point(self, datm_tuple: DatmFiles, datm_syr, datm_eyr, datm_streams_file):
        """
        Create all of a DATM dataset at a point.
        """
        logging.info("----------------------------------------------------------------------")
        logging.info("Creating DATM files at %s, %s", self.plon.__str__(), self.plat.__str__())

        # --  create data files
        infile = []
        outfile = []
        solarfiles = []
        precfiles = []
        tpqwfiles = []
        for year in range(datm_syr, datm_eyr + 1):
            ystr = str(year)
            for month in range(1, 13):
                mstr = str(month)
                if month < 10:
                    mstr = "0" + mstr

                dtag = ystr + "-" + mstr

                fsolar = os.path.join(datm_tuple.indir, datm_tuple.dir_solar,
                                      "{}{}.nc".format(datm_tuple.tag_solar, dtag))
                fsolar2 = "{}{}.{}.nc".format(datm_tuple.tag_solar, self.tag, dtag)
                fprecip = os.path.join(datm_tuple.indir, datm_tuple.dir_prec,
                                       "{}{}.nc".format(datm_tuple.tag_prec, dtag))
                fprecip2 = "{}{}.{}.nc".format(datm_tuple.tag_prec, self.tag, dtag)
                ftpqw = os.path.join(datm_tuple.indir, datm_tuple.dir_tpqw,
                                     "{}{}.nc".format(datm_tuple.tag_tpqw, dtag))
                ftpqw2 = "{}{}.{}.nc".format(datm_tuple.tag_tpqw, self.tag, dtag)

                outdir = os.path.join(self.output_dir, datm_tuple.outdir)
                infile += [fsolar, fprecip, ftpqw]
                outfile += [os.path.join(outdir, fsolar2),
                            os.path.join(outdir, fprecip2),
                            os.path.join(outdir, ftpqw2)]
                solarfiles.append(
                    os.path.join("${}".format(USRDAT_DIR), datm_tuple.outdir, fsolar2))
                precfiles.append(
                    os.path.join("${}".format(USRDAT_DIR), datm_tuple.outdir, fprecip2))
                tpqwfiles.append(
                    os.path.join("${}".format(USRDAT_DIR), datm_tuple.outdir, ftpqw2))

        nm = len(infile)
        for n in range(nm):
            logging.debug(outfile[n])
            file_in = infile[n]
            file_out = outfile[n]
            self.extract_datm_at(file_in, file_out)

        logging.info("All DATM files are created in: %s", datm_tuple.outdir)

        # write to user_nl_datm_streams if specified
        if self.create_user_mods:
            with open(datm_streams_file, "a") as file:
                self.write_datm_streams_lines(datm_tuple.name_solar, solarfiles, file)
                self.write_datm_streams_lines(datm_tuple.name_prec, precfiles, file)
                self.write_datm_streams_lines(datm_tuple.name_tpqw, tpqwfiles, file)
