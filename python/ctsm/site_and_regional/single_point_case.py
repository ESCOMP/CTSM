import os

import numpy as np
import xarray as xr

from datetime import date

from ctsm.site_and_regional.base_case import BaseCase, USRDAT_DIR


class SinglePointCase(BaseCase):
    """
    A case to encapsulate single point cases.
    ...
    Attributes
    ----------
    plat : float
        latitude
    plon : float
        longitude
    site_name: str -- default = None
        Site name
    overwrite_single_pft : bool
        flag to overwrite surface data with one uniform plant functional type
    dominant_pft: int
        index of plant functional type to set to 100% cover if overwrite_single_pft = True
    zero_nonveg_landunits : bool
        flag to set surface data to all natural vegetation (100% NATVEG, 0% other)
    uniform_snowpack
        flag to set the the surface data STD_ELEV to 0.0
    no_saturation_excess : bool
        flag to set the surface data FMAX to 0.0
    output_dir : str
        main output directory to write subset files to
    tag : str
        ending tag for output file naming
    fdomain_in : str
        file name of input domain file to subset
    fdomain_out : str
        file name of output subset domain domain file
    fluse_in : str
        file name of input land use file to subset
    fluse_out : str
        file name of output subset land use file
    fsurf_in : str
        file name of input surface data file to subset
    fsurf_out : str
        file name of output subset surface data file
    fdatmdomain_in : str
        file name of input DATM domain file to subset
    fdatmdomain_out : str
        file name of output subset DATM domain file
    datm_syr : int
        starting year for subset DATM data
    datm_eyr : int
        ending year for subset DATM data
    dir_tpqw : str
        input directory for TPQW DATM data
    dir_prec : str
        input directory for precipitation DATM data
    dir_solar : str
        input directory for solar DATM data
    tag_tpqw : str
        tag (file naming convention) for input TPQW DATM data
    tag_prec : str
        tag (file naming convention) for input precipitation DATM data
    tag_solar : str
        tag (file naming convention) for input solar DATM data
    name_tpqw : str
        stream name for TPQW DATM data
    name_prec : str
        stream name for precipitation DATM data
    name_solar : str
        stream name for solar DATM data
    dir_output_datm : str
        directory to write subset DATM data to (default to within main output directory)
    datm_stream_file : str
        file name of usr_nl_datm_streams file to write to for user_mods creation

    Methods
    -------
    create_tag:
        create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.
    create_fileout_name:
        creates a file name from a basename and a specified tag
    create_domain_at_point:
        Create domain file at a single point.
    create_landuse_at_point:
        Create landuse file at a single point.
    create_surfdata_at_point:
        Create surface dataset at a single point.
    create_datmdomain_at_point:
        Create DATM domain file at a single point.
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
            no_saturation_excess,
            output_dir,
    ):
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm, create_user_mods)
        self.plat = plat
        self.plon = plon
        self.site_name = site_name
        self.overwrite_single_pft = overwrite_single_pft
        self.dominant_pft = dominant_pft
        self.zero_nonveg_landunits = zero_nonveg_landunits
        self.uniform_snowpack = uniform_snowpack
        self.no_saturation_excess = no_saturation_excess
        self.output_dir = output_dir
        self.tag = None
        self.fdomain_in = None
        self.fdomain_out = None
        self.fluse_in = None
        self.fluse_out = None
        self.fsurf_in = None
        self.fsurf_out = None
        self.fdatmdomain_in = None
        self.fdatmdomain_out = None
        self.datm_syr = None
        self.datm_eyr = None
        self.name_tpqw = None
        self.name_prec = None
        self.name_solar = None
        self.dir_output_datm = None
        self.dir_tpqw = None
        self.tag_tpqw = None
        self.dir_prec = None
        self.tag_prec = None
        self.dir_input_datm = None
        self.tag_solar = None
        self.dir_solar = None
        self.datm_streams_file = None

    def create_tag(self):
        if self.site_name:
            self.tag = self.site_name
        else:
            self.tag = "{}_{}".format(str(self.plon), str(self.plat))

    @staticmethod
    def create_fileout_name(filename, tag):

        basename = os.path.basename(filename)
        items = basename.split("_")
        today = date.today()
        today_string = today.strftime("%y%m%d")
        new_string = "{}_{}_c{}.nc".format("_".join([items[0]] + items[2:7]), tag, today_string)

        return new_string

    def create_domain_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating domain file at ", self.plon, self.plat)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fdomain_in, "xc", "yc", "ni", "nj")
        # extract gridcell closest to plon/plat
        f3 = f2.sel(ni=self.plon, nj=self.plat, method="nearest")
        # expand dimensions
        f3 = f3.expand_dims(["nj", "ni"])

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fdomain_in

        wfile = os.path.join(self.output_dir, self.fdomain_out)
        f3.to_netcdf(path=wfile, mode="w")
        print("Successfully created file (fdomain_out) " + self.fdomain_out)
        f2.close()
        f3.close()

    def create_landuse_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating landuse file at ", self.plon, self.plat)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        # extract gridcell closest to plon/plat
        f3 = f2.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")

        # expand dimensions
        f3 = f3.expand_dims(["lsmlat", "lsmlon"])
        # specify dimension order
        # f3 = f3.transpose('time','lat','lon')
        f3 = f3.transpose(u"time", u"cft", u"natpft", u"lsmlat", u"lsmlon")
        # f3['YEAR'] = f3['YEAR'].squeeze()

        # revert expand dimensions of YEAR
        year = np.squeeze(np.asarray(f3["YEAR"]))
        x = xr.DataArray(year, coords={"time": f3["time"]}, dims="time", name="YEAR")
        x.attrs["units"] = "unitless"
        x.attrs["long_name"] = "Year of PFT data"
        f3["YEAR"] = x

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fluse_in

        wfile = os.path.join(self.output_dir, self.fluse_out)
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode="w")
        print("Successfully created file (fluse_out) " + self.fluse_out)
        f2.close()
        f3.close()

    def create_surfdata_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating surface dataset file at ", self.plon, self.plat)
        # create 1d coordinate variables to enable sel() method
        filename = os.path.join(self.output_dir, self.fsurf_in)
        f2 = self.create_1d_coord(filename, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        # extract gridcell closest to plon/plat
        f3 = f2.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")
        # expand dimensions
        f3 = f3.expand_dims(["lsmlat", "lsmlon"]).copy(deep=True)

        # update the plon and plat to match the surface data
        # we do this so that if we create user_mods the PTS_LON and PTS_LAT in CIME match
        # the surface data coordinates - which is required
        self.plat = f3.coords["lsmlat"].values[0]
        self.plon = f3.coords["lsmlon"].values[0]

        # modify surface data properties
        if self.overwrite_single_pft:
            f3["PCT_NAT_PFT"][:, :, :] = 0
            f3["PCT_NAT_PFT"][:, :, self.dominant_pft] = 100
        if self.zero_nonveg_landunits:
            f3["PCT_NATVEG"][:, :] = 100
            f3["PCT_CROP"][:, :] = 0
            f3["PCT_LAKE"][:, :] = 0.0
            f3["PCT_WETLAND"][:, :] = 0.0
            f3["PCT_URBAN"][:, :, ] = 0.0
            f3["PCT_GLACIER"][:, :] = 0.0
        if self.uniform_snowpack:
            f3["STD_ELEV"][:, :] = 20.0
        if self.no_saturation_excess:
            f3["FMAX"][:, :] = 0.0

        # specify dimension order
        # f3 = f3.transpose(u'time', u'cft', u'natpft', u'lsmlat', u'lsmlon')
        f3 = f3.transpose(
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
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fsurf_in
        del f3.attrs["History_Log"]
        # mode 'w' overwrites file
        f3.to_netcdf(path=self.fsurf_out, mode="w")
        print("Successfully created file (fsurf_out): " + self.fsurf_out)
        f2.close()
        f3.close()

    def create_datmdomain_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating DATM domain file at ", self.plon, self.plat)
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(self.fdatmdomain_in, "xc", "yc", "ni", "nj")
        # extract gridcell closest to plon/plat
        f3 = f2.sel(ni=self.plon, nj=self.plat, method="nearest")
        # expand dimensions
        f3 = f3.expand_dims(["nj", "ni"])
        wfile = self.fdatmdomain_out
        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fdatmdomain_in
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode="w")
        print("Successfully created file (fdatmdomain_out): " + self.fdatmdomain_out)
        f2.close()
        f3.close()

    def extract_datm_at(self, file_in, file_out):
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(file_in, "LONGXY", "LATIXY", "lon", "lat")
        # extract gridcell closest to plon/plat
        f3 = f2.sel(lon=self.plon, lat=self.plat, method="nearest")
        # expand dimensions
        f3 = f3.expand_dims(["lat", "lon"])
        # specify dimension order
        f3 = f3.transpose(u"scalar", "time", "lat", "lon")

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = file_in
        # mode 'w' overwrites file
        f3.to_netcdf(path=file_out, mode="w")
        print("Successfully created file: " + file_out)
        f2.close()
        f3.close()

    def write_shell_commands(self, file):
        """
        writes out xml commands commands to a file (i.e. shell_commands) for single-point runs

        file - file connection to shell_commands file
        """

        # write_to_file surrounds text with newlines
        self.write_to_file("! Change below line if you move the subset data directory", file)
        self.write_to_file("./xmlchange {}={}".format(USRDAT_DIR, self.output_dir), file)
        self.write_to_file("./xmlchange PTS_LON=".format(str(self.plon)), file)
        self.write_to_file("./xmlchange PTS_LAT=".format(str(self.plat)), file)
        self.write_to_file("./xmlchange MPILIB=mpi-serial", file)
        file.close()

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

    def create_datm_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating DATM files at ", self.plon, self.plat, ".")

        # --  create data files
        infile = []
        outfile = []
        solarfiles = []
        precfiles = []
        tpqwfiles = []
        for y in range(self.datm_syr, self.datm_eyr + 1):
            ystr = str(y)
            for m in range(1, 13):
                mstr = str(m)
                if m < 10:
                    mstr = "0" + mstr

                dtag = ystr + "-" + mstr

                fsolar = os.path.join(self.dir_input_datm, self.dir_solar, "{}{}.nc".format(self.tag_solar, dtag))
                fsolar2 = "{}{}.{}.nc".format(self.tag_solar, self.tag, dtag)
                fprecip = os.path.join(self.dir_input_datm, self.dir_solar, "{}{}.nc".format(self.tag_prec, dtag))
                fprecip2 = "{}{}.{}.nc".format(self.tag_prec, self.tag, dtag)
                ftpqw = os.path.join(self.dir_input_datm, self.dir_tpqw, "{}{}.nc".format(self.tag_tpqw, dtag))
                ftpqw2 = "{}{}.{}.nc".format(self.tag_tpqw, self.tag, dtag)

                outdir = os.path.join(self.output_dir, self.dir_output_datm)
                infile += [fsolar, fprecip, ftpqw]
                outfile += [os.path.join(outdir, fsolar2),
                            os.path.join(outdir, fprecip2),
                            os.path.join(outdir, ftpqw2)]
                solarfiles.append(os.path.join("${}".format(USRDAT_DIR), self.dir_output_datm, fsolar2))
                precfiles.append(os.path.join("${}".format(USRDAT_DIR), self.dir_output_datm, fprecip2))
                tpqwfiles.append(os.path.join("${}".format(USRDAT_DIR), self.dir_output_datm, ftpqw2))

        nm = len(infile)
        for n in range(nm):
            print(outfile[n])
            file_in = infile[n]
            file_out = outfile[n]
            self.extract_datm_at(file_in, file_out)

        print("All DATM files are created in: " + self.output_dir)

        # write to user_nl_datm_streams if specified
        if self.create_user_mods:
            with open(self.datm_streams_file, "a") as file:
                self.write_datm_streams_lines(self.name_solar, solarfiles, file)
                self.write_datm_streams_lines(self.name_prec, precfiles, file)
                self.write_datm_streams_lines(self.name_tpqw, tpqwfiles, file)
