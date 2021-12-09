from ctsm.site_and_regional.base_case import BaseCase
import os
import numpy as np
import xarray as xr
from datetime import date


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
        no_saturation_excess,
        output_dir
    ):
        super().__init__(create_domain, create_surfdata, create_landuse, create_datm)
        self.plat = plat
        self.plon = plon
        self.site_name = site_name
        self.overwrite_single_pft = overwrite_single_pft
        self.dominant_pft = dominant_pft
        self.zero_nonveg_landunits = zero_nonveg_landunits
        self.uniform_snowpack = uniform_snowpack
        self.no_saturation_excess = no_saturation_excess
        self.output_dir = output_dir

    def create_tag(self):
        if self.site_name:
            self.tag = self.site_name
        else:
            self.tag = str(self.plon) + "_" + str(self.plat)

    @staticmethod
    def create_fileout_name(filename, tag):

        basename = os.path.basename(filename)
        items = basename.split("_")
        today = date.today()
        today_string = today.strftime("%y%m%d")
        new_string = (
            items[0]
            + "_"
            + items[2]
            + "_"
            + items[3]
            + "_"
            + items[4]
            + "_"
            + items[5]
            + "_"
            + items[6]
            + "_"
            + tag
            + "_c"
            + today_string
            + ".nc"
        )
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
        print("Successfully created file (fdomain_out)" + self.fdomain_out)
        f2.close()
        f3.close()

    def create_landuse_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating landuse file at ", self.plon, self.plat, ".")
        # create 1d coordinate variables to enable sel() method
        f2 = self.create_1d_coord(
            self.fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
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
        x = xr.DataArray(
            year, coords={"time": f3["time"]}, dims="time", name="YEAR")
        x.attrs["units"] = "unitless"
        x.attrs["long_name"] = "Year of PFT data"
        f3["YEAR"] = x

        # update attributes
        self.update_metadata(f3)
        f3.attrs["Created_from"] = self.fluse_in

        wfile = os.path.join(self.output_dir, self.fluse_out)
        # mode 'w' overwrites file
        f3.to_netcdf(path=wfile, mode="w")
        print("Successfully created file (luse_out)" + self.fluse_out, ".")
        f2.close()
        f3.close()

    def create_surfdata_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating surface dataset file at ", self.plon, self.plat, ".")
        # create 1d coordinate variables to enable sel() method
        filename = os.path.join(self.output_dir, self.fsurf_in)
        f2 = self.create_1d_coord(
            filename, "LONGXY", "LATIXY", "lsmlon", "lsmlat")
        # extract gridcell closest to plon/plat
        f3 = f2.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")
        # expand dimensions
        f3 = f3.expand_dims(["lsmlat", "lsmlon"]).copy(deep=True)

        # update the plon and plat to match the surface data
        self.plat = f3.coords['lsmlat'].values[0]
        self.plon = f3.coords['lsmlon'].values[0]

        # modify surface data properties
        if self.overwrite_single_pft:
            f3["PCT_NAT_PFT"][:, :, :] = 0
            f3["PCT_NAT_PFT"][:, :, self.dominant_pft] = 100
        if self.zero_nonveg_landunits:
            f3["PCT_NATVEG"][:, :] = 100
            f3["PCT_CROP"][:, :] = 0
            f3["PCT_LAKE"][:, :] = 0.0
            f3["PCT_WETLAND"][:, :] = 0.0
            f3["PCT_URBAN"][
                :,
                :,
            ] = 0.0
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
        print("Successfully created file (fsurf_out) :" + self.fsurf_out)
        f2.close()
        f3.close()

    def create_datmdomain_at_point(self):
        print("----------------------------------------------------------------------")
        print("Creating DATM domain file at ", self.plon, self.plat, ".")
        # create 1d coordinate variables to enable sel() method
        filename = self.fdatmdomain_in
        f2 = self.create_1d_coord(filename, "xc", "yc", "ni", "nj")
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
        print("Successfully created file (fdatmdomain_out) :" + self.fdatmdomain_out)
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
        print("Successfully created file :" + file_out)
        f2.close()
        f3.close()

    def write_shell_commands(self, file):
        # writes out shell commands for single-point runs

        file.write('! Change below line if you move the subset data directory')
        file.write('\n' + './xmlchange CLM_USRDAT_DIR=' + self.out_dir + '\n')
        file.write('\n' + "./xmlchange PTS_LON=" + str(self.plon) + '\n')
        file.write('\n' + "./xmlchange PTS_LAT=" + str(self.plat) + '\n')
        file.write('\n' + "./xmlchange MPILIB=mpi-serial" + '\n')
        file.close()

    def write_nl_commands(self, streamname, file):
        line_mapalgo = streamname + ':mapalgo=none'
        line_meshfile = streamname + ':meshfile=none'

        file.write("\n" + line_meshfile + "\n")
        file.write("\n" + line_mapalgo + "\n")

    def create_datm_at_point(self, create_user_mods, datm_streams_file):
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

                fsolar = os.path.join(self.dir_solar, self.tag_solar + dtag + ".nc")
                fsolar2 = os.path.join(self.tag_solar + self.tag + "." + dtag + ".nc")
                fprecip = os.path.join(self.dir_input_datm, self.dir_prec, self.tag_prec + dtag + ".nc")
                fprecip2 = os.path.join(self.tag_prec + self.tag + "." + dtag + ".nc")
                ftpqw = os.path.join(self.dir_input_datm, self.dir_tpqw, self.tag_tpqw + dtag + ".nc")
                ftpqw2 = os.path.join(self.tag_tpqw + self.tag + "." + dtag + ".nc")

                outdir = os.path.join(self.output_dir, self.dir_output_datm)
                infile += [fsolar, fprecip, ftpqw]
                outfile += [os.path.join(outdir, fsolar2), os.path.join(outdir, fprecip2),
                            os.path.join(outdir, ftpqw2)]
                solarfiles.append(os.path.join("$CLM_USRDAT_DIR", self.dir_output_dtam, fsolar2))
                precfiles.append(os.path.join("$CLM_USRDAT_DIR", self.dir_output_dtam, fprecip2))
                tpqwfiles.append(os.path.join("$CLM_USRDAT_DIR", self.dir_output_dtam, ftpqw2))

        nm = len(infile)
        for n in range(nm):
            print(outfile[n])
            file_in = infile[n]
            file_out = outfile[n]
            self.extract_datm_at(file_in, file_out)

        print("All DATM files are created in: " + outdir)

        # write to user_nl_datm_streams if specified
        if create_user_mods:
            solarfile_line = self.name_solar + ':datafiles=' + ','.join(solarfiles)
            precfile_line = self.name_prec + ':datafiles=' + ','.join(precfiles)
            tpqwfile_line = self.name_tpqw + ':datafiles=' + ','.join(tpqwfiles)

            with open(datm_streams_file, 'a') as user_nl_file:
                user_nl_file.write('\n' + solarfile_line + '\n')
                self.write_nl_commands(self.name_solar, user_nl_file)
                user_nl_file.write('\n' + precfile_line + '\n')
                self.write_nl_commands(self.name_prec, user_nl_file)
                user_nl_file.write('\n' + tpqwfile_line + '\n')
                self.write_nl_commands(self.name_tpqw, user_nl_file)
