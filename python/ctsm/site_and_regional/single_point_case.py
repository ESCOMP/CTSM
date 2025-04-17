"""
This module includes the definition for SinglePointCase class.
"""

# -- Import libraries
# -- Import Python Standard Libraries
import logging
import os
import argparse

# -- 3rd party libraries
import numpy as np
import xarray as xr

# -- import local classes for this script
from ctsm.site_and_regional.base_case import BaseCase, USRDAT_DIR, DatmFiles
from ctsm.utils import add_tag_to_filename, ensure_iterable

logger = logging.getLogger(__name__)

NAT_PFT = 15  # natural pfts
NUM_PFT = 17  # for runs with generic crops
MAX_PFT = 78  # for runs with explicit crops

# -- constants to represent months of year
FIRST_MONTH = 1
LAST_MONTH = 12


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
    dom_pft : int
        dominant pft type for this single point (None if not specified)
    evenly_split_cropland : bool
        flag for splitting cropland evenly among all crop types
    pct_pft : list
        weight or percentage of each pft.
    cth : list
        canopy top height (m)
    cbh : list
        canopy bottom height (m)
    num_pft : list
        total number of pfts for surface dataset (if crop 78 pft, else 16 pft)
    uni_snow : bool
        flag for creating datasets using uniform snowpack
    saturation_excess : bool
        flag for making dataset using saturation excess
    overwrite : bool
        flag for over-writing files if they already exist

    Methods
    -------
    create_tag
        create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.

    create_domain_at_point
        Create domain file at a single point.

    create_landuse_at_point:
        Create landuse file at a single point.

    modify_surfdata_atpoint:
        Modify surface dataset based on combination of user choices.

    create_surfdata_at_point:
        Create surface dataset at a single point.

    create_datmdomain_at_point:
        Create DATM domain file at a single point.

    extract_datm_at:
        Extract DATM for one file at a single point.

    create_datm_at_point:
        Extract all DATM data at a single point.
    """

    # pylint: disable=too-many-instance-attributes
    # the ones we have are useful

    def __init__(
        self,
        plat,
        plon,
        *,
        site_name,
        create_domain,
        create_surfdata,
        create_landuse,
        create_datm,
        create_user_mods,
        dom_pft,
        evenly_split_cropland,
        pct_pft,
        num_pft,
        cth,
        cbh,
        include_nonveg,
        uni_snow,
        cap_saturation,
        out_dir,
        overwrite,
    ):
        super().__init__(
            create_domain=create_domain,
            create_surfdata=create_surfdata,
            create_landuse=create_landuse,
            create_datm=create_datm,
            create_user_mods=create_user_mods,
            overwrite=overwrite,
        )
        self.plat = plat
        self.plon = plon
        self.site_name = site_name
        self.dom_pft = dom_pft
        self.evenly_split_cropland = evenly_split_cropland
        self.pct_pft = pct_pft
        self.num_pft = num_pft
        self.cth = cth
        self.cbh = cbh
        self.include_nonveg = include_nonveg
        self.uni_snow = uni_snow
        self.cap_saturation = cap_saturation
        self.out_dir = out_dir

        self.create_tag()
        self.check_dom_pft()
        # self.check_nonveg()
        self.check_pct_pft()

    def create_tag(self):
        """
        Create a tag for single point which is the site name
        or the "lon-lat" format if the site name does not exist.
        """
        if self.site_name:
            self.tag = self.site_name
        else:
            self.tag = "{}_{}".format(str(self.plon), str(self.plat))

    def check_dom_pft(self):
        """
        A function to sanity check values in dom_pft:

        - Compare dom_pft (values if more than one) with num_pft:
          i.e. If dom_pft is 18 without crop it fails.

        - Check for mixed land-units:
          If we have more than one dom_pft, they should be in the
          same range.
          e.g. If users specified multiple dom_pft, they should be
          either in :
            - 0 - NAT_PFT-1 range
            or
            - NAT_PFT - MAX_PFT range
            - give an error: mixed land units not possible

        -------------
        Raises:
            Error (ArgumentTypeError):
                If any dom_pft is bigger than MAX_PFT.
            Error (ArgumentTypeError):
                If any dom_pft is less than 1.
            Error (ArgumentTypeError):
                If mixed land units are chosen.
                dom_pft values are both in range of (0 - NAT_PFT-1) and (NAT_PFT - MAX_PFT).


        """

        if self.dom_pft is None:
            logger.warning(
                "No dominant pft type is chosen. "
                "If you want to choose a dominant pft type, please use --dompft flag."
            )
        else:
            min_dom_pft = min(self.dom_pft)
            max_dom_pft = max(self.dom_pft)

            # -- check dom_pft values should be between 0-MAX_PFT
            if min_dom_pft < 0 or max_dom_pft > MAX_PFT:
                err_msg = "values for --dompft should be between 1 and 78."
                raise argparse.ArgumentTypeError(err_msg)

            # -- check dom_pft vs num_pft
            if max_dom_pft > self.num_pft:
                err_msg = "Please use --crop flag when --dompft is above 16."
                raise argparse.ArgumentTypeError(err_msg)

            # -- check dom_pft vs MAX_pft
            if self.num_pft - 1 < max_dom_pft < NUM_PFT:
                logger.info(
                    "WARNING, you trying to run with generic crops (16 PFT surface dataset)"
                )

            # -- check if all dom_pft are in the same range:
            if min_dom_pft < NAT_PFT <= max_dom_pft:
                err_msg = (
                    "You are subsetting using mixed land units that have both "
                    "natural pfts and crop cfts. Check your surface dataset. "
                )
                raise argparse.ArgumentTypeError(err_msg)

    def check_nonveg(self):
        """
        A function to check at least one of the following arguments is given:
        --include-nonveg
        --dompft DOMPFT

        Basically, this function raises an error
        when zero out non veg land units (by default true) and not provide a dominant pft:

        The user can run ./subset_data using:
        ./subset_data point --dompft
        ./subset_data point --include-nonveg
        ./subset_data point --dompft --include-nonveg

        But this will raise an error:
        ./subset_data point

        By default include_nonveg = False, which means that it zeros out the non-veg landunits.
        """

        if not self.include_nonveg:
            if self.dom_pft is None:
                err_msg = """
                \n
                By default, this will zero out non-veg land units.
                To include non-veg land units, you need to specify --include-nonveg flag.
                To zero-out non-veg land units, you need to specify --dompft.

                You should specify at least one of the following arguments:
                --dompft DOMPFT
                --include-nonveg
                """
                raise argparse.ArgumentTypeError(err_msg)

    def check_pct_pft(self):
        """
        A function to error check pct_pft and calculate it if necessary.

        If the user gives dom_pft and pct_pft :
        - Check if length of dom_pft and pct_pft matches.
          For example, --dompft 8 --pctpft 0.4 0.6 should give an error.

        - Check if the sum of pct_pft is equal to 100% or 1.
          For example, --dompft 8 14 --pctpft 0.6 0.9 should give an error.

        - If the sum of pct_pft is 1, convert it to % (multiply by 100)

        If the user gives one or more dom_pft but no pct_pft, assume equal pct_pft:
        - pct_pft = 100 / number of given dom_pft
          For example, if two dom_pft (s) are given, each of them is 50%.

        """

        # -- if both dom_pft and pct_pft is given:
        if self.dom_pft and self.pct_pft:

            # -- check if the same number of values are given
            if len(self.dom_pft) != len(self.pct_pft):
                err_msg = "Please provide the same number of inputs for --dompft and --pctpft."
                raise argparse.ArgumentTypeError(err_msg)

            # -- check if the sum of pct_pft is equal to 1 or 100
            if sum(self.pct_pft) != 1 and sum(self.pct_pft) != 100:
                err_msg = "Sum of --pctpft values should be equal to 1 or 100."
                raise argparse.ArgumentTypeError(err_msg)

            # -- convert franction to percentage
            if sum(self.pct_pft) == 1:
                self.pct_pft = [pct * 100 for pct in self.pct_pft]

        # -- if the user did not give --pctpft at all (assume equal percentage)
        elif self.dom_pft:
            pct = 100 / len(self.dom_pft)
            self.pct_pft = [pct for pft in self.dom_pft]

        # -- if the user only gave --pctpft with no --dompft
        elif self.pct_pft:
            err_msg = """
                      \n
                      --pctpft is specfied without --dompft.
                      Please specify your dominant pft by --dompft.
                      """
            raise argparse.ArgumentTypeError(err_msg)

        logger.info(" - dominant pft(s) : %s", self.dom_pft)
        logger.info(" - percentage of dominant pft(s) : %s", self.pct_pft)

    def create_domain_at_point(self, indir, file):
        """
        Create domain file for this SinglePointCase class.
        """
        logger.info("----------------------------------------------------------------------")
        logger.info("Creating domain file at %s, %s.", str(self.plon), str(self.plat))

        # specify files
        fdomain_in = os.path.join(indir, file)
        fdomain_out = add_tag_to_filename(fdomain_in, self.tag)
        logger.info("fdomain_in:  %s", fdomain_in)
        logger.info("fdomain_out: %s", os.path.join(self.out_dir, fdomain_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fdomain_in, "xc", "yc", "ni", "nj")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(ni=self.plon, nj=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["nj", "ni"])

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fdomain_in

        wfile = os.path.join(self.out_dir, fdomain_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fdomain_out) %s", wfile)
        f_in.close()
        f_out.close()

    def create_landuse_at_point(self, indir, file, user_mods_dir):
        """
        Create landuse file at a single point.
        """
        logger.info("----------------------------------------------------------------------")
        logger.info(
            "Creating land use file at %s, %s.",
            str(self.plon),
            str(self.plat),
        )

        # specify files
        fluse_in = os.path.join(indir, file)
        fluse_out = add_tag_to_filename(fluse_in, self.tag, replace_res=True)
        logger.info("fluse_in:  %s", fluse_in)
        logger.info("fluse_out: %s", os.path.join(self.out_dir, fluse_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fluse_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["lsmlat", "lsmlon"])

        # specify dimension order
        f_out = f_out.transpose("time", "cft", "natpft", "lsmlat", "lsmlon", "numurbl")

        # revert expand dimensions of YEAR
        year = np.squeeze(np.asarray(f_out["YEAR"]))
        temp_xr = xr.DataArray(year, coords={"time": f_out["time"]}, dims="time", name="YEAR")
        temp_xr.attrs["units"] = "unitless"
        temp_xr.attrs["long_name"] = "Year of PFT data"
        f_out["YEAR"] = temp_xr

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fluse_in

        wfile = os.path.join(self.out_dir, fluse_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fluse_out), %s", wfile)
        f_in.close()
        f_out.close()

        # write to user_nl_clm data if specified
        if self.create_user_mods:
            with open(os.path.join(user_mods_dir, "user_nl_clm"), "a") as nl_clm:
                line = "flanduse_timeseries = '${}'".format(os.path.join(USRDAT_DIR, fluse_out))
                self.write_to_file(line, nl_clm)

    def modify_surfdata_atpoint(self, f_orig):
        """
        Function to modify surface dataset based on the user flags chosen.
        """
        f_mod = f_orig.copy(deep=True)

        # -- modify surface data properties
        if self.dom_pft is not None:
            max_dom_pft = max(self.dom_pft)
            # -- First initialize everything:
            if max_dom_pft < NAT_PFT:
                f_mod["PCT_NAT_PFT"][:, :, :] = 0
            else:
                f_mod["PCT_CFT"][:, :, :] = 0

            # Do we need to initialize these here?
            # Because we set them in include_nonveg
            # f_mod["PCT_NATVEG"][:, :] = 0
            # f_mod["PCT_CROP"][:, :] = 0

            # -- loop over all dom_pft and pct_pft
            iterable_length = len(self.dom_pft)
            cth_to_zip = ensure_iterable(self.cth, iterable_length)
            cbh_to_zip = ensure_iterable(self.cbh, iterable_length)
            zip_pfts = zip(self.dom_pft, self.pct_pft, cth_to_zip, cbh_to_zip)
            for dom_pft, pct_pft, cth, cbh in zip_pfts:
                if cth is not None:
                    f_mod["MONTHLY_HEIGHT_TOP"][:, :, :, dom_pft] = cth
                    f_mod["MONTHLY_HEIGHT_BOT"][:, :, :, dom_pft] = cbh
                if dom_pft < NAT_PFT:
                    f_mod["PCT_NAT_PFT"][:, :, dom_pft] = pct_pft
                else:
                    dom_pft = dom_pft - NAT_PFT
                    f_mod["PCT_CFT"][:, :, dom_pft] = pct_pft

        # -------------------------------
        # By default include_nonveg=False
        # When we use --include-nonveg we turn it to True
        # Therefore by default we are hitting the following if:

        if not self.include_nonveg:
            logger.info("Zeroing out non-vegetation land units in the surface data.")
            f_mod["PCT_LAKE"][:, :] = 0.0
            f_mod["PCT_WETLAND"][:, :] = 0.0
            f_mod["PCT_URBAN"][:, :, :] = 0.0
            f_mod["PCT_GLACIER"][:, :] = 0.0
            f_mod["PCT_OCEAN"][:, :] = 0.0

            if self.dom_pft is not None:
                max_dom_pft = max(self.dom_pft)
                if max_dom_pft < NAT_PFT:
                    f_mod["PCT_NATVEG"][:, :] = 100
                    f_mod["PCT_CROP"][:, :] = 0
                else:
                    f_mod["PCT_NATVEG"][:, :] = 0
                    f_mod["PCT_CROP"][:, :] = 100
            else:
                # -- recalculate percentages after zeroing out non-veg landunits
                # -- so they add up to 100%.
                tot_pct = f_mod["PCT_CROP"] + f_mod["PCT_NATVEG"]
                f_mod["PCT_CROP"] = f_mod["PCT_CROP"] / tot_pct * 100
                f_mod["PCT_NATVEG"] = f_mod["PCT_NATVEG"] / tot_pct * 100

        if self.evenly_split_cropland:
            f_mod["PCT_CFT"][:, :, :] = 100.0 / f_mod["PCT_CFT"].shape[2]

        else:
            logger.info(
                "You chose --include-nonveg --> \
                Do not zero non-vegetation land units in the surface data."
            )

        if self.uni_snow:
            f_mod["STD_ELEV"][:, :] = 20.0
        if self.cap_saturation:
            f_mod["FMAX"][:, :] = 0.0

        return f_mod

    def create_surfdata_at_point(self, indir, file, user_mods_dir, specify_fsurf_out):
        """
        Create surface data file at a single point.
        """
        # pylint: disable=too-many-statements
        logger.info("----------------------------------------------------------------------")
        logger.info(
            "Creating surface dataset file at %s, %s",
            str(self.plon),
            str(self.plat),
        )

        # specify file
        fsurf_in = os.path.join(indir, file)
        if specify_fsurf_out is None:
            fsurf_out = add_tag_to_filename(fsurf_in, self.tag, replace_res=True)
        else:
            fsurf_out = specify_fsurf_out
        logger.info("fsurf_in:  %s", fsurf_in)
        logger.info("fsurf_out: %s", os.path.join(self.out_dir, fsurf_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fsurf_in, "LONGXY", "LATIXY", "lsmlon", "lsmlat")

        # extract gridcell closest to plon/plat
        f_tmp = f_in.sel(lsmlon=self.plon, lsmlat=self.plat, method="nearest")

        # expand dimensions
        f_tmp = f_tmp.expand_dims(["lsmlat", "lsmlon"]).copy(deep=True)

        f_out = self.modify_surfdata_atpoint(f_tmp)

        # specify dimension order
        f_out = f_out.transpose(
            "time",
            "cft",
            "lsmpft",
            "natpft",
            "nglcec",
            "nglcecp1",
            "nlevsoi",
            "nlevurb",
            "numrad",
            "numurbl",
            "lsmlat",
            "lsmlon",
        )

        # update lsmlat and lsmlon to match site specific instead of the nearest point
        # we do this so that if we create user_mods the PTS_LON and PTS_LAT in CIME match
        # the surface data coordinates - which is required
        f_out["lsmlon"] = np.atleast_1d(self.plon)
        f_out["lsmlat"] = np.atleast_1d(self.plat)
        f_out["LATIXY"][:, :] = self.plat
        f_out["LONGXY"][:, :] = self.plon

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fsurf_in

        wfile = os.path.join(self.out_dir, fsurf_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fsurf_out) %s", wfile)
        f_in.close()
        f_tmp.close()
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
        logger.info("----------------------------------------------------------------------")
        logger.info(
            "Creating DATM domain file at %s, %s",
            str(self.plon),
            str(self.plat),
        )

        # specify files
        fdatmdomain_in = os.path.join(datm_tuple.indir, datm_tuple.fdomain_in)
        datm_file = add_tag_to_filename(fdatmdomain_in, self.tag)
        fdatmdomain_out = os.path.join(datm_tuple.outdir, datm_file)
        logger.info("fdatmdomain_in:  %s", fdatmdomain_in)
        logger.info("fdatmdomain out: %s", os.path.join(self.out_dir, fdatmdomain_out))

        # create 1d coordinate variables to enable sel() method
        f_in = self.create_1d_coord(fdatmdomain_in, "xc", "yc", "ni", "nj")

        # extract gridcell closest to plon/plat
        f_out = f_in.sel(ni=self.plon, nj=self.plat, method="nearest")

        # expand dimensions
        f_out = f_out.expand_dims(["nj", "ni"])

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = fdatmdomain_in

        wfile = os.path.join(self.out_dir, fdatmdomain_out)
        self.write_to_netcdf(f_out, wfile)
        logger.info("Successfully created file (fdatmdomain_out) : %s", wfile)
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
        f_out = f_out.transpose("scalar", "time", "lat", "lon")

        # update attributes
        self.update_metadata(f_out)
        f_out.attrs["Created_from"] = file_in

        self.write_to_netcdf(f_out, file_out)
        logger.info("Successfully created file : %s", file_out)
        f_in.close()
        f_out.close()

    def write_shell_commands(self, file, datm_syr, datm_eyr):
        """
        writes out xml commands commands to a file (i.e. shell_commands) for single-point runs
        """
        # write_to_file surrounds text with newlines
        with open(file, "w") as nl_file:
            self.write_to_file("# Change below line if you move the subset data directory", nl_file)
            self.write_to_file("./xmlchange {}={}".format(USRDAT_DIR, self.out_dir), nl_file)
            self.write_to_file("./xmlchange PTS_LON={}".format(str(self.plon)), nl_file)
            self.write_to_file("./xmlchange PTS_LAT={}".format(str(self.plat)), nl_file)
            self.write_to_file("./xmlchange MPILIB=mpi-serial", nl_file)
            if self.create_datm:
                self.write_to_file(f"./xmlchange DATM_YR_ALIGN={datm_syr}", nl_file)
                self.write_to_file(f"./xmlchange DATM_YR_START={datm_syr}", nl_file)
                self.write_to_file(f"./xmlchange DATM_YR_END={datm_eyr}", nl_file)

    def write_datm_streams_lines(self, streamname, datmfiles, file):
        """
        writes out lines for the user_nl_datm_streams file for a specific DATM stream
        for using subset DATM data at a single point

        streamname - stream name (e.g. TPQW)
        datmfiles - comma-separated list (str) of DATM file names
        file - file connection to user_nl_datm_streams file
        """
        self.write_to_file("{}:datafiles={}".format(streamname, ",".join(datmfiles)), file)
        self.write_to_file("{}:mapalgo=none".format(streamname), file)
        self.write_to_file("{}:meshfile=none".format(streamname), file)

    def create_datm_at_point(self, datm_tuple: DatmFiles, datm_syr, datm_eyr, datm_streams_file):
        """
        Create all of a DATM dataset at a point.
        """
        logger.info("----------------------------------------------------------------------")
        logger.info("Creating DATM files at %s, %s", str(self.plon), str(self.plat))

        # --  create data files
        infile = []
        outfile = []
        solarfiles = []
        precfiles = []
        tpqwfiles = []
        for year in range(datm_syr, datm_eyr + 1):
            ystr = str(year)
            for month in range(FIRST_MONTH, LAST_MONTH + 1):
                mstr = str(month)
                if month < 10:
                    mstr = "0" + mstr

                dtag = ystr + "-" + mstr

                fsolar = os.path.join(
                    datm_tuple.indir,
                    datm_tuple.dir_solar,
                    "{}{}.nc".format(datm_tuple.tag_solar, dtag),
                )
                fsolar2 = "{}{}.{}.nc".format(datm_tuple.tag_solar, self.tag, dtag)
                fprecip = os.path.join(
                    datm_tuple.indir,
                    datm_tuple.dir_prec,
                    "{}{}.nc".format(datm_tuple.tag_prec, dtag),
                )
                fprecip2 = "{}{}.{}.nc".format(datm_tuple.tag_prec, self.tag, dtag)
                ftpqw = os.path.join(
                    datm_tuple.indir,
                    datm_tuple.dir_tpqw,
                    "{}{}.nc".format(datm_tuple.tag_tpqw, dtag),
                )
                ftpqw2 = "{}{}.{}.nc".format(datm_tuple.tag_tpqw, self.tag, dtag)

                outdir = os.path.join(self.out_dir, datm_tuple.outdir)
                infile += [fsolar, fprecip, ftpqw]
                outfile += [
                    os.path.join(outdir, fsolar2),
                    os.path.join(outdir, fprecip2),
                    os.path.join(outdir, ftpqw2),
                ]
                solarfiles.append(
                    os.path.join("${}".format(USRDAT_DIR), datm_tuple.outdir, fsolar2)
                )
                precfiles.append(
                    os.path.join("${}".format(USRDAT_DIR), datm_tuple.outdir, fprecip2)
                )
                tpqwfiles.append(os.path.join("${}".format(USRDAT_DIR), datm_tuple.outdir, ftpqw2))

        for idx, out_f in enumerate(outfile):
            logger.debug(out_f)
            self.extract_datm_at(infile[idx], out_f)

        logger.info("All DATM files are created in: %s", datm_tuple.outdir)

        # write to user_nl_datm_streams if specified
        if self.create_user_mods:
            with open(datm_streams_file, "a") as file:
                self.write_datm_streams_lines(datm_tuple.name_solar, solarfiles, file)
                self.write_datm_streams_lines(datm_tuple.name_prec, precfiles, file)
                self.write_datm_streams_lines(datm_tuple.name_tpqw, tpqwfiles, file)
