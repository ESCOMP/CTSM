#! /usr/bin/env python3

"""
|------------------------------------------------------------------|
|---------------------  Instructions  -----------------------------|
|------------------------------------------------------------------|
This is a wrapper script for running CTSM simulation for one or more
neon sites.

This script is only for neon site and we will develop a more general
code later.

This script first creates and builds a generic base case.
Next, it will clone the base_case for different neon sites and run
types to reduce the need to build ctsm everytime.

This script will do the following:
    1) Create a generic base case for cloning.
    2) Make the case for the specific neon site(s).
    3) Make changes to the case, for:
        a. AD spinup
	b. post-AD spinup
        c. transient
    	#---------------
    	d. SASU or Matrix spinup
    4) Build and submit the case.

-------------------------------------------------------------------
Instructions for running using conda python environments:

../../py_env_create
conda activate ctsm_py

-------------------------------------------------------------------
To see the available options:
    ./run_neon.py --help
-------------------------------------------------------------------
"""
# TODO (NS)
# - [ ]
# - [ ] Case dependency and the ability to check case status
# - [ ] If Case dependency works we don't need finidat given explicilty for post-ad and transient.

# - [ ] checkout_externals instead of using env varaiable
# - [ ] wget the fields available and run for those available

# - [ ] Matrix spin-up if (SASU) Eric merged it in
# - [ ] Make sure both AD and SASU are not on at the same time

# - [ ] Make sure CIME and other dependencies is checked out.


# Import libraries
import argparse
import datetime
import glob
import logging
import os
import re
import shutil
import sys
import time
import pandas as pd

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import parse_isoduration
from ctsm.download_utils import download_file

from standard_script_setup import *

logger = logging.getLogger(__name__)


def get_parser(args, description, valid_neon_sites):
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    CIME.utils.setup_standard_logging_options(parser)

    parser.print_usage = parser.print_help

    parser.add_argument(
        "--neon-sites",
        help="4-letter neon site code.",
        action="store",
        required=False,
        choices=valid_neon_sites + ["all"],
        dest="neon_sites",
        default=["OSBS"],
        nargs="+",
    )

    parser.add_argument(
        "--base-case",
        help="""
                Root Directory of base case build
                [default: %(default)s]
                """,
        action="store",
        dest="base_case_root",
        type=str,
        required=False,
        default=None,
    )

    parser.add_argument(
        "--output-root",
        help="""
                Root output directory of cases
                [default: %(default)s]
                """,
        action="store",
        dest="output_root",
        type=str,
        required=False,
        default="CIME_OUTPUT_ROOT as defined in cime",
    )

    parser.add_argument(
        "--overwrite",
        help="""
                overwrite existing case directories
                [default: %(default)s]
                """,
        action="store_true",
        dest="overwrite",
        required=False,
        default=False,
    )

    parser.add_argument(
        "--setup-only",
        help="""
                Only setup the requested cases, do not build or run
                [default: %(default)s]
                """,
        action="store_true",
        dest="setup_only",
        required=False,
        default=False,
    )

    parser.add_argument(
        "--rerun",
        help="""
                If the case exists but does not appear to be complete, restart it.
                [default: %(default)s]
                """,
        action="store_true",
        dest="rerun",
        required=False,
        default=False,
    )

    parser.add_argument(
        "--no-batch",
        help="""
                Run locally, do not use batch queueing system (if defined for Machine)
                [default: %(default)s]
                """,
        action="store_true",
        dest="no_batch",
        required=False,
        default=False,
    )

    parser.add_argument(
        "--run-type",
        help="""
                        Type of run to do
                        [default: %(default)s]
                        """,
        choices=["ad", "postad", "transient", "sasu"],
        default="transient",
    )

    parser.add_argument(
        "--prism",
        help="""
                Uses the PRISM reanaylsis precipitation data for the site instead of the NEON data
                (only available over Continental US)
                """,
        action="store_true",
        dest="prism",
        required=False,
        default=False,
    )

    parser.add_argument(
        "--experiment",
        help="""
                Appends the case name with string for model experiment
                """,
        action="store",
        dest="experiment",
        type=str,
        required=False,
        default=None,
    )

    parser.add_argument(
        "--run-length",
        help="""
                How long to run (modified ISO 8601 duration)
                [default: %(default)s]
                """,
        required=False,
        type=str,
        default="0Y",
    )

    parser.add_argument(
        "--start-date",
        help="""
                Start date for running CTSM simulation in ISO format.
                [default: %(default)s]
                (currently non-functional)
                """,
        action="store",
        dest="start_date",
        required=False,
        type=datetime.date.fromisoformat,
        default=datetime.datetime.strptime("2018-01-01", "%Y-%m-%d"),
    )

    parser.add_argument(
        "--end-date",
        help="""
                End date for running CTSM simulation in ISO format.
                [default: %(default)s]
                """,
        action="store",
        dest="end_date",
        required=False,
        type=datetime.date.fromisoformat,
        default=datetime.datetime.strptime("2021-01-01", "%Y-%m-%d"),
    )

    parser.add_argument(
        "--run-from-postad",
        help="""
                        For transient runs only - should we start from the postad spinup or finidat?
                        By default start from finidat, if this flag is used the postad run must be available.
                        """,
        action="store_true",
        required=False,
        default=False,
    )
    parser.add_argument(
        "--neon-version",
        help="""
                Neon data version to use for this simulation.
                [default: use the latest data available]
                """,
        action="store",
        dest="user_version",
        required=False,
        type=str,
        choices=["v1", "v2", "v3"],
    )

    args = CIME.utils.parse_args_and_handle_standard_logging_options(args, parser)

    if "all" in args.neon_sites:
        neon_sites = valid_neon_sites
    else:
        neon_sites = args.neon_sites
        for site in neon_sites:
            if site not in valid_neon_sites:
                raise ValueError("Invalid site name {}".format(site))

    if "CIME_OUTPUT_ROOT" in args.output_root:
        args.output_root = None

    if args.run_length == "0Y":
        if args.run_type == "ad":
            run_length = "100Y"
        elif args.run_type == "postad":
            run_length = "100Y"
        else:
            # The transient run length is set by cdeps atm buildnml to
            # the last date of the available tower data
            # this value is not used
            run_length = "4Y"
    else:
        run_length = args.run_length

    run_length = parse_isoduration(run_length)
    base_case_root = None
    if args.base_case_root:
        base_case_root = os.path.abspath(args.base_case_root)

    # Reduce output level for this script unless --debug or
    # --verbose is provided on the command line
    if not args.debug and not args.verbose:
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.WARN)

    return (
        neon_sites,
        args.output_root,
        args.run_type,
        args.experiment,
        args.prism,
        args.overwrite,
        run_length,
        base_case_root,
        args.run_from_postad,
        args.setup_only,
        args.no_batch,
        args.rerun,
        args.user_version,
    )


class NeonSite:
    """
    A class for encapsulating neon sites.

    ...

    Attributes
    ----------

    Methods
    -------
    """

    def __init__(self, name, start_year, end_year, start_month, end_month, finidat):
        self.name = name
        self.start_year = int(start_year)
        self.end_year = int(end_year)
        self.start_month = int(start_month)
        self.end_month = int(end_month)
        self.cesmroot = path_to_ctsm_root()
        self.finidat = finidat

    def __str__(self):
        return (
            str(self.__class__) + "\n" + "\n".join((str(item) + " = " for item in (self.__dict__)))
        )

    def build_base_case(
        self, cesmroot, output_root, res, compset, overwrite=False, setup_only=False
    ):
        """
        Function for building a base_case to clone.
        To spend less time on building ctsm for the neon cases,
        all the other cases are cloned from this case

        Args:
        self:
            The NeonSite object
        base_root (str):
            root of the base_case CIME
        res (str):
            base_case resolution or gridname
        compset (str):
            base case compset
        overwrite (bool) :
            Flag to overwrite the case if exists
        """
        print("---- building a base case -------")
        self.base_case_root = output_root
        user_mods_dirs = [os.path.join(cesmroot, "cime_config", "usermods_dirs", "NEON", self.name)]
        if not output_root:
            output_root = os.getcwd()
        case_path = os.path.join(output_root, self.name)

        logger.info("base_case_name : {}".format(self.name))
        logger.info("user_mods_dir  : {}".format(user_mods_dirs[0]))

        if overwrite and os.path.isdir(case_path):
            print("Removing the existing case at: {}".format(case_path))
            shutil.rmtree(case_path)

        with Case(case_path, read_only=False) as case:
            if not os.path.isdir(case_path):
                print("---- creating a base case -------")

                case.create(
                    case_path,
                    cesmroot,
                    compset,
                    res,
                    run_unsupported=True,
                    answer="r",
                    output_root=output_root,
                    user_mods_dirs=user_mods_dirs,
                    driver="nuopc",
                )

                print("---- base case created ------")

                # --change any config for base_case:
                # case.set_value("RUN_TYPE","startup")
                print("---- base case setup ------")
                case.case_setup()
            else:
                # For existing case check that the compset name is correct
                existingcompname = case.get_value("COMPSET")
                match = re.search("^HIST", existingcompname, flags=re.IGNORECASE)
                if re.search("^HIST", compset, flags=re.IGNORECASE) is None:
                    expect(
                        match is None,
                        "Existing base case is a historical type and should not be  -- rerun with the --orverwrite option",
                    )
                else:
                    expect(
                        match is not None,
                        "Existing base case should be a historical type and is not -- rerun with the --orverwrite option",
                    )
                # reset the case
                case.case_setup(reset=True)
            case_path = case.get_value("CASEROOT")

            if setup_only:
                return case_path

            print("---- base case build ------")
            print("--- This may take a while and you may see WARNING messages ---")
            # always walk through the build process to make sure it's up to date.
            t0 = time.time()
            build.case_build(case_path, case=case)
            t1 = time.time()
            total = t1 - t0
            print("Time required to building the base case: {} s.".format(total))
            # update case_path to be the full path to the base case
        return case_path

    def diff_month(self):
        """
        Determine difference between two dates in months
        """
        d1 = datetime.datetime(self.end_year, self.end_month, 1)
        d2 = datetime.datetime(self.start_year, self.start_month, 1)
        return (d1.year - d2.year) * 12 + d1.month - d2.month

    def get_batch_query(self, case):
        """
        Function for querying the batch queue query command for a case, depending on the
        user's batch system.

        Args:
        case:
            case object
        """

        if case.get_value("BATCH_SYSTEM") == "none":
            return "none"
        return case.get_value("batch_query")

    def run_case(
        self,
        base_case_root,
        run_type,
        prism,
        run_length,
        user_version,
        overwrite=False,
        setup_only=False,
        no_batch=False,
        rerun=False,
        experiment=False,
    ):
        """
        Run case.

        Args:
        self
        base_case_root: str, opt
            file path of base case
        run_type: str, opt
            transient, post_ad, or ad case, default transient
        prism: bool, opt
            if True, use PRISM precipitation, default False
        run_length: str, opt
            length of run, default '4Y'
        user_version: str, opt
            default 'latest'
        overwrite: bool, opt
            default False
        setup_only: bool, opt
            default False; if True, set up but do not run case
        no_batch: bool, opt
            default False
        rerun: bool, opt
            default False
        experiment: str, opt
            name of experiment, default False
        """
        user_mods_dirs = [
            os.path.join(self.cesmroot, "cime_config", "usermods_dirs", "NEON", self.name)
        ]
        expect(
            os.path.isdir(base_case_root),
            "Error base case does not exist in {}".format(base_case_root),
        )
        # -- if user gives a version:
        if user_version:
            version = user_version
        else:
            version = "latest"

        print("using this version:", version)

        if experiment is not None:
            self.name = self.name + "." + experiment
        case_root = os.path.abspath(os.path.join(base_case_root, "..", self.name + "." + run_type))

        rundir = None
        if os.path.isdir(case_root):
            if overwrite:
                print("---- removing the existing case -------")
                shutil.rmtree(case_root)
            elif rerun:
                with Case(case_root, read_only=False) as case:
                    rundir = case.get_value("RUNDIR")
                    # For existing case check that the compset name is correct
                    existingcompname = case.get_value("COMPSET")
                    match = re.search("^HIST", existingcompname, flags=re.IGNORECASE)
                    if re.search("^HIST", compset, flags=re.IGNORECASE) is None:
                        expect(
                            match is None,
                            "Existing base case is a historical type and should not be  -- rerun with the --orverwrite option",
                        )
                    else:
                        expect(
                            match is not None,
                            "Existing base case should be a historical type and is not -- rerun with the --orverwrite option",
                        )
                    if os.path.isfile(os.path.join(rundir, "ESMF_Profile.summary")):
                        print("Case {} appears to be complete, not rerunning.".format(case_root))
                    elif not setup_only:
                        print("Resubmitting case {}".format(case_root))
                        case.submit(no_batch=no_batch)
                        print("-----------------------------------")
                        print("Successfully submitted case!")
                        batch_query = self.get_batch_query(case)
                        if batch_query != "none":
                            print(f"Use {batch_query} to check its run status")
                    return
            else:
                logger.warning("Case already exists in {}, not overwritting.".format(case_root))
                return

        if run_type == "postad":
            adcase_root = case_root.replace(".postad", ".ad")
            if not os.path.isdir(adcase_root):
                logger.warning("postad requested but no ad case found in {}".format(adcase_root))
                return

        if not os.path.isdir(case_root):
            # read_only = False should not be required here
            with Case(base_case_root, read_only=False) as basecase:
                print("---- cloning the base case in {}".format(case_root))
                #
                # EBK: 11/05/2022 -- Note keeping the user_mods_dirs argument is important. Although
                # it causes some of the user_nl_* files to have duplicated inputs. It also ensures
                # that the shell_commands file is copied, as well as taking care of the DATM inputs.
                # See https://github.com/ESCOMP/CTSM/pull/1872#pullrequestreview-1169407493
                #
                basecase.create_clone(case_root, keepexe=True, user_mods_dirs=user_mods_dirs)

        with Case(case_root, read_only=False) as case:
            if run_type != "transient":
                # in order to avoid the complication of leap years,
                # we always set the run_length in units of days.
                case.set_value("STOP_OPTION", "ndays")
                case.set_value("REST_OPTION", "end")
            case.set_value("CONTINUE_RUN", False)
            case.set_value("NEONVERSION", version)
            if prism:
                case.set_value("CLM_USRDAT_NAME", "NEON.PRISM")

            if run_type == "ad":
                case.set_value("CLM_FORCE_COLDSTART", "on")
                case.set_value("CLM_ACCELERATED_SPINUP", "on")
                case.set_value("RUN_REFDATE", "0018-01-01")
                case.set_value("RUN_STARTDATE", "0018-01-01")
                case.set_value("RESUBMIT", 1)
                case.set_value("STOP_N", run_length)

            else:
                case.set_value("CLM_FORCE_COLDSTART", "off")
                case.set_value("CLM_ACCELERATED_SPINUP", "off")
                case.set_value("RUN_TYPE", "hybrid")

            if run_type == "postad":
                self.set_ref_case(case)
                case.set_value("STOP_N", run_length)

            # For transient cases STOP will be set in the user_mod_directory
            if run_type == "transient":
                if self.finidat:
                    case.set_value("RUN_TYPE", "startup")
                else:
                    if not self.set_ref_case(case):
                        return
                case.set_value("CALENDAR", "GREGORIAN")
                case.set_value("RESUBMIT", 0)
                case.set_value("STOP_OPTION", "nmonths")

            if not rundir:
                rundir = case.get_value("RUNDIR")

            self.modify_user_nl(case_root, run_type, rundir)

            case.create_namelists()
            # explicitly run check_input_data
            case.check_all_input_data()
            if not setup_only:
                case.submit(no_batch=no_batch)
                print("-----------------------------------")
                print("Successfully submitted case!")
                batch_query = self.get_batch_query(case)
                if batch_query != "none":
                    print(f"Use {batch_query} to check its run status")

    def set_ref_case(self, case):
        """
        Set an existing case as the reference case, eg for use with spinup.
        """
        rundir = case.get_value("RUNDIR")
        case_root = case.get_value("CASEROOT")
        if case_root.endswith(".postad"):
            ref_case_root = case_root.replace(".postad", ".ad")
            root = ".ad"
        else:
            ref_case_root = case_root.replace(".transient", ".postad")
            root = ".postad"
        if not os.path.isdir(ref_case_root):
            logger.warning(
                "ERROR: spinup must be completed first, could not find directory {}".format(
                    ref_case_root
                )
            )
            return False

        with Case(ref_case_root) as refcase:
            refrundir = refcase.get_value("RUNDIR")
        case.set_value("RUN_REFDIR", refrundir)
        case.set_value("RUN_REFCASE", os.path.basename(ref_case_root))
        refdate = None
        for reffile in glob.iglob(refrundir + "/{}{}.clm2.r.*.nc".format(self.name, root)):
            m = re.search("(\d\d\d\d-\d\d-\d\d)-\d\d\d\d\d.nc", reffile)
            if m:
                refdate = m.group(1)
            symlink_force(reffile, os.path.join(rundir, os.path.basename(reffile)))
        logger.info("Found refdate of {}".format(refdate))
        if not refdate:
            logger.warning("Could not find refcase for {}".format(case_root))
            return False

        for rpfile in glob.iglob(refrundir + "/rpointer*"):
            safe_copy(rpfile, rundir)
        if not os.path.isdir(os.path.join(rundir, "inputdata")) and os.path.isdir(
            os.path.join(refrundir, "inputdata")
        ):
            symlink_force(os.path.join(refrundir, "inputdata"), os.path.join(rundir, "inputdata"))

        case.set_value("RUN_REFDATE", refdate)
        if case_root.endswith(".postad"):
            case.set_value("RUN_STARTDATE", refdate)
        # NOTE: if start options are set, RUN_STARTDATE should be modified here
        return True

    def modify_user_nl(self, case_root, run_type, rundir):
        """
        Modify user namelist. If transient, include finidat in user_nl;
        Otherwise, adjust user_nl to include different mfilt, nhtfrq, and variables in hist_fincl1.
        """
        user_nl_fname = os.path.join(case_root, "user_nl_clm")
        user_nl_lines = None
        if run_type == "transient":
            if self.finidat:
                user_nl_lines = [
                    "finidat = '{}/inputdata/lnd/ctsm/initdata/{}'".format(rundir, self.finidat)
                ]
        else:
            user_nl_lines = [
                "hist_fincl2 = ''",
                "hist_mfilt = 20",
                "hist_nhtfrq = -8760",
                "hist_empty_htapes = .true.",
                "hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC', 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'",
            ]

        if user_nl_lines:
            with open(user_nl_fname, "a") as fd:
                for line in user_nl_lines:
                    fd.write("{}\n".format(line))


def check_neon_listing(valid_neon_sites):
    """
    A function to download and parse neon listing file.
    """
    listing_file = "listing.csv"
    url = "https://storage.neonscience.org/neon-ncar/listing.csv"

    download_file(url, listing_file)
    available_list = parse_neon_listing(listing_file, valid_neon_sites)
    return available_list


def parse_neon_listing(listing_file, valid_neon_sites):
    """
    A function to parse neon listing file
    and find neon sites with the dates
    where data is available.

    Args:
        listing_file (str): downloaded listing file

    Returns:
        available_list :
            list of neon_site objects that is found
            on the downloaded listing file.
    """

    # pd.set_option("display.max_rows", None, "display.max_columns", None)

    available_list = []

    df = pd.read_csv(listing_file)

    # check for finidat files for transient run
    finidatlist = df[df["object"].str.contains("lnd/ctsm")]

    # -- filter lines with atm/cdep
    df = df[df["object"].str.contains("atm/cdeps/")]

    # -- split the object str to extract site name
    df = df["object"].str.split("/", expand=True)

    # -- groupby site name
    grouped_df = df.groupby(8)
    for key, item in grouped_df:
        # -- check if it is a valid neon site
        if any(key in x for x in valid_neon_sites):
            site_name = key
            tmp_df = grouped_df.get_group(key)

            # -- filter files only ending with YYYY-MM.nc
            tmp_df = tmp_df[tmp_df[9].str.contains("\d\d\d\d-\d\d.nc")]

            # -- find all the data versions
            # versions = tmp_df[7].unique()
            # print ("all versions available for ", site_name,":", *versions)
            latest_version = tmp_df[7].iloc[-1]
            # print ("latests version available for ", site_name,":", latest_version)

            tmp_df = tmp_df[tmp_df[7].str.contains(latest_version)]
            # -- remove .nc from the file names
            tmp_df[9] = tmp_df[9].str.replace(".nc", "", regex=False)

            tmp_df2 = tmp_df[9].str.split("-", expand=True)

            # ignore any prefix in file name and just get year
            tmp_df2[0] = tmp_df2[0].str.slice(-4)

            # -- figure out start_year and end_year
            start_year = tmp_df2[0].iloc[0]
            end_year = tmp_df2[0].iloc[-1]

            # -- figure out start_month and end_month
            start_month = tmp_df2[1].iloc[0]
            end_month = tmp_df2[1].iloc[-1]

            logger.debug("Valid neon site " + site_name + " found!")
            logger.debug("File version {}".format(latest_version))
            logger.debug("start_year={}".format(start_year))
            logger.debug("end_year={}".format(end_year))
            logger.debug("start_month={}".format(start_month))
            logger.debug("end_month={}".format(end_month))
            finidat = None
            for line in finidatlist["object"]:
                if site_name in line:
                    finidat = line.split(",")[0].split("/")[-1]

            neon_site = NeonSite(site_name, start_year, end_year, start_month, end_month, finidat)
            logger.debug(neon_site)
            available_list.append(neon_site)

    return available_list


def main(description):
    """
    Determine valid neon sites. Make an output directory if it does not exist.
    Loop through requested sites and run CTSM at that site.
    """
    cesmroot = path_to_ctsm_root()
    # Get the list of supported neon sites from usermods
    valid_neon_sites = glob.glob(
        os.path.join(cesmroot, "cime_config", "usermods_dirs", "NEON", "[!d]*")
    )
    valid_neon_sites = sorted([v.split("/")[-1] for v in valid_neon_sites])

    (
        site_list,
        output_root,
        run_type,
        experiment,
        prism,
        overwrite,
        run_length,
        base_case_root,
        run_from_postad,
        setup_only,
        no_batch,
        rerun,
        user_version,
    ) = get_parser(sys.argv, description, valid_neon_sites)

    if output_root:
        logger.debug("output_root : " + output_root)
        if not os.path.exists(output_root):
            os.makedirs(output_root)

    # -- check neon listing file for available data:
    available_list = check_neon_listing(valid_neon_sites)

    # =================================
    # -- all neon sites can be cloned from one generic case
    # -- so no need to define a base_case for every site.

    res = "CLM_USRDAT"
    if run_type == "transient":
        compset = "IHist1PtClm51Bgc"
    else:
        compset = "I1PtClm51Bgc"

    # --  Looping over neon sites

    for neon_site in available_list:
        if neon_site.name in site_list:
            if run_from_postad:
                neon_site.finidat = None
            if not base_case_root:
                base_case_root = neon_site.build_base_case(
                    cesmroot, output_root, res, compset, overwrite, setup_only
                )
            logger.info("-----------------------------------")
            logger.info("Running CTSM for neon site : {}".format(neon_site.name))
            neon_site.run_case(
                base_case_root,
                run_type,
                prism,
                run_length,
                user_version,
                overwrite,
                setup_only,
                no_batch,
                rerun,
                experiment,
            )
