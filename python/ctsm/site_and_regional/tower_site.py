"""
This module includes the definition for the TowerSite class,
which has NeonSite and Plumber2Site child classes. This class defines common
functionalities that are in both NeonSite and Plumber2Site classes.
"""
# -- Import libraries

# -- standard libraries
import os.path
import glob
import logging
import re
import shutil
import sys
import time

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root
from ctsm.utils import abort

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


# pylint: disable=too-many-instance-attributes
class TowerSite:
    """
    Parent class to NeonSite and Plumber2Site classes.
    ...
    Attributes
    ----------
    Methods
    -------
    """

    def __init__(self, name, start_year, end_year, start_month, end_month, finidat):
        """
        Initializes TowerSite with the given arguments.
        Parameters
        ----------
        """
        self.name = name
        self.start_year = int(start_year)
        self.end_year = int(end_year)
        self.start_month = int(start_month)
        self.end_month = int(end_month)
        self.cesmroot = path_to_ctsm_root()
        self.finidat = finidat

    def __str__(self):
        """
        Converts ingredients of the TowerSite to string for printing.
        """
        return "{}\n{}".format(
            str(self.__class__),
            "\n".join(
                (
                    "{} = {}".format(str(key), str(self.__dict__[key]))
                    for key in sorted(self.__dict__)
                )
            ),
        )

    def build_base_case(
        self, cesmroot, output_root, res, compset, user_mods_dirs, overwrite=False, setup_only=False
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
        # pylint: disable=attribute-defined-outside-init
        self.base_case_root = output_root
        # pylint: enable=attribute-defined-outside-init
        if not output_root:
            output_root = os.getcwd()
        case_path = os.path.join(output_root, self.name)

        logger.info("base_case_name : %s", self.name)
        logger.info("user_mods_dir  : %s", user_mods_dirs[0])

        if overwrite and os.path.isdir(case_path):
            print("Removing the existing case at: {}".format(case_path))
            if os.getcwd() == case_path:
                abort("Trying to remove the directory tree that we are in")

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
                        """Existing base case is a historical type and should not be
                        --rerun with the --overwrite option""",
                    )
                else:
                    expect(
                        match is not None,
                        """Existing base case should be a historical type and is not
                        --rerun with the --overwrite option""",
                    )
                # reset the case
                case.case_setup(reset=True)
            case_path = case.get_value("CASEROOT")

            if setup_only:
                return case_path

            print("---- base case build ------")
            print("--- This may take a while and you may see WARNING messages ---")
            # always walk through the build process to make sure it's up to date.
            initial_time = time.time()
            build.case_build(case_path, case=case)
            end_time = time.time()
            total = end_time - initial_time
            print("Time required to building the base case: {} s.".format(total))
            # update case_path to be the full path to the base case
        return case_path

    # pylint: disable=no-self-use
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

    def modify_user_nl(self, case_root, run_type, rundir, site_lines=None):
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
            ] + site_lines

        if user_nl_lines:
            with open(user_nl_fname, "a") as nl_file:
                for line in user_nl_lines:
                    nl_file.write("{}\n".format(line))

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
                "ERROR: spinup must be completed first, could not find directory %s", ref_case_root
            )
            return False

        with Case(ref_case_root) as refcase:
            refrundir = refcase.get_value("RUNDIR")
        case.set_value("RUN_REFDIR", refrundir)
        case.set_value("RUN_REFCASE", os.path.basename(ref_case_root))
        refdate = None
        for reffile in glob.iglob(refrundir + "/{}{}.clm2.r.*.nc".format(self.name, root)):
            m_searched = re.search(r"(\d\d\d\d-\d\d-\d\d)-\d\d\d\d\d.nc", reffile)
            if m_searched:
                refdate = m_searched.group(1)
            symlink_force(reffile, os.path.join(rundir, os.path.basename(reffile)))
        logger.info("Found refdate of %s", refdate)
        if not refdate:
            logger.warning("Could not find refcase for %s", case_root)
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

    # pylint: disable=too-many-statements
    # TODO: This code should be broken up into smaller pieces
    def run_case(
        self,
        base_case_root,
        run_type,
        prism,
        run_length,
        user_version,
        tower_type,
        user_mods_dirs,
        overwrite,
        setup_only,
        no_batch,
        rerun,
        experiment,
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

        if (experiment is not False) and (experiment is not None):
            self.name = self.name + "." + experiment
        case_root = os.path.abspath(os.path.join(base_case_root, "..", self.name + "." + run_type))

        rundir = None
        if os.path.isdir(case_root):
            if overwrite:
                print("---- removing the existing case -------")
                if os.getcwd() == case_root:
                    abort("Trying to remove the directory tree that we are in")

                shutil.rmtree(case_root)
            elif rerun:
                with Case(case_root, read_only=False) as case:
                    rundir = case.get_value("RUNDIR")
                    # For existing case check that the compset name is correct
                    existingcompname = case.get_value("COMPSET")
                    match = re.search("^HIST", existingcompname, flags=re.IGNORECASE)
                    # pylint: disable=undefined-variable
                    if re.search("^HIST", compset, flags=re.IGNORECASE) is None:
                        expect(
                            match is None,
                            """Existing base case is a historical type and should not be
                            --rerun with the --overwrite option""",
                        )
                    # pylint: enable=undefined-variable
                    else:
                        expect(
                            match is not None,
                            """Existing base case should be a historical type and is not
                            --rerun with the --overwrite option""",
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
                logger.warning("Case already exists in %s, not overwritting", case_root)
                return
        if run_type == "postad":
            adcase_root = case_root.replace(".postad", ".ad")
            if not os.path.isdir(adcase_root):
                logger.warning("postad requested but no ad case found in %s", adcase_root)
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
            if tower_type == "NEON":
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
                case.case_setup()
                self.set_ref_case(case)
                case.set_value("STOP_N", run_length)

            # For transient cases STOP will be set in the user_mod_directory
            if run_type == "transient":
                case.case_setup()
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
