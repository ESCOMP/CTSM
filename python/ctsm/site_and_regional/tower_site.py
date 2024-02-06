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

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


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
        # pylint: disable=attribute-defined-outside-init
        self.base_case_root = output_root
        # pylint: enable=attribute-defined-outside-init
        user_mods_dirs = [os.path.join(cesmroot, "cime_config", "usermods_dirs", "NEON", self.name)]
        if not output_root:
            output_root = os.getcwd()
        case_path = os.path.join(output_root, self.name)

        logger.info("base_case_name : %s", self.name)
        logger.info("user_mods_dir  : %s", user_mods_dirs[0])

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
            build.case_build(
                case_path, case=case
            )  # TODO: this causes issues if run from tower_site.py
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
