"""
This module contains the NeonSite class and class functions which are used in run_neon.py
"""

# Import libraries
import glob
import logging
import os
import re
import shutil
import sys
import time

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# -- import local classes for this script
from ctsm.site_and_regional.base_site import BaseSite
from ctsm.site_and_regional.tower_site import TowerSite

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.path_utils import path_to_ctsm_root

from CIME import build
from CIME.case import Case
from CIME.utils import safe_copy, expect, symlink_force

logger = logging.getLogger(__name__)


# pylint: disable=too-many-instance-attributes
class NeonSite(TowerSite):
    """
    A class for encapsulating neon sites.
    """

    def __init__(self, name, start_year, end_year, start_month, end_month, finidat):
        super().__init__(name, start_year, end_year, start_month, end_month, finidat)

    # pylint: disable=too-many-statements
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
                """hist_fincl1 = 'TOTECOSYSC', 'TOTECOSYSN', 'TOTSOMC', 'TOTSOMN', 'TOTVEGC',
                                 'TOTVEGN', 'TLAI', 'GPP', 'CPOOL', 'NPP', 'TWS', 'H2OSNO'""",
            ]

        if user_nl_lines:
            with open(user_nl_fname, "a") as nl_file:
                for line in user_nl_lines:
                    nl_file.write("{}\n".format(line))
