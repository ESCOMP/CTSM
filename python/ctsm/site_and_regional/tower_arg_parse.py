"""
Argument parser to use throughout run_tower.py
"""

import argparse
import logging
import os
import sys

# Get the ctsm util tools and then the cime tools.
_CTSM_PYTHON = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position, import-error, unused-import, wrong-import-order
from ctsm import add_cime_to_path
from ctsm.utils import parse_isoduration
from CIME.utils import parse_args_and_handle_standard_logging_options
from CIME.utils import setup_standard_logging_options


# TODO: Refactor to shorten this and remove this pylint disable!
# pylint: disable=too-many-statements
def get_parser(args, description, valid_neon_sites, valid_plumber_sites):
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    setup_standard_logging_options(parser)

    parser.print_usage = parser.print_help

    parser.add_argument(
        "--neon-sites",
        help="4-letter neon site code.",
        action="store",
        required=False,
        choices=valid_neon_sites + ["all"] + [None],
        dest="neon_sites",
        default=None,
        nargs="+",
    )

    parser.add_argument(
        "--plumber-sites",
        help="six character PLUMBER2 site code (eg, AR-SLu)",
        action="store",
        required=False,
        choices=valid_plumber_sites + ["all"] + [None],
        dest="plumber_sites",
        default=None,
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
        "--no-input-data-check",
        "--no-check-input-data",
        help="""
                Don't check for input data. Implies --setup-only.
                [default: %(default)s]
                """,
        action="store_true",
        dest="no_input_data_check",
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
        choices=["ad", "postad", "transient"],  # , "sasu"],
        # default ad for plumber, transient for neon
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

    parser.add_argument(
        "--xmlchange",
        help="""
                Any xmlchanges (e.g., CLM_CO2_TYPE=constant,CCSM_CO2_PPMV=500)
                [default: %(default)s]
                """,
        required=False,
        type=str,
        default=None,
    )

    args = parse_args_and_handle_standard_logging_options(args, parser)

    if args.neon_sites:
        if "all" in args.neon_sites:
            neon_sites = valid_neon_sites
        else:
            neon_sites = args.neon_sites
            for site in neon_sites:
                if site not in valid_neon_sites:
                    raise ValueError("Invalid neon site name {}".format(site))
        if args.run_type is None:
            args.run_type = "transient"
    else:
        neon_sites = None
    if args.plumber_sites:
        if "all" in args.plumber_sites:
            plumber_sites = valid_plumber_sites
        else:
            plumber_sites = args.plumber_sites
            for site in plumber_sites:
                if site not in valid_plumber_sites:
                    raise ValueError("Invalid plumber site name {}".format(site))
        if (
            not args.neon_sites
        ):  # default to not changing neon behavior if both neon and plumber are present
            if args.run_type is None:
                args.run_type = "ad"
    else:
        plumber_sites = None

    if "CIME_OUTPUT_ROOT" in args.output_root:
        args.output_root = None

    base_case_root = None
    if args.base_case_root:
        base_case_root = os.path.abspath(args.base_case_root)
        if not os.path.exists(base_case_root):
            raise ValueError("Base case root does not exist: {}".format(base_case_root))

    # Reduce output level for this script unless --debug or
    # --verbose is provided on the command line
    if not args.debug and not args.verbose:
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.WARN)

    # --no-input-data-check implies --setup-only
    if args.no_input_data_check and not args.setup_only:
        args.setup_only = True

    return (
        neon_sites,
        plumber_sites,
        args.output_root,
        args.run_type,
        args.experiment,
        args.prism,
        args.overwrite,
        base_case_root,
        args.run_from_postad,
        args.setup_only,
        args.no_batch,
        args.rerun,
        args.no_input_data_check,
        args.user_version,
        args.xmlchange,
    )
