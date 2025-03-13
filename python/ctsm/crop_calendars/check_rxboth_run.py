"""
Check the results of a run with prescribed sowing dates and maturity requirements
"""

import sys
import argparse
import glob
import os
import numpy as np

# Import the CTSM Python utilities.
# sys.path.insert() is necessary for RXCROPMATURITY to work. The fact that it's calling this script
# in the RUN phase seems to require the python/ directory to be manually added to path.
_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)
import ctsm.crop_calendars.cropcal_module as cc  # pylint: disable=wrong-import-position
from ctsm.crop_calendars.check_rx_obeyed import (  # pylint: disable=wrong-import-position
    check_rx_obeyed,
)
from ctsm.crop_calendars.check_constant_vars import (  # pylint: disable=wrong-import-position
    check_constant_vars,
)


def main(argv):
    """
    Main method: Check the results of a run with prescribed sowing dates and maturity requirements
    """
    # Set arguments
    parser = argparse.ArgumentParser(description="ADD DESCRIPTION HERE")
    parser.add_argument(
        "-d", "--directory", help="Directory with CLM output history files", required=True
    )
    parser.add_argument(
        "--rx_sdates_file", "--rx-sdates-file", help="Prescribed sowing dates file", required=True
    )
    parser.add_argument(
        "--rx_gdds_file",
        "--rx-gdds-file",
        help="Prescribed maturity requirements file",
        required=True,
    )
    parser.add_argument(
        "-y1",
        "--first_usable_year",
        "--first-usable-year",
        type=int,
        help="First usable year in the outputs",
        required=True,
    )
    parser.add_argument(
        "-yN",
        "--last_usable_year",
        "--last-usable-year",
        type=int,
        help="Last usable year in the outputs",
        required=True,
    )
    args = parser.parse_args(argv)

    # Note that _PERHARV will be stripped off upon import
    my_vars = [
        "GRAINC_TO_FOOD_PERHARV",
        "GRAINC_TO_FOOD_ANN",
        "SDATES",
        "SDATES_PERHARV",
        "SYEARS_PERHARV",
        "HDATES",
        "HYEARS",
        "GDDHARV_PERHARV",
        "GDDACCUM_PERHARV",
        "HUI_PERHARV",
        "SOWING_REASON_PERHARV",
        "HARVEST_REASON_PERHARV",
    ]

    any_bad = False

    annual_outfiles = glob.glob(os.path.join(args.directory, "*.clm2.h1.*.nc"))

    # These should be constant in a Prescribed Calendars (rxboth) run, as long as the inputs were
    # static.
    case = {
        "const_vars": ["SDATES", "GDDHARV"],
        "rx_sdates_file": args.rx_sdates_file,
        "rx_gdds_file": args.rx_gdds_file,
    }

    case["ds"], any_bad_import_output = cc.import_output(
        annual_outfiles,
        my_vars=my_vars,
        year_1=args.first_usable_year,
        year_n=args.last_usable_year,
        throw_errors=False,
    )
    any_bad = any_bad or any_bad_import_output

    _, any_bad_check_const_vars = check_constant_vars(
        case["ds"], case, ignore_nan=True, verbose=True, throw_error=True
    )
    any_bad = any_bad or any_bad_check_const_vars

    # Import GGCMI sowing and harvest dates, and check sims
    casename = "Prescribed Calendars"
    gdd_min = None
    if "rx_sdates_file" in case:
        if case["rx_sdates_file"]:
            case["rx_sdates_ds"] = cc.import_rx_dates("sdate", case["rx_sdates_file"], case["ds"])
        if case["rx_gdds_file"]:
            case["rx_gdds_ds"] = cc.import_rx_dates("gdd", case["rx_gdds_file"], case["ds"])

        # Equalize lons/lats
        lonlat_tol = 1e-4
        for ds_name in ["rx_sdates_ds", "rx_gdds_ds"]:
            if ds_name in case:
                for coord_name in ["lon", "lat"]:
                    max_diff_orig = np.max(
                        np.abs(case[ds_name][coord_name].values - case["ds"][coord_name].values)
                    )
                    if max_diff_orig > lonlat_tol:
                        raise RuntimeError(
                            f"{ds_name} {coord_name} values differ too much ({max_diff_orig} > "
                            + f"{lonlat_tol})"
                        )
                    if max_diff_orig > 0:
                        case[ds_name] = case[ds_name].assign_coords(
                            {coord_name: case["ds"][coord_name].values}
                        )
                        max_diff = np.max(
                            np.abs(case[ds_name][coord_name].values - case["ds"][coord_name].values)
                        )
                        print(f"{ds_name} {coord_name} max_diff {max_diff_orig} → {max_diff}")
                    else:
                        print(f"{ds_name} {coord_name} max_diff {max_diff_orig}")

        # Check
        if case["rx_sdates_file"]:
            sdate_not_obeyed = check_rx_obeyed(
                case["ds"].vegtype_str.values,
                case["rx_sdates_ds"].isel(time=0),
                case["ds"],
                casename,
                "SDATES",
            )
            any_bad = any_bad or sdate_not_obeyed
        if case["rx_gdds_file"]:
            gdds_not_obeyed = check_rx_obeyed(
                case["ds"].vegtype_str.values,
                case["rx_gdds_ds"].isel(time=0),
                case["ds"],
                casename,
                "GDDHARV",
                gdd_min=gdd_min,
            )
            any_bad = any_bad or gdds_not_obeyed

    if any_bad:
        msg = "\n   ".join(
            [
                "Unexpected behavior in rxboth run:",
                f"any_bad_import_output:    {any_bad_import_output}",
                f"any_bad_check_const_vars: {any_bad_check_const_vars}",
                f"sdate_not_obeyed:         {sdate_not_obeyed}",
                f"gdds_not_obeyed:          {gdds_not_obeyed}",
            ]
        )
        raise RuntimeError(msg)


if __name__ == "__main__":
    main(sys.argv[1:])
