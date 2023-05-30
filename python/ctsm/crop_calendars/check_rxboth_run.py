# %% Setup

import numpy as np
import sys, argparse
import cropcal_module as cc
import glob, os


def main(argv):
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
    myVars = [
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

    annual_outfiles = glob.glob(os.path.join(args.directory, "*.clm2.h1.*.nc"))

    # These should be constant in a Prescribed Calendars (rxboth) run, as long as the inputs were
    # static.
    case = {
        "constantVars": ["SDATES", "GDDHARV"],
        "rx_sdates_file": args.rx_sdates_file,
        "rx_gdds_file": args.rx_gdds_file,
    }

    case["ds"] = cc.import_output(
        annual_outfiles,
        myVars=myVars,
        y1=args.first_usable_year,
        yN=args.last_usable_year,
    )
    cc.check_constant_vars(case["ds"], case, ignore_nan=True, verbose=True, throw_error=True)

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
        for v in ["rx_sdates_ds", "rx_gdds_ds"]:
            if v in case:
                for l in ["lon", "lat"]:
                    max_diff_orig = np.max(np.abs(case[v][l].values - case["ds"][l].values))
                    if max_diff_orig > lonlat_tol:
                        raise RuntimeError(
                            f"{v} {l} values differ too much ({max_diff_orig} > {lonlat_tol})"
                        )
                    elif max_diff_orig > 0:
                        case[v] = case[v].assign_coords({l: case["ds"][l].values})
                        max_diff = np.max(np.abs(case[v][l].values - case["ds"][l].values))
                        print(f"{v} {l} max_diff {max_diff_orig} → {max_diff}")
                    else:
                        print(f"{v} {l} max_diff {max_diff_orig}")

        # Check
        if case["rx_sdates_file"]:
            cc.check_rx_obeyed(
                case["ds"].vegtype_str.values,
                case["rx_sdates_ds"].isel(time=0),
                case["ds"],
                casename,
                "SDATES",
            )
        if case["rx_gdds_file"]:
            cc.check_rx_obeyed(
                case["ds"].vegtype_str.values,
                case["rx_gdds_ds"].isel(time=0),
                case["ds"],
                casename,
                "GDDHARV",
                gdd_min=gdd_min,
            )


if __name__ == "__main__":
    main(sys.argv[1:])
