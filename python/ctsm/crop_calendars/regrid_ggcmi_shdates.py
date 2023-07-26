from subprocess import run
import os
import shutil
import glob
import argparse
import sys


def run_and_check(cmd):
    result = run(
        cmd,
        shell=True,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Trouble running `{result.args}` in shell:\n{result.stdout}\n{result.stderr}"
        )


# Functionized because these are shared by process_ggcmi_shdates
def define_arguments(parser):
    # Required
    parser.add_argument(
        "-rr",
        "--regrid-resolution",
        help="Target CLM resolution, to be saved in output filenames.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-rt",
        "--regrid-template-file",
        help="Template file to be used in regridding of inputs.",
        type=str,
        required=True,
    )
    return parser


def main(
    regrid_resolution, regrid_template_file_in, regrid_input_directory, regrid_output_directory
):

    # Ensure we can call necessary shell scripts
    for cmd in ["ncks", "ncrename", "ncpdq", "cdo"]:
        run_and_check(f"{cmd} --help")

    os.chdir(regrid_input_directory)
    if not os.path.exists(regrid_output_directory):
        os.makedirs(regrid_output_directory)

    templatefile = os.path.join(regrid_output_directory, "template.nc")

    # For some reason, doing ncks -v directly doesn't work. Have to copy it over first.
    regrid_template_file = os.path.join(
        regrid_output_directory,
        os.path.basename(regrid_template_file_in),
    )
    shutil.copyfile(regrid_template_file_in, regrid_template_file)

    if os.path.exists(templatefile):
        os.remove(templatefile)

    for v in ["planting_day", "maturity_day", "growing_season_length"]:
        run_and_check(f"ncks -A -v area '{regrid_template_file}' '{templatefile}'")
        run_and_check(f"ncrename -v area,{v} '{templatefile}'")
    os.remove(regrid_template_file)
    run_and_check(f"ncpdq -O -h -a -lat '{templatefile}' '{templatefile}'")

    input_files = glob.glob("*nc4")
    input_files.sort()
    for f in input_files:
        print(f[0:6])
        f2 = os.path.join(regrid_output_directory, f)
        f3 = f2.replace(".nc4", f"_nninterp-{regrid_resolution}.nc4")

        if os.path.exists(f3):
            os.remove(f3)

        # Sometimes cdo fails for no apparent reason. In testing this never happened more than 3x in a row.
        try:
            run_and_check(f"cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")
        except:
            try:
                run_and_check(f"cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")
            except:
                try:
                    run_and_check(f"cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")
                except:
                    run_and_check(f"cdo -L -remapnn,'{templatefile}' -setmisstonn '{f}' '{f3}'")

    os.remove(templatefile)


if __name__ == "__main__":
    ###############################
    ### Process input arguments ###
    ###############################
    parser = argparse.ArgumentParser(
        description="Regrids raw sowing and harvest date files provided by GGCMI to a target CLM resolution."
    )

    # Define arguments
    parser = define_arguments(parser)
    parser.add_argument(
        "-i",
        "--regrid-input-directory",
        help="Directory containing the raw GGCMI sowing/harvest date files.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--regrid-output-directory",
        help="Directory where regridded output files should be saved.",
        type=str,
        required=True,
    )

    # Get arguments
    args = parser.parse_args(sys.argv[1:])

    ###########
    ### Run ###
    ###########
    main(
        args.regrid_resolution,
        os.path.realpath(args.regrid_template_file),
        os.path.realpath(args.regrid_input_directory),
        os.path.realpath(args.regrid_output_directory),
    )
