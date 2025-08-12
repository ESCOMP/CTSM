"""
Tool for changing parameters on CTSM paramfile
"""
import os

from ctsm.args_utils import comma_separated_list
from ctsm.param_utils.paramfile_shared import paramfile_parser_setup

PFTNAME_VAR = "pftname"


def check_arguments(args):
    """
    Check arguments to set_paramfile
    """
    if not os.path.exists(args.input):
        raise FileNotFoundError(args.input)

    # Avoid potentially overwriting canonical files
    if os.path.exists(args.output):
        raise FileExistsError(args.output)


def get_arguments():
    """
    Parse command-line arguments for setting variables on a netCDF file.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:
            - input: Path to the input netCDF file
            - output: Path to the output netCDF file
            - pft: Optional list of PFT names whose values you want to change
    """
    parser, pft_flags = paramfile_parser_setup(
        "Print values of one or more variables from a netCDF file."
    )
    parser.add_argument("-o", "--output", required=True, help="Output netCDF file")
    parser.add_argument(
        *pft_flags,
        help="Comma-separated list of PFTs whose values you want to change (only applies to PFT-specific variables)",
        type=comma_separated_list,
    )
    args = parser.parse_args()

    check_arguments(args)

    return args


def main():
    """
    Main entry point for the script.
    Parses arguments, opens a netCDF file, makes changes, and saves a new netCDF file.
    """
    args = get_arguments()


if __name__ == "__main__":
    main()
