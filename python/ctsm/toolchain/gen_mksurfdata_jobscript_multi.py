"""
gen_mksurfdata_jobscript_multi.py generates a jobscript for running the
mksurfdata executable to generate many fsurdat files at once. For detailed
instructions, see README.
"""

import os
import sys
import logging

from ctsm.utils import abort
from ctsm.toolchain.gen_mksurfdata_namelist import main as main_nml
from ctsm.ctsm_logging import process_logging_args
from ctsm.toolchain.gen_mksurfdata_jobscript_single import base_get_parser
from ctsm.toolchain.gen_mksurfdata_jobscript_single import check_parser_args
from ctsm.toolchain.gen_mksurfdata_jobscript_single import write_runscript_part1
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_mpirun


logger = logging.getLogger(__name__)

valid_scenarios = [
    "global-potveg",
    "global-present",
    "global-present-low-res",
    "global-present-ultra-hi-res",
    "global-hist-1850-f19",
    "global-hist-1850-f45",
    "crop-tropics-present",
    "crop",
    "crop-global-present",
    "crop-global-present-low-res",
    "crop-global-present-ne16",
    "crop-global-present-ne30",
    "crop-global-present-ne120",
    "crop-global-present-mpasa480",
    "crop-global-present-nldas",
    "crop-global-1850",
    "crop-global-1850-low-res",
    "crop-global-1850-ne16",
    "crop-global-1850-ne30",
    "crop-global-1850-ne120",
    "crop-global-1850-mpasa480",
    "crop-global-hist",
    "crop-global-hist-low-res",
    "crop-global-hist-ne16",
    "crop-global-hist-ne30",
    "crop-global-hist-f09",
    "crop-global-SSP1-1.9-f09",
    "crop-global-SSP1-2.6-f09",
    "crop-global-SSP2-4.5-f09",
    "crop-global-SSP2-4.5-f19",
    "crop-global-SSP2-4.5-f10",
    "crop-global-SSP2-4.5-f45",
    "crop-global-SSP2-4.5-ne0np4",
    "crop-global-SSP2-4.5-ne3",
    "crop-global-SSP2-4.5-ne16",
    "crop-global-SSP2-4.5-ne30",
    "crop-global-SSP2-4.5-hcru",
    "crop-global-SSP2-4.5-C96",
    "crop-global-SSP2-4.5-mpasa120",
    "crop-global-SSP3-7.0-f09",
    "crop-global-SSP4-3.4-f09",
    "crop-global-SSP4-6.0-f09",
    "crop-global-SSP5-3.4-f09",
    "crop-global-SSP5-8.5-f09",
]


def get_parser():
    """
    Get parser object for this script.
    """
    parser = base_get_parser(default_js_name="mksurfdata_jobscript_multi.sh")

    parser.add_argument(
        "--scenario",
        help="""scenario""",
        choices=valid_scenarios,
        action="store",
        dest="scenario",
        required=True,
    )

    return parser


def write_runscript(
    *,
    args,
    scenario,
    jobscript_file,
    number_of_nodes,
    tasks_per_node,
    account,
    walltime,
    target_list,
    resolution_dict,
    dataset_dict,
    runfile,
):
    """
    Write run script
    """
    # --------------------------
    # Write batch header (part 1)
    # --------------------------
    name = f"mksrf_{scenario}"
    attribs = write_runscript_part1(
        number_of_nodes=number_of_nodes,
        tasks_per_node=tasks_per_node,
        machine=args.machine,
        account=account,
        walltime=walltime,
        runfile=runfile,
        descrip=scenario,
        name=name,
    )
    # --------------------------
    # Obtain mpirun command from env_mach_specific.xml
    # --------------------------
    (executable, mksurfdata_path, env_mach_path) = get_mpirun(args, attribs)

    # Run env_mach_specific.sh to control the machine dependent
    # environment including the paths to compilers and libraries
    # external to cime such as netcdf
    runfile.write(". " + env_mach_path + "\n")
    check = "if [ $? != 0 ]; then echo 'Error running env_specific_script'; exit -4; fi"
    runfile.write(f"{check} \n")
    for target in target_list:
        res_set = dataset_dict[target][1]
        if res_set not in resolution_dict:
            abort(f"Resolution is not in the resolution_dict: {res_set}")
        for res in resolution_dict[res_set]:
            namelist = f"{scenario}_{res}.namelist"
            command = os.path.join(os.getcwd(), "gen_mksurfdata_namelist")
            command = command + " " + dataset_dict[target][0] + " " + res
            command = command + " --silent"
            command = command + f" --namelist {namelist}"
            print(f"command is {command}")
            sys.argv = [x for x in command.split(" ") if x]
            main_nml()
            print(f"generated namelist {namelist}")
            output = f"{executable} {mksurfdata_path} < {namelist}"
            runfile.write(f"{output} \n")
            check = f"if [ $? != 0 ]; then echo 'Error running resolution {res}'; exit -4; fi"
            runfile.write(f"{check} \n")
            runfile.write(f"echo Successfully ran resolution {res}\n")

    runfile.write(f"echo Successfully ran {jobscript_file}\n")


def main():
    """
    See docstring at the top.
    """
    # --------------------------
    # Obtain input args
    # --------------------------
    args = get_parser().parse_args()
    process_logging_args(args)
    check_parser_args(args)
    scenario = args.scenario
    jobscript_file = args.jobscript_file
    number_of_nodes = args.number_of_nodes
    tasks_per_node = args.tasks_per_node
    account = args.account
    walltime = args.walltime

    # --------------------------
    # Determine target list
    # --------------------------
    target_list = [scenario]

    # --------------------------
    # Error checking
    # --------------------------
    for scenario_list in target_list:
        if scenario_list not in valid_scenarios:
            abort("Input scenario is NOT in valid_scenarios")
    # --------------------------
    # Determine resolution sets that are referenced in commands
    # TODO slevis: When new resolutions become supported in ccs_config, the
    # first entry will change to
    # "standard_res_no_crop": [
    # "0.9x1.25",
    # "1.9x2.5",
    # "mpasa60",
    # "mpasa60-3conus",
    # "mpasa60-3centralUS",
    # ],
    # --------------------------
    resolution_dict = {
        "standard_res_no_crop": ["0.9x1.25", "1.9x2.5", "mpasa60"],
        "f09": ["0.9x1.25"],
        "f19": ["1.9x2.5"],
        "hcru": ["360x720cru"],
        "C96": ["C96"],
        "mpasa120": ["mpasa120"],
        "f10": ["10x15"],
        "f45": ["4x5"],
        "low_res_no_crop": ["4x5", "10x15"],
        "ultra_hi_res_no_crop": ["mpasa15", "mpasa3p75"],
        "standard_res": ["360x720cru", "0.9x1.25", "1.9x2.5", "C96", "mpasa120"],
        "standard_res_no_f09": ["360x720cru", "1.9x2.5", "C96", "mpasa120"],
        "low_res": ["4x5", "10x15", "ne3np4.pg3"],
        "mpasa480": ["mpasa480"],
        "nldas_res": ["0.125nldas2"],
        "5x5_amazon": ["5x5_amazon"],
        "ne3": ["ne3np4.pg3"],
        "ne16": ["ne16np4.pg3"],
        "ne30": ["ne30np4.pg3", "ne30np4.pg2", "ne30np4"],
        "ne0np4": [
            "ne0np4.ARCTICGRIS.ne30x8",
            "ne0np4.ARCTIC.ne30x4",
            "ne0np4CONUS.ne30x8",
            "ne0np4.POLARCAP.ne30x4",
        ],
        "ne120": [
            "ne0np4.ARCTICGRIS.ne30x8",
            "ne0np4.ARCTIC.ne30x4",
            "ne0np4CONUS.ne30x8",
            "ne0np4.POLARCAP.ne30x4",
            "ne120np4.pg3",
        ],
    }

    # --------------------------
    # Determine commands for each target list
    # --------------------------
    dataset_dict = {
        "global-potveg": (
            "--start-year 1850 --end-year 1850 --nocrop --potveg               --res",
            "f09",
        ),
        "global-present": (
            "--start-year 2000 --end-year 2000 --nocrop                        --res",
            "standard_res_no_crop",
        ),
        "global-present-low-res": (
            "--start-year 2000 --end-year 2000 --nocrop                        --res",
            "low_res_no_crop",
        ),
        "global-present-ultra-hi-res": (
            "--start-year 2000 --end-year 2000 --nocrop                        --res",
            "ultra_hi_res_no_crop",
        ),
        "global-hist-1850-f19": (
            "--start-year 1850 --end-year 2023 --nocrop --res",
            "f19",
        ),
        "global-hist-1850-f45": (
            "--start-year 1850 --end-year 2023 --nocrop --res",
            "f45",
        ),
        "crop-tropics-present": (
            "--start-year 2000 --end-year 2000                                 --res",
            "5x5_amazon",
        ),
        "crop-global-present": (
            "--start-year 2000 --end-year 2000                                 --res",
            "standard_res",
        ),
        "crop-global-present-low-res": (
            "--start-year 2000 --end-year 2000                                 --res",
            "low_res",
        ),
        "crop-global-present-ne16": (
            "--start-year 2000 --end-year 2000                                 --res",
            "ne16",
        ),
        "crop-global-present-ne30": (
            "--start-year 2000 --end-year 2000                                 --res",
            "ne30",
        ),
        "crop-global-present-ne120": (
            "--start-year 2000 --end-year 2000                                 --res",
            "ne120",
        ),
        "crop-global-present-mpasa480": (
            "--start-year 2000 --end-year 2000                                 --res",
            "mpasa480",
        ),
        "crop-global-present-nldas": (
            "--start-year 2000 --end-year 2000                                 --res",
            "nldas_res",
        ),
        "crop-global-1850": (
            "--start-year 1850 --end-year 1850                                 --res",
            "standard_res",
        ),
        "crop-global-1850-low-res": (
            "--start-year 1850 --end-year 1850                                 --res",
            "low_res",
        ),
        "crop-global-1850-ne16": (
            "--start-year 1850 --end-year 1850                                 --res",
            "ne16",
        ),
        "crop-global-1850-ne30": (
            "--start-year 1850 --end-year 1850                                 --res",
            "ne30",
        ),
        "crop-global-1850-ne120": (
            "--start-year 1850 --end-year 1850                                 --res",
            "ne120",
        ),
        "crop-global-1850-mpasa480": (
            "--start-year 1850 --end-year 1850                                 --res",
            "mpasa480",
        ),
        "crop-global-hist": (
            "--start-year 1850 --end-year 2023 --nosurfdata                    --res",
            "standard_res_no_f09",
        ),
        "crop-global-hist-low-res": (
            "--start-year 1850 --end-year 2023 --nosurfdata                    --res",
            "low_res",
        ),
        "crop-global-hist-ne16": (
            "--start-year 1850 --end-year 2023 --nosurfdata                    --res",
            "ne16",
        ),
        "crop-global-hist-ne30": (
            "--start-year 1850 --end-year 2023 --nosurfdata                    --res",
            "ne30",
        ),
        "crop-global-hist-f09": (
            "--start-year 1700 --end-year 2023                                 --res",
            "f09",
        ),
        "crop-global-SSP1-1.9-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP1-1.9 --res",
            "f09",
        ),
        "crop-global-SSP1-2.6-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP1-2.6 --res",
            "f09",
        ),
        "crop-global-SSP2-4.5-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "f09",
        ),
        "crop-global-SSP2-4.5-hcru": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "hcru",
        ),
        "crop-global-SSP2-4.5-f19": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "f19",
        ),
        "crop-global-SSP2-4.5-f10": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "f10",
        ),
        "crop-global-SSP2-4.5-f45": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "f45",
        ),
        "crop-global-SSP2-4.5-ne0np4": (
            "--start-year 1979 --end-year 2026 --ssp-rcp SSP2-4.5 --res",
            "ne0np4",
        ),
        "crop-global-SSP2-4.5-ne3": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "ne3",
        ),
        "crop-global-SSP2-4.5-ne30": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "ne30",
        ),
        "crop-global-SSP2-4.5-ne16": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "ne16",
        ),
        "crop-global-SSP2-4.5-C96": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "C96",
        ),
        "crop-global-SSP2-4.5-mpasa120": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res",
            "mpasa120",
        ),
        "crop-global-SSP3-7.0-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP3-7.0 --res",
            "f09",
        ),
        "crop-global-SSP4-3.4-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP4-3.4 --res",
            "f09",
        ),
        "crop-global-SSP4-6.0-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP4-6.0 --res",
            "f09",
        ),
        "crop-global-SSP5-3.4-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP5-3.4 --res",
            "f09",
        ),
        "crop-global-SSP5-8.5-f09": (
            "--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP5-8.5 --res",
            "f09",
        ),
    }

    # --------------------------
    # Write run script
    # --------------------------
    with open(jobscript_file, "w", encoding="utf-8") as runfile:

        write_runscript(
            args=args,
            scenario=scenario,
            jobscript_file=jobscript_file,
            number_of_nodes=number_of_nodes,
            tasks_per_node=tasks_per_node,
            account=account,
            walltime=walltime,
            target_list=target_list,
            resolution_dict=resolution_dict,
            dataset_dict=dataset_dict,
            runfile=runfile,
        )

    print(f"echo Successfully created jobscript {jobscript_file}\n")
