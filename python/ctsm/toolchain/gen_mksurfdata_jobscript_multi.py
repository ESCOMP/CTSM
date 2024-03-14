"""
gen_mksurfdata_jobscript_multi.py generates a jobscript for running the
mksurfdata executable to generate many fsurdat files at once. For detailed
instructions, see README.
"""
import os
import sys
import argparse

from ctsm.path_utils import path_to_ctsm_root
from ctsm.toolchain.gen_mksurfdata_namelist import main as main_nml
from ctsm.utils import abort

valid_scenarios = [
    "global-potveg",
    "global-present",
    "global-present-low-res",
    "global-present-ultra-hi-res",
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
    "crop-global-SSP1-1.9-f09",
    "crop-global-SSP1-2.6-f09",
    "crop-global-SSP2-4.5-f09",
    "crop-global-SSP2-4.5-f19",
    "crop-global-SSP2-4.5-f10",
    "crop-global-SSP2-4.5-f45",
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
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        "-v",
        "--verbose",
        help="ncrease output verbosity",
        action="store_true",
    )
    parser.add_argument(
        "--account",
        help="""account number (default P93300606)""",
        action="store",
        dest="account",
        required=False,
        default="P93300641",
    )
    parser.add_argument(
        "--bld-path",
        help="""Path to build directory for mksurfdata_esmf""",
        action="store",
        dest="bld_path",
        default=os.path.join(path_to_ctsm_root(), "tools", "mksurfdata_esmf", "tool_bld"),
    )
    parser.add_argument(
        "--number-of-nodes",
        help="""number of derecho nodes requested (required)""",
        action="store",
        dest="number_of_nodes",
        required=True,
    )
    parser.add_argument(
        "--tasks-per-node",
        help="""number of mpi tasks per node for derecho requested""",
        action="store",
        dest="tasks_per_node",
        required=False,
        default="128",
    )
    parser.add_argument(
        "--walltime",
        help="""Wallclock time for job submission default is 12:00:00)""",
        action="store",
        dest="walltime",
        required=False,
        default="12:00:00",
    )
    parser.add_argument(
        "--queue",
        help="""Queue to submit to)""",
        action="store",
        dest="queue",
        required=False,
        default="main",
    )
    parser.add_argument(
        "--scenario",
        help="""scenario""",
        choices=valid_scenarios,
        action="store",
        dest="scenario",
        required=True,
    )
    parser.add_argument(
        "--jobscript-file",
        help="""output jobscript file to be submitted with qsub
                [default: %(default)s]""",
        action="store",
        dest="jobscript_file",
        required=False,
        default="mksurfdata_jobscript_multi.sh",
    )
    return parser


def write_runscript(
    scenario,
    jobscript_file,
    number_of_nodes,
    tasks_per_node,
    account,
    walltime,
    queue,
    target_list,
    resolution_dict,
    dataset_dict,
    env_specific_script,
    mksurfdata,
    runfile,
):
    """
    Write run script
    """
    runfile.write("#!/bin/bash \n")
    runfile.write(f"#PBS -A {account} \n")
    runfile.write(f"#PBS -N mksrf_{scenario} \n")
    runfile.write("#PBS -j oe \n")
    runfile.write("#PBS -k eod \n")
    runfile.write("#PBS -S /bin/bash \n")
    runfile.write(f"#PBS -q {queue} \n")
    runfile.write(f"#PBS -l walltime={walltime} \n")
    runfile.write(
        "#PBS -l select="
        + f"{number_of_nodes}:ncpus={tasks_per_node}:mpiprocs={tasks_per_node}:mem=218GB \n"
    )
    runfile.write(
        f"# This is a batch script to run a set of resolutions for mksurfdata_esmf {scenario}\n"
    )
    runfile.write(
        "# NOTE: THIS SCRIPT IS AUTOMATICALLY GENERATED "
        + "SO IN GENERAL YOU SHOULD NOT EDIT it!!\n\n"
    )
    runfile.write("\n")

    # Run env_mach_specific.sh to control the machine dependent
    # environment including the paths to compilers and libraries
    # external to cime such as netcdf
    runfile.write(". " + env_specific_script + "\n")
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
            output = f"time mpiexec {mksurfdata} < {namelist}"
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
    scenario = args.scenario
    jobscript_file = args.jobscript_file
    number_of_nodes = args.number_of_nodes
    tasks_per_node = args.tasks_per_node
    account = args.account
    walltime = args.walltime
    queue = args.queue

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
        "low_res": ["4x5", "10x15", "ne3np4.pg3"],
        "mpasa480": ["mpasa480"],
        "nldas_res": ["0.125nldas2"],
        "5x5_amazon": ["5x5_amazon"],
        "ne3": ["ne3np4.pg3"],
        "ne16": ["ne16np4.pg3"],
        "ne30": ["ne30np4.pg3", "ne30np4.pg2", "ne30np4"],
        "ne120": [
            "ne0np4.ARCTICGRIS.ne30x8",
            "ne0np4.ARCTIC.ne30x4",
            "ne0np4CONUS.ne30x8",
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
        "crop-tropics-present": (
            "--start-year 2000 --end-year 2000                                 --res",
            "5x5_amazon",
        ),
        "crop-global-present": (
            "--start-year 2000 --end-year 2000 --vic                           --res",
            "standard_res",
        ),
        "crop-global-present-low-res": (
            "--start-year 2000 --end-year 2000 --vic                           --res",
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
            # TODO slevis: --hirespft uses old data for now, so keep out
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
            "--start-year 1850 --end-year 2015 --nosurfdata                    --res",
            "standard_res",
        ),
        "crop-global-hist-low-res": (
            "--start-year 1850 --end-year 2015 --nosurfdata                    --res",
            "low_res",
        ),
        "crop-global-hist-ne16": (
            "--start-year 1850 --end-year 2015 --nosurfdata                    --res",
            "ne16",
        ),
        "crop-global-hist-ne30": (
            "--start-year 1850 --end-year 2015 --nosurfdata                    --res",
            "ne30",
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
    # TODO Here, reuse code from gen_mksurfdata_jobscript_single
    # that's found in the section titled "Obtain mpirun command ..."
    # --------------------------

    # --------------------------
    # Make sure files exist or exit
    # --------------------------
    if not os.path.exists(args.bld_path):
        print(
            args.bld_path
            + " directory does NOT exist -- build mksurdata_esmf before running this script --"
            + " using ./gen_mksurfdata_build.sh"
        )
        sys.exit(1)

    env_specific_script = os.path.join(args.bld_path, ".env_mach_specific.sh")
    if not os.path.exists(env_specific_script):
        print(env_specific_script + " does NOT exist")
        sys.exit(1)
    mksurfdata = os.path.join(args.bld_path, "mksurfdata")
    if not os.path.exists(mksurfdata):
        print(mksurfdata + " does NOT exist")
        sys.exit(1)
    # --------------------------
    # Write run script
    # --------------------------
    with open(jobscript_file, "w", encoding="utf-8") as runfile:

        write_runscript(
            scenario,
            jobscript_file,
            number_of_nodes,
            tasks_per_node,
            account,
            walltime,
            queue,
            target_list,
            resolution_dict,
            dataset_dict,
            env_specific_script,
            mksurfdata,
            runfile,
        )

    print(f"echo Successfully created jobscript {jobscript_file}\n")
