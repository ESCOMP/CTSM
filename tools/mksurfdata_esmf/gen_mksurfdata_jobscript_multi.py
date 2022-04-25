#!/usr/bin/env python3

import sys, os, shutil
import logging
import argparse, textwrap
import subprocess
from datetime import datetime

valid_scenarios=["all",
                 "standard",
                 "global-present",
                 "global-present-T42",
                 "global-present-nldas",
                 "tropics",
                 "crop-tropics-present",
                 "crop-tropics-historical",
                 "crop-tropics-transient",
                 "crop",
                 "crop-global-present",
                 "crop-global-present-ne16np4",
                 "crop-global-present-ne120np4",
                 "crop-global-present-0.125",
                 "crop-global-present-f05",
                 "crop-global-historical",
                 "crop-global-historical-f05",
                 "crop-global-historical-ne120np4",
                 "crop-global-transient",
                 "crop-global-transient-ne120np4",
                 "crop-global-transient-f05",
                 "crop-global-future",
                 "crop-global-SSP1-2.6",
                 "crop-global-SSP3-7.0",
                 "crop-global-SSP5-3.4",
                 "crop-global-SSP2-4.5",
                 "crop-global-SSP1-1.9",
                 "crop-global-SSP4-3.4",
                 "crop-global-SSP4-6.0",
                 "crop-global-SSP5-8.5"]

def get_parser():
    """
    Get parser object for this script.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help

    parser.add_argument(
        '-v', '--verbose',
        help="ncrease output verbosity",
        action="store_true",
    )
    parser.add_argument(
        "--account",
        help="""account number (default P93300606)""",
        action="store",
        dest="account",
        required=False,
        default="P93300606"
    )
    parser.add_argument(
        "--number-of-nodes",
        help="""number of cheyenne nodes requested (required)""",
        action="store",
        dest="number_of_nodes",
        required=True,
    )
    parser.add_argument(
        "--tasks-per-node",
        help="""number of mpi tasks per node for cheyenne requested (default is 12)""",
        action="store",
        dest="tasks_per_node",
        required=False,
        default="12",
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
        help="""output jobscript file to be submitted on cheyenne (default is mksurfdata_jobscript_multi)]""",
        action="store",
        dest="jobscript_file",
        required=False,
        default="mksurfdata_jobscript_multi"
    )
    return parser

def main ():

    # --------------------------
    # Obtain input args
    # --------------------------
    args = get_parser().parse_args()
    scenario = args.scenario
    jobscript_file = args.jobscript_file
    number_of_nodes = args.number_of_nodes
    tasks_per_node = args.tasks_per_node
    account = args.account

    # --------------------------
    # Determine target list 
    # --------------------------
    target_list = []
    if scenario == 'all':
        target_list = ["global-present",
                       "global-present-T42",
                       "global-present-nldas",
                       "crop-global-present",
                       "crop-global-present-ne16np4",
                       "crop-global-present-ne120np4",
                       "crop-global-historical",
                       "crop-global-transient",
                       "crop-tropics-present",
                       "crop-global-SSP1-2.6", 
                       "crop-global-SSP3-7.0",
                       "crop-global-SSP5-3.4",
                       "crop-global-SSP2-4.5",
                       "crop-global-SSP1-1.9",
                       "crop-global-SSP4-3.4",
                       "crop-global-SSP4-6.0",
                       "crop-global-SSP5-8.5"]
    elif scenario == 'standard':
        target_list = ["global-present",
                       "global-present-T42",
                       "global-present-nldas"]
    elif scenario == "crop":
        target_list = ["crop-global-present",
                       "crop-global-present-ne16np4",
                       "crop-global-present-ne120np4",
                       "crop-global-historical",
                       "crop-global-transient"]
    elif scenario == "crop-global-future":
        target_list = ["crop-global-SSP1-2.6", 
                       "crop-global-SSP3-7.0",
                       "crop-global-SSP5-3.4",
                       "crop-global-SSP2-4.5",
                       "crop-global-SSP1-1.9",
                       "crop-global-SSP4-3.4",
                       "crop-global-SSP4-6.0",
                       "crop-global-SSP5-8.5"]
    elif scenario == "tropics":
        target_list = ["crop-tropics-present"]
    else:
        target_list = [scenario]

    # --------------------------
    # Determine resolution sets that are referenced in commands
    # --------------------------
    resolution_dict = {
        "standard_res_no_crop" : ["0.9x1.25","1.9x2.5","10x15"],
        "standard_res"         : ['0.9x1.25','1.9x2.5','10x15','4x5','C96',
                                  'ne30np4','ne30np4.pg2','ne30np4.pg3','ne30np4.pg4','ne120np4.pg3',
                                  'ne0np4.ARCTICGRIS.ne30x8','ne0np4.ARCTIC.ne30x4','ne0np4CONUS.ne30x8'],
        "future_res"           : ["0.9x1.25","1.9x2.5","10x15"],
        "trans_res"            : ['0.9x1.25','1.9x2.5','10x15',
                                  'ne30np4','ne0np4.ARCTICGRIS.ne30x8','ne0np4.ARCTIC.ne30x4','ne0np4CONUS.ne30x8','ne120np4'],
        "T42_res"              : ['T42'],
        "nldas_res"            : ['0.125nldas2'],
        "5x5_amazon_res"       : ['5x5_amazon'],
        "ne16np4_res"          : ['ne120np4'],
        "ne120np4_res"         : ['ne120np4'],
    }

    # --------------------------
    # Determine commands for each target list
    # --------------------------
    dataset_dict={"global-present"                  : ("--start-year 2000 --end-year 2000 --nocrop --vic                  --res", "standard_res_no_crop"),
                  "global-present-T42"              : ("--start-year 2000 --end-year 2000 --nocrop --vic                  --res", "T42_res"),
                  "global-present-nldas"            : ("--start-year 2000 --end-year 2000 --nocrop --vic                  --res", "nldas_res"),
                  "crop-tropics-present"            : ("--start-year 2000 --end-year 2000                                 --res", "5x5_amazon_res"),
                  "crop-global-present"             : ("--start-year 2000 --end-year 2000                                 --res", "standard_res"),
                  "crop-global-present-ne16np4"     : ("--start-year 2000 --end-year 2000                                 --res", "ne16np4_res"),
                  "crop-global-present-ne120np4"    : ("--start-year 2000 --end-year 2000                                 --res", "ne120np4_res"),
                  "crop-global-present-0.125"       : ("--start-year 2000 --end-year 2000 --hirespft                      --res", "nldas_res"),
                  "crop-global-historical"          : ("--start-year 1850 --end-year 1850 --ssp-rcp SSP5-8.5              --res", "standard_res"),
                  "crop-global-historical-ne120np4" : ("--start-year 1850 --end-year 1850 --ssp-rcp SSP5-8.5              --res", "ne120np4_res"),
                  "crop-global-transient"           : ("--start-year 1850 --end-year 2000 --nosurfdata                    --res", "trans_res"),
                  "crop-global-transient-ne120np4"  : ("--start-year 1850 --end-year 2000 --nosurfdata                    --res", "ne120np4_res"),
                  "crop-global-SSP1-1.9"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP1-1.9 --res", "future_res"),
                  "crop-global-SSP1-2.6"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP1-2.6 --res", "future_res"),
                  "crop-global-SSP2-4.5"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP2-4.5 --res", "future_res"),
                  "crop-global-SSP3-7.0"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP3-7.0 --res", "future_res"),
                  "crop-global-SSP4-3.4"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP4-3.4 --res", "future_res"),
                  "crop-global-SSP4-6.0"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP4-6.0 --res", "future_res"),
                  "crop-global-SSP5-3.4"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP4-3.4 --res", "future_res"),
                  "crop-global-SSP5-8.5"            : ("--start-year 1850 --end-year 2100 --nosurfdata --ssp-rcp SSP5-8.5 --res", "future_res")
                  }

    # --------------------------
    # Write run script
    # --------------------------
    with open(jobscript_file, "w",encoding='utf-8') as runfile:

        runfile.write('#!/bin/bash \n')
        runfile.write(f"#PBS -A {account} \n")
        runfile.write('#PBS -N mksurfdata \n')
        runfile.write('#PBS -j oe \n')
        runfile.write('#PBS -q regular \n')
        runfile.write('#PBS -l walltime=30:00 \n')
        runfile.write(f"#PBS -l select={number_of_nodes}:ncpus=36:mpiprocs={tasks_per_node} \n")
        runfile.write("\n")

        np = int(tasks_per_node) * int(number_of_nodes)

        for target in target_list:
            res_set = dataset_dict[target][1]
            for res in resolution_dict[res_set]:
                command = os.path.join(os.getcwd(), "gen_mksurfdata_namelist.py")
                command = command + " " + dataset_dict[target][0] + " " + res
                print (f"command is {command}")
                commands = [x for x in command.split(' ') if x]
                run_cmd = subprocess.run(commands, check=True, capture_output=True)
                if run_cmd.returncode != 0:
                    print ("Error in calling gen_mksurfdata_namelist.py")
                    sys.exit(1)
                output = run_cmd.stdout.decode('utf-8').strip()
                namelist = output.split(' ')[-1]
                print (f"generated namelist {namelist}")
                output = f"mpiexec_mpt -p \"%g:\" -np {np} ./src/mksurfdata < {namelist}"
                runfile.write(f"{output} \n")

    print (f"Successfully created jobscript {jobscript_file}")
    sys.exit(0)

if __name__ == "__main__":
    main()
