#!/usr/bin/env python3

import sys, os, shutil
import logging
import argparse, textwrap
import subprocess
from datetime import datetime

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
        help="increase output verbosity",
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
        help="""number of mpi tasks per node for cheyenne requested (required)""",
        action="store",
        dest="tasks_per_node",
        required=True,
    )
    parser.add_argument(
        "--machine",
        help="""currently this recognizes cheyenne, casper, izumi (default
                cheyenne); this needs to be a cime machine, i.e. a machine
                that has been ported to cime where you can build a cime model;
                for details see the README in this directory""",
        action="store",
        dest="machine",
        required=False,
        choices=['cheyenne', 'casper', 'izumi'],
        default='cheyenne'
    )
    parser.add_argument(
        "--namelist-file",
        help="""input namelist file (required)""",
        action="store",
        dest="namelist_file",
        required=True,
    )
    parser.add_argument(
        "--jobscript-file",
        help="""output jobscript file to be submitted on cheyenne (default is mksurfdata_jobscript_single)]""",
        action="store",
        dest="jobscript_file",
        required=False,
        default="mksurfdata_jobscript_single"
    )
    return parser

def main ():

    # --------------------------
    # Obtain input args
    # --------------------------
    args = get_parser().parse_args()
    namelist_file = args.namelist_file
    jobscript_file = args.jobscript_file
    number_of_nodes = args.number_of_nodes
    tasks_per_node = args.tasks_per_node
    machine = args.machine
    account = args.account

    # --------------------------
    # Write run script
    # --------------------------
    with open(jobscript_file, "w",encoding='utf-8') as runfile:
        runfile.write('#!/bin/bash \n')
        runfile.write('#PBS -N mksurfdata \n')
        runfile.write('#PBS -j oe \n')
        runfile.write('#PBS -l walltime=30:00 \n')
        if machine == 'cheyenne':
            runfile.write(f"#PBS -A {account} \n")
            runfile.write('#PBS -q regular \n')
            runfile.write(f"#PBS -l select={number_of_nodes}:ncpus=36:mpiprocs={tasks_per_node} \n")
        elif machine == 'casper':
            runfile.write(f"#PBS -A {account} \n")
            runfile.write('#PBS -q casper \n')
            runfile.write(f"#PBS -l select={number_of_nodes}:ncpus=12:mpiprocs={tasks_per_node}:mem=80GB \n")
        elif machine == 'izumi':
            runfile.write('#PBS -q medium \n')
            runfile.write(f'#PBS -l nodes={number_of_nodes}:ppn={tasks_per_node} -r n \n')
        runfile.write("\n")

        tool_path = os.path.dirname(os.path.abspath(__file__))

        runfile.write(f'cd {tool_path} \n')
        runfile.write("\n")
        runfile.write('. ./bld/.env_mach_specific.sh \n')
        runfile.write("\n")

        np = int(tasks_per_node) * int(number_of_nodes)

        if machine == 'cheyenne':
            output = f"mpiexec_mpt -p \"%g:\" -np {np} ./bld/mksurfdata < {namelist_file}"
        elif machine == 'casper':
            output = f"mpiexec -np {np} ./bld/mksurfdata < {namelist_file}"
        elif machine == 'izumi':
            output = f"mpirun -np {np} ./bld/mksurfdata < {namelist_file}"
        runfile.write(f"{output} \n")

    print (f"Successfully created jobscript {jobscript_file}")
    sys.exit(0)

if __name__ == "__main__":
    main()
