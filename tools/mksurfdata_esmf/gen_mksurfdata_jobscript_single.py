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
        "--mpi-tasks",
        help="""number of mpi tasks requested (required)""",
        action="store",
        dest="mpi_tasks",
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
    mpi_tasks = args.mpi_tasks
    tasks_per_node = args.tasks_per_node
    account = args.account

    # --------------------------
    # Write run script
    # --------------------------
    with open(jobscript_file, "w",encoding='utf-8') as runfile:

        np = int(mpi_tasks)
        nodes = int(np / int(tasks_per_node))

        runfile.write('#!/bin/bash \n')
        runfile.write(f"#PBS -A {account} \n")
        runfile.write('#PBS -N mksurfdata \n')
        runfile.write('#PBS -j oe \n')
        runfile.write('#PBS -q regular \n')
        runfile.write('#PBS -l walltime=30:00 \n')
        runfile.write(f"#PBS -l select={nodes}:ncpus=36:mpiprocs={tasks_per_node} \n")

        runfile.write("\n")
        runfile.write("export TMPDIR=/glade/scratch/$USER/temp \n")
        runfile.write("mkdir -p $TMPDIR \n")
        runfile.write("\n")

        output = f"mpiexec_mpt -p \"%g:\" -np {np} ./src/mksurfdata < {namelist_file}"
        runfile.write(f"{output} \n")

    print (f"Successfully created jobscript {jobscript_file}")
    sys.exit(0)

if __name__ == "__main__":
    main()
