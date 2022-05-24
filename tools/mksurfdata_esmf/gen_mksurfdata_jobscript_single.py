#!/usr/bin/env python3
"""
gen_mksurfdata_jobscript_single.py generates a jobscript for running the
mksurfdata executable to generate a single fsurdat file. For detailed
instructions, see README.
"""
import sys
import argparse

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
        help="""account number (default: %(default)s)""",
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
        help="""default: %(default)s""",
        action="store",
        dest="machine",
        required=False,
        choices=['cheyenne', 'casper'],  # TODO Plan on adding izumi here
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
        help="""output jobscript file to be submitted on cheyenne (default: %(default)s)""",
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
        runfile.write('# Edit the batch directives for your batch system \n')
        runfile.write('# Below are the batch directives used on cheyenne \n')
        runfile.write(f"#PBS -A {account} \n")
        runfile.write('#PBS -N mksurfdata \n')
        runfile.write('#PBS -j oe \n')
        runfile.write('#PBS -l walltime=30:00 \n')
        if machine == 'cheyenne':
            runfile.write('#PBS -q regular \n')
            runfile.write(f"#PBS -l select={number_of_nodes}:ncpus=36:mpiprocs={tasks_per_node} \n")
        elif machine == 'casper':
            runfile.write('#PBS -q casper \n')
            runfile.write(f"#PBS -l select={number_of_nodes}:ncpus=12:mpiprocs={tasks_per_node} \n")
        runfile.write("\n")

        n_p = int(tasks_per_node) * int(number_of_nodes)

        # Run env_mach_specific.sh to control the machine dependent environment
        # including the paths to compilers and libraries external to cime such
        # as netcdf
        runfile.write('. ./.env_mach_specific.sh \n')
        runfile.write('# Edit the mpirun command to use the MPI executable ' \
                      'on your system and the arguments it requires \n')
        if machine == 'cheyenne':
            output = f"mpiexec_mpt -p \"%g:\" -np {n_p} ./bld/mksurfdata < {namelist_file}"
        elif machine == 'casper':
            output = f"mpiexec -np {n_p} ./bld/mksurfdata < {namelist_file}"
        runfile.write(f"{output} \n")

    print (f"Successfully created jobscript {jobscript_file}")
    sys.exit(0)

if __name__ == "__main__":
    main()
