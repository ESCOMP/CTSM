#!/usr/bin/env python3
"""
gen_mksurfdata_jobscript_single.py generates a jobscript for running the
mksurfdata executable to generate a single fsurdat file. For detailed
instructions, see README.
"""
import os
import sys
import argparse
import logging

_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                            os.pardir,
                            os.pardir,
                            'python')
sys.path.insert(1, _CTSM_PYTHON)

from ctsm import add_cime_to_path  # pylint: disable=unused-import
from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from CIME.XML.env_mach_specific import EnvMachSpecific
from CIME.BuildTools.configure import FakeCase

logger = logging.getLogger(__name__)

def get_parser():
    """
    Get parser object for this script.
    """
    # set up logging allowing user control
    setup_logging_pre_config()

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help
    add_logging_args(parser)

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
        help="""currently this recognizes cheyenne, casper, izumi (default
                %(default)s); this needs to be a cime machine, i.e. a machine
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
        help="""output jobscript file to be submitted on cheyenne (default: %(default)s)""",
        action="store",
        dest="jobscript_file",
        required=False,
        default="mksurfdata_jobscript_single"
    )
    return parser

def main ():
    """
    See docstring at the top.
    """
    # --------------------------
    # Obtain input args
    # --------------------------
    args = get_parser().parse_args()
    process_logging_args(args)
    namelist_file = args.namelist_file
    jobscript_file = args.jobscript_file
    number_of_nodes = args.number_of_nodes
    tasks_per_node = args.tasks_per_node
    machine = args.machine
    account = args.account

    # --------------------------
    # Write run script (part 1)
    # --------------------------
    with open(jobscript_file, "w",encoding='utf-8') as runfile:
        runfile.write('#!/bin/bash \n')
        runfile.write('# Edit the batch directives for your batch system \n')
        runfile.write(f'# Below are default batch directives for {machine} \n')
        runfile.write('#PBS -N mksurfdata \n')
        runfile.write('#PBS -j oe \n')
        runfile.write('#PBS -k eod \n')
        runfile.write('#PBS -S /bin/bash \n')
        if machine == 'cheyenne':
            attribs = {'mpilib': 'default'}
            runfile.write('#PBS -l walltime=30:00 \n')
            runfile.write(f"#PBS -A {account} \n")
            runfile.write('#PBS -q regular \n')
            runfile.write(f"#PBS -l select={number_of_nodes}:ncpus=tasks_per_node}:mpiprocs={tasks_per_node} \n")
        elif machine == 'casper':
            attribs = {'mpilib': 'default'}
            runfile.write('#PBS -l walltime=1:00:00 \n')
            runfile.write(f"#PBS -A {account} \n")
            runfile.write('#PBS -q casper \n')
            runfile.write(f'#PBS -l select={number_of_nodes}:ncpus={tasks_per_node}:' \
                          f'mpiprocs={tasks_per_node}:mem=80GB \n')
        elif machine == 'izumi':
            attribs = {'mpilib': 'mvapich2'}
            runfile.write('#PBS -l walltime=2:00:00 \n')
            runfile.write('#PBS -q medium \n')
            runfile.write(f'#PBS -l nodes={number_of_nodes}:ppn={tasks_per_node},mem=555GB -r n \n')
            tool_path = os.path.dirname(os.path.abspath(__file__))
            runfile.write("\n")
            runfile.write(f'cd {tool_path} \n')

        runfile.write("\n")

        # --------------------------
        # Obtain mpirun command from env_mach_specific.xml
        # --------------------------
        bld_path = './tool_bld'
        # Get the ems_file object with standalone_configure=True
        # and the fake_case object with mpilib=attribs['mpilib']
        # so as to use the get_mpirun function pointing to fake_case
        ems_file = EnvMachSpecific(bld_path, standalone_configure=True)
        fake_case = FakeCase(compiler=None, mpilib=attribs['mpilib'],
                             debug=False, comp_interface=None)
        total_tasks = int(tasks_per_node) * int(number_of_nodes)
        cmd = ems_file.get_mpirun(fake_case, attribs, job='name',
                                  overrides = {"total_tasks" : total_tasks,})
        # cmd is a tuple:
        # cmd[0] contains the mpirun command (eg mpirun, mpiexe, etc) as string
        # cmd[1] contains a list of strings that we append as options to cmd[0]
        # The replace function removes unnecessary characters that appear in
        # some such options
        executable = f'{cmd[0]} {" ".join(cmd[1])}'.replace('ENV{', ''). \
                                                    replace('}', '')

        mksurfdata_path = os.path.join(bld_path, 'mksurfdata')
        env_mach_path = os.path.join(bld_path, '.env_mach_specific.sh')

        # --------------------------
        # Write run script (part 2)
        # --------------------------
        runfile.write('# Run env_mach_specific.sh to control the machine ' \
                      'dependent environment including the paths to ' \
                      'compilers and libraries external to cime such as netcdf')
        runfile.write(f'\n. {env_mach_path}\n')
        check = f'if [ $? != 0 ]; then echo "Error running env_mach_specific"; exit -4; fi'
        runfile.write(f"{check} \n")
        runfile.write('# Edit the mpirun command to use the MPI executable ' \
                      'on your system and the arguments it requires \n')
        output = f'{executable} {mksurfdata_path} < {namelist_file}'
        runfile.write(f"{output} \n")
        logger.info('run command is %s', output)

        check = f'if [ $? != 0 ]; then echo "Error running resolution {res}"; exit -4; fi'
        runfile.write(f"{check} \n")
        runfile.write(f"echo Successfully ran resolution\n")

    print (f"echo Successfully created jobscript {jobscript_file}\n")
    sys.exit(0)

if __name__ == "__main__":
    main()
