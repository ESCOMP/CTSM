"""
gen_mksurfdata_jobscript_single.py generates a jobscript for running the
mksurfdata executable to generate a single fsurdat file. For detailed
instructions, see README.
"""

import os
import argparse
import logging


from ctsm import add_cime_to_path  # pylint: disable=unused-import
from ctsm.ctsm_logging import setup_logging_pre_config, add_logging_args, process_logging_args
from ctsm.utils import abort
from ctsm.path_utils import path_to_ctsm_root
from CIME.XML.env_mach_specific import (  # pylint: disable=import-error,wrong-import-order
    EnvMachSpecific,
)
from CIME.BuildTools.configure import FakeCase  # pylint: disable=import-error,wrong-import-order

logger = logging.getLogger(__name__)


def base_get_parser(default_js_name="mksurfdata_jobscript_single.sh"):
    """
    Get parser object for the gen_mksurfdata_jobscript scripts
    """
    # set up logging allowing user control
    setup_logging_pre_config()

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.print_usage = parser.print_help
    add_logging_args(parser)

    default_account = os.environ.get("ACCOUNT")
    if default_account is None:
        default_account = "P93300641"
    parser.add_argument(
        "--account",
        help="""account number (default: %(default)s)""",
        action="store",
        dest="account",
        required=False,
        default=default_account,
    )
    parser.add_argument(
        "--number-of-nodes",
        help="""number of derecho nodes requested (required)""",
        action="store",
        dest="number_of_nodes",
        type=int,
        required=True,
    )
    parser.add_argument(
        "--bld-path",
        help="""Path to build directory for mksurfdata_esmf""",
        action="store",
        dest="bld_path",
        default=os.path.join(path_to_ctsm_root(), "tools", "mksurfdata_esmf", "tool_bld"),
    )
    parser.add_argument(
        "--tasks-per-node",
        help="""number of mpi tasks per node for derecho requested (required)""",
        action="store",
        dest="tasks_per_node",
        type=int,
        required=False,
        default="128",
    )
    parser.add_argument(
        "--machine",
        help="""currently this recognizes derecho, casper, izumi (default
                %(default)s); this needs to be a cime machine, i.e. a machine
                that has been ported to cime where you can build a cime model;
                for details see the README in this directory""",
        action="store",
        dest="machine",
        required=False,
        choices=["derecho", "casper", "izumi"],
        default="derecho",
    )
    parser.add_argument(
        "--jobscript-file",
        help="""output jobscript file to be submitted with qsub (default: %(default)s)""",
        action="store",
        dest="jobscript_file",
        required=False,
        default=default_js_name,
    )
    parser.add_argument(
        "--walltime",
        help="""Wallclock time for job submission default is 12:00:00)""",
        action="store",
        dest="walltime",
        required=False,
        default="12:00:00",
    )

    return parser


def get_parser():
    """
    Get parser object for this script.
    """
    parser = base_get_parser()
    parser.add_argument(
        "--namelist-file",
        help="""input namelist file (required)""",
        action="store",
        dest="namelist_file",
        required=True,
    )
    return parser


def check_parser_args(args):
    """Checking for the argument parser values"""
    if args.number_of_nodes < 1:
        abort("Input argument --number_of_nodes is zero or negative and needs to be positive")
    if args.tasks_per_node < 1:
        abort("Input argument --tasks_per_node is zero or negative and needs to be positive")
    if not os.path.exists(args.bld_path):
        abort("Input Build path (" + args.bld_path + ") does NOT exist, aborting")

    mksurfdata_path = os.path.join(args.bld_path, "mksurfdata")
    if not os.path.exists(mksurfdata_path):
        abort(
            "mksurfdata_esmf executable ("
            + mksurfdata_path
            + ") does NOT exist in the bld-path, aborting"
        )
    env_mach_path = os.path.join(args.bld_path, ".env_mach_specific.sh")
    if not os.path.exists(env_mach_path):
        abort(
            "Environment machine specific file ("
            + env_mach_path
            + ") does NOT exist in the bld-path, aborting"
        )


def write_runscript_part1(
    *,
    number_of_nodes,
    tasks_per_node,
    machine,
    account,
    walltime,
    runfile,
    descrip="input namelist",
    name="mksurfdata",
):
    """
    Write run script (part 1) Batch headers
    """
    runfile.write("#!/bin/bash\n")
    runfile.write("# Edit the batch directives for your batch system\n")
    runfile.write(f"# Below are default batch directives for {machine}\n")
    runfile.write(f"#PBS -N {name}\n")
    runfile.write("#PBS -j oe\n")
    runfile.write("#PBS -k eod\n")

    runfile.write("#PBS -S /bin/bash\n")
    attribs = None
    if machine == "derecho":
        attribs = {"mpilib": "default"}
        runfile.write(f"#PBS -l walltime={walltime}\n")
        runfile.write(f"#PBS -A {account}\n")
        runfile.write("#PBS -q main\n")
        ncpus = 128
        runfile.write(
            "#PBS -l select="
            + f"{number_of_nodes}:ncpus={ncpus}:mpiprocs={tasks_per_node}:mem=218GB\n"
        )
    elif machine == "casper":
        attribs = {"mpilib": "default"}
        ncpus = 36
        runfile.write(f"#PBS -l walltime={walltime}\n")
        runfile.write(f"#PBS -A {account}\n")
        runfile.write("#PBS -q casper\n")
        runfile.write(
            f"#PBS -l select={number_of_nodes}:ncpus={tasks_per_node}:"
            f"mpiprocs={tasks_per_node}:mem=80GB\n"
        )
    elif machine == "izumi":
        attribs = {"mpilib": "mvapich2"}
        ncpus = 48
        runfile.write(f"#PBS -l walltime={walltime}\n")
        runfile.write("#PBS -q medium\n")
        runfile.write(f"#PBS -l nodes={number_of_nodes}:ppn={tasks_per_node},mem=555GB -r n\n")
        tool_path = os.path.dirname(os.path.abspath(__file__))
        runfile.write("\n")
        runfile.write(f"cd {tool_path}\n")

    runfile.write("\n")
    runfile.write(
        f"# This is a batch script to run a set of resolutions for mksurfdata_esmf {descrip}\n"
    )
    runfile.write(
        "# NOTE: THIS SCRIPT IS AUTOMATICALLY GENERATED "
        + "SO IN GENERAL YOU SHOULD NOT EDIT it!!\n\n"
    )

    # Make sure tasks_per_node doesn't exceed the number of cpus per node
    if tasks_per_node > ncpus:
        abort("Number of tasks per node exceeds the number of processors per node on this machine")
    return attribs


def get_mpirun(args, attribs):
    """
    Get the mpirun command for this machine
    This requires a working env_mach_specific.xml file in the build directory
    """
    bld_path = args.bld_path
    # Get the ems_file object with standalone_configure=True
    # and the fake_case object with mpilib=attribs['mpilib']
    # so as to use the get_mpirun function pointing to fake_case
    ems_file = EnvMachSpecific(bld_path, standalone_configure=True)
    fake_case = FakeCase(compiler=None, mpilib=attribs["mpilib"], debug=False, comp_interface=None)
    total_tasks = int(args.tasks_per_node) * int(args.number_of_nodes)
    cmd = ems_file.get_mpirun(
        fake_case,
        attribs,
        job="name",
        exe_only=True,
        overrides={
            "total_tasks": total_tasks,
        },
    )
    # cmd is a tuple:
    # cmd[0] contains the mpirun command (eg mpirun, mpiexe, etc) as string
    # cmd[1] contains a list of strings that we append as options to cmd[0]
    # The replace function removes unnecessary characters that appear in
    # some such options
    executable = f'time {cmd[0]} {" ".join(cmd[1])}'.replace("ENV{", "").replace("}", "")

    mksurfdata_path = os.path.join(bld_path, "mksurfdata")
    env_mach_path = os.path.join(bld_path, ".env_mach_specific.sh")

    return (executable, mksurfdata_path, env_mach_path)


def write_runscript_part2(namelist_file, runfile, executable, mksurfdata_path, env_mach_path):
    """
    Write run script (part 2)
    """
    runfile.write(
        "# Run env_mach_specific.sh to control the machine "
        "dependent environment including the paths to "
        "compilers and libraries external to cime such as netcdf"
    )
    runfile.write(f"\n. {env_mach_path}\n")
    check = 'if [ $? != 0 ]; then echo "Error running env_mach_specific script"; exit -4; fi'
    runfile.write(f"{check} \n")
    runfile.write(
        "# Edit the mpirun command to use the MPI executable "
        "on your system and the arguments it requires \n"
    )
    output = f"{executable} {mksurfdata_path} < {namelist_file}"
    runfile.write(f"{output} \n")
    logger.info("run command is %s", output)

    check = f'if [ $? != 0 ]; then echo "Error running for namelist  {namelist_file}"; exit -4; fi'
    runfile.write(f"{check} \n")
    runfile.write("echo Successfully ran resolution\n")


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
    namelist_file = args.namelist_file
    jobscript_file = args.jobscript_file
    number_of_nodes = args.number_of_nodes
    tasks_per_node = args.tasks_per_node
    machine = args.machine
    account = args.account
    walltime = args.walltime

    # --------------------------
    # Write to file
    # --------------------------
    with open(jobscript_file, "w", encoding="utf-8") as runfile:
        # --------------------------
        # Write batch header (part 1)
        # --------------------------
        attribs = write_runscript_part1(
            number_of_nodes=number_of_nodes,
            tasks_per_node=tasks_per_node,
            machine=machine,
            account=account,
            walltime=walltime,
            runfile=runfile,
        )
        # --------------------------
        # Obtain mpirun command from env_mach_specific.xml
        # --------------------------
        (executable, mksurfdata_path, env_mach_path) = get_mpirun(args, attribs)
        # --------------------------
        # Write commands to run
        # --------------------------
        write_runscript_part2(namelist_file, runfile, executable, mksurfdata_path, env_mach_path)

    print(f"echo Successfully created jobscript {jobscript_file}\n")
