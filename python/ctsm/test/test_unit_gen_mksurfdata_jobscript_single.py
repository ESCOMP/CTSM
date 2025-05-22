#!/usr/bin/env python3

"""
Unit tests for gen_mksurfdata_jobscript_single.py subroutines:
"""

import unittest
import os
import sys
import shutil

import tempfile

from ctsm import unit_testing
from ctsm.path_utils import path_to_ctsm_root
from ctsm.path_utils import path_to_cime
from ctsm.os_utils import run_cmd_output_on_error
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_parser
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_mpirun
from ctsm.toolchain.gen_mksurfdata_jobscript_single import check_parser_args
from ctsm.toolchain.gen_mksurfdata_jobscript_single import write_runscript_part1


def add_args(machine, nodes, tasks):
    """add arguments to sys.argv"""
    args_to_add = [
        "--machine",
        machine,
        "--number-of-nodes",
        str(nodes),
        "--tasks-per-node",
        str(tasks),
    ]
    for item in args_to_add:
        sys.argv.append(item)


def create_empty_file(filename):
    """create an empty file"""
    os.system("touch " + filename)


# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


# pylint: disable=protected-access
# pylint: disable=too-many-instance-attributes
class TestFGenMkSurfJobscriptSingle(unittest.TestCase):
    """Tests the gen_mksurfdata_jobscript_single subroutines"""

    def setUp(self):
        """Setup for trying out the methods"""
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._testinputs_path = testinputs_path
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)
        self._account = "ACCOUNT_NUMBER"
        self._jobscript_file = "output_jobscript"
        self._output_compare = """#!/bin/bash
# Edit the batch directives for your batch system
# Below are default batch directives for derecho
#PBS -N mksurfdata
#PBS -j oe
#PBS -k eod
#PBS -S /bin/bash
#PBS -l walltime=12:00:00
#PBS -A ACCOUNT_NUMBER
#PBS -q main
#PBS -l select=1:ncpus=128:mpiprocs=64:mem=218GB

# This is a batch script to run a set of resolutions for mksurfdata_esmf input namelist
# NOTE: THIS SCRIPT IS AUTOMATICALLY GENERATED SO IN GENERAL YOU SHOULD NOT EDIT it!!

"""
        self._bld_path = os.path.join(self._tempdir, "tools_bld")
        os.makedirs(self._bld_path)
        self.assertTrue(os.path.isdir(self._bld_path))
        self._nlfile = os.path.join(self._tempdir, "namelist_file")
        create_empty_file(self._nlfile)
        self.assertTrue(os.path.exists(self._nlfile))
        self._mksurf_exe = os.path.join(self._bld_path, "mksurfdata")
        create_empty_file(self._mksurf_exe)
        self.assertTrue(os.path.exists(self._mksurf_exe))
        self._env_mach = os.path.join(self._bld_path, ".env_mach_specific.sh")
        create_empty_file(self._env_mach)
        self.assertTrue(os.path.exists(self._env_mach))
        sys.argv = [
            "gen_mksurfdata_jobscript_single",
            "--bld-path",
            self._bld_path,
            "--namelist-file",
            self._nlfile,
            "--jobscript-file",
            self._jobscript_file,
            "--account",
            self._account,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def assertFileContentsEqual(self, expected, filepath, msg=None):
        """Asserts that the contents of the file given by 'filepath' are equal to
        the string given by 'expected'. 'msg' gives an optional message to be
        printed if the assertion fails.

        Copied from test_unit_job_launcher_no_batch should go to utils!"""

        with open(filepath, "r") as myfile:
            contents = myfile.read()

        self.assertEqual(expected, contents, msg=msg)

    def test_simple_derecho_args(self):
        """test simple derecho arguments"""
        machine = "derecho"
        nodes = 1
        tasks = 64
        add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        check_parser_args(args)
        with open(self._jobscript_file, "w", encoding="utf-8") as runfile:
            attribs = write_runscript_part1(
                number_of_nodes=nodes,
                tasks_per_node=tasks,
                machine=machine,
                account=self._account,
                walltime=args.walltime,
                runfile=runfile,
            )
            self.assertEqual({"mpilib": "default"}, attribs, msg="attribs not as expected")

        self.assertFileContentsEqual(self._output_compare, self._jobscript_file)

    def test_derecho_mpirun(self):
        """
        test derecho mpirun. This would've helped caught a problem we ran into
        It will also be helpful when sumodules are updated to guide to solutions
        to problems
        """
        machine = "derecho"
        nodes = 4
        tasks = 128
        add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        check_parser_args(args)
        self.assertEqual(machine, args.machine)
        self.assertEqual(tasks, args.tasks_per_node)
        self.assertEqual(nodes, args.number_of_nodes)
        self.assertEqual(self._account, args.account)
        # Create the env_mach_specific.xml file needed for get_mpirun
        # This will catch problems with our usage of CIME objects
        # Doing this here will also catch potential issues in the gen_mksurfdata_build script
        configure_path = os.path.join(path_to_cime(), "CIME", "scripts", "configure")
        self.assertTrue(os.path.exists(configure_path))
        options = " --macros-format CMake --silent --compiler intel --machine " + machine
        cmd = configure_path + options
        cmd_list = cmd.split()
        run_cmd_output_on_error(
            cmd=cmd_list, errmsg="Trouble running configure", cwd=self._bld_path
        )
        self.assertTrue(os.path.exists(self._env_mach))
        expected_attribs = {"mpilib": "default"}
        with open(self._jobscript_file, "w", encoding="utf-8") as runfile:
            attribs = write_runscript_part1(
                number_of_nodes=nodes,
                tasks_per_node=tasks,
                machine=machine,
                account=self._account,
                walltime=args.walltime,
                runfile=runfile,
            )
            self.assertEqual(attribs, expected_attribs)
            (executable, mksurfdata_path, env_mach_path) = get_mpirun(args, attribs)
            expected_exe = "time mpibind "
            self.assertEqual(executable, expected_exe)
            self.assertEqual(mksurfdata_path, self._mksurf_exe)
            self.assertEqual(env_mach_path, self._env_mach)

    def test_too_many_tasks(self):
        """test trying to use too many tasks"""
        machine = "derecho"
        nodes = 1
        tasks = 129
        add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        check_parser_args(args)
        with open(self._jobscript_file, "w", encoding="utf-8") as runfile:
            with self.assertRaisesRegex(
                SystemExit,
                "Number of tasks per node exceeds the number of processors per node"
                + " on this machine",
            ):
                write_runscript_part1(
                    number_of_nodes=nodes,
                    tasks_per_node=tasks,
                    machine=machine,
                    account=self._account,
                    walltime=args.walltime,
                    runfile=runfile,
                )

    def test_zero_tasks(self):
        """test for fail on zero tasks"""
        machine = "derecho"
        nodes = 5
        tasks = 0
        add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        with self.assertRaisesRegex(
            SystemExit,
            "Input argument --tasks_per_node is zero or negative and needs to be positive",
        ):
            check_parser_args(args)

    def test_bld_build_path(self):
        """test for bad build path"""
        machine = "derecho"
        nodes = 10
        tasks = 64
        add_args(machine, nodes, tasks)
        # Remove the build path directory
        shutil.rmtree(self._bld_path, ignore_errors=True)
        args = get_parser().parse_args()
        with self.assertRaisesRegex(SystemExit, "Input Build path .+ does NOT exist, aborting"):
            check_parser_args(args)

    def test_mksurfdata_exist(self):
        """test fails if mksurfdata does not exist"""
        machine = "derecho"
        nodes = 10
        tasks = 64
        add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        os.remove(self._mksurf_exe)
        with self.assertRaisesRegex(SystemExit, "mksurfdata_esmf executable "):
            check_parser_args(args)

    def test_env_mach_specific_exist(self):
        """test fails if the .env_mach_specific.sh file does not exist"""
        machine = "derecho"
        nodes = 10
        tasks = 64
        add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        os.remove(self._env_mach)
        with self.assertRaisesRegex(SystemExit, "Environment machine specific file"):
            check_parser_args(args)

    def test_bad_machine(self):
        """test bad machine name"""
        machine = "zztop"
        nodes = 1
        tasks = 64
        add_args(machine, nodes, tasks)
        with self.assertRaises(SystemExit):
            get_parser().parse_args()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
