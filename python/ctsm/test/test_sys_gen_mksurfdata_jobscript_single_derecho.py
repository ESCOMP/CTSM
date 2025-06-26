#!/usr/bin/env python3

"""
System tests for gen_mksurfdata_jobscript_single.py subroutines on Derecho
"""

import unittest
import os

from ctsm import unit_testing
from ctsm.test_gen_mksurfdata_jobscript_single_parent import TestFGenMkSurfJobscriptSingleParent
from ctsm.path_utils import path_to_cime
from ctsm.os_utils import run_cmd_output_on_error
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_parser
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_mpirun
from ctsm.toolchain.gen_mksurfdata_jobscript_single import check_parser_args
from ctsm.toolchain.gen_mksurfdata_jobscript_single import write_runscript_part1


# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


# pylint: disable=protected-access
# pylint: disable=too-many-instance-attributes
class TestFGenMkSurfJobscriptSingleDerecho(TestFGenMkSurfJobscriptSingleParent):
    """Tests the gen_mksurfdata_jobscript_single subroutines on Derecho"""

    def test_derecho_mpirun(self):
        """
        test derecho mpirun. This would've helped caught a problem we ran into
        It will also be helpful when sumodules are updated to guide to solutions
        to problems
        """
        machine = "derecho"
        nodes = 4
        tasks = 128
        unit_testing.add_machine_node_args(machine, nodes, tasks)
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


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
