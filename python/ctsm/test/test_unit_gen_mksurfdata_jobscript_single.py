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
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_parser
from ctsm.toolchain.gen_mksurfdata_jobscript_single import write_runscript_part1

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


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
        self._bld_path = os.path.join(self._tempdir, "tools_bld")
        os.makedirs(self._bld_path)
        self._nlfile = os.path.join(self._tempdir, "namelist_file")
        self._sys_argv = [
            "gen_mksurfdata_jobscript_single",
            "--bld-path",
            self._bld_path,
            "--namelist-file",
            self._nlfile,
            "--jobscript-file",
            self._jobscript_file,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def add_args(self, machine, nodes, tasks):
        """ " add arguments"""
        self._sys_argv.append(
            ["--machine", machine, "--number-of-nodes", nodes, "--tasks-per-node", tasks]
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_simple_derecho_args(self):
        """test simple derecho arguments"""
        machine = "derecho"
        nodes = 1
        tasks = 64
        self.add_args(machine, nodes, tasks)
        args = get_parser()
        with open(self._jobscript_file, "w", encoding="utf-8") as runfile:
            attribs = write_runscript_part1(nodes, tasks, machine, self._account, runfile)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
