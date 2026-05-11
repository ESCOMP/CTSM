#!/usr/bin/env python3

"""System tests for gen_mksurfdata_jobscript_single"""

import os

import unittest
import tempfile
import shutil
import sys

from ctsm.path_utils import path_to_ctsm_root
from ctsm.toolchain.gen_mksurfdata_jobscript_single import main
from ctsm import unit_testing

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysGenMkSurfJSSingle(unittest.TestCase):
    """System tests for gen_mksurfdata_jobscript_single"""

    def setUp(self):
        """Setp temporary directory to make the files in"""
        self._original_wd = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)
        self.outfile = "jobscript.sh"
        self.namelist = "res.namelist"
        bld_path = os.path.join(
            path_to_ctsm_root(), "python", "ctsm", "test", "testinputs", "mksurfdata_esmf_bld"
        )
        sys.argv = [
            "gen_mksurfdata_jobscript_single",
            "--bld-path",
            bld_path,
            "--namelist",
            self.namelist,
            "--jobscript-file",
            self.outfile,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._original_wd)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def createJS(self, nodes, tasks_per_node, option_list=None):
        """
        Create a JobScript by sending a list of options in
        """
        if option_list is None:
            option_list = []
        if len(option_list) > 1:
            sys.argv.extend(option_list)
        sys.argv.extend(
            [
                "--number-of-nodes",
                nodes,
                "--tasks-per-node",
                tasks_per_node,
            ]
        )
        print(sys.argv)
        main()
        self.assertTrue(os.path.exists(self.outfile), "Output jobscript file should exist")

    def test_simple_jobscript_single(self):
        """
        Test that a standard simple namelist works
        """
        self.createJS(nodes="4", tasks_per_node="12")

    def test_casper_jobscript_single(self):
        """
        Test that a standard simple namelist works for casper
        """
        opt_list = ["--machine", "casper"]
        self.createJS(nodes="4", tasks_per_node="12", option_list=opt_list)

    def test_izumi_jobscript_single(self):
        """
        Test that a standard simple namelist works for asper
        """
        opt_list = ["--machine", "izumi"]
        self.createJS(nodes="4", tasks_per_node="12", option_list=opt_list)

    def test_bad_bld_path(self):
        """
        Test aborts if the input bld-path does NOT exist
        """
        with self.assertRaisesRegex(SystemExit, "Input Build path"):
            self.createJS(nodes="4", tasks_per_node="12", option_list=["--bld-path", "zztop"])

    def test_neg_nodes(self):
        """
        Test aborts if the input node count is negative
        """
        with self.assertRaisesRegex(
            SystemExit,
            "Input argument --number_of_nodes is zero or negative and needs to be positive",
        ):
            self.createJS(nodes="-4", tasks_per_node="12")

    def test_neg_tasks(self):
        """
        Test aborts if the input tasks_per_node is zero or negative
        """
        with self.assertRaisesRegex(
            SystemExit,
            "Input argument --tasks_per_node is zero or negative and needs to be positive",
        ):
            self.createJS(nodes="4", tasks_per_node="0")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
