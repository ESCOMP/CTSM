#!/usr/bin/env python3

"""System tests for gen_mksurfdata_jobscript_single

"""

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
        self._tempdir = tempfile.mkdtemp()
        os.chdir(self._tempdir)
        self.outfile = "jobscript.sh"
        self.namelist = "res.namelist"
        sys.argv = [
            "gen_mksurfdata_jobscript_single",
            "--namelist",
            self.namelist,
            "--jobscript-file",
            self.outfile,
        ]

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def createJS(self, nodes, tasks_per_node, options=[]):
        """
        Create a JobScript by sending a list of options in
        """
        if len(options) > 1:
            sys.argv.extend(options)
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
        # pylint: disable=no-self-use
        self.createJS(nodes="4", tasks_per_node="12")


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
