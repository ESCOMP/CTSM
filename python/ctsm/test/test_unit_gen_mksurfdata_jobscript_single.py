#!/usr/bin/env python3

"""
Unit tests for gen_mksurfdata_jobscript_single.py subroutines:
"""

import unittest
import os
import argparse
import sys
import shutil

import tempfile

from ctsm import unit_testing
from ctsm.path_utils import path_to_ctsm_root
from ctsm.toolchain.gen_mksurfdata_jobscript_single import get_parser
from ctsm.toolchain.gen_mksurfdata_jobscript_single import check_parser_args
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
        self._output_compare = """#!/bin/bash 
# Edit the batch directives for your batch system 
# Below are default batch directives for derecho 
#PBS -N mksurfdata 
#PBS -j oe 
#PBS -k eod 
#PBS -S /bin/bash 
#PBS -l walltime=59:00 
#PBS -A ACCOUNT_NUMBER 
#PBS -q main 
#PBS -l select=1:ncpus=64:mpiprocs=64 

"""
        self._bld_path = os.path.join(self._tempdir, "tools_bld")
        os.makedirs(self._bld_path)
        self._nlfile = os.path.join(self._tempdir, "namelist_file")
        sys.argv = [
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
        args_to_add = ["--machine", machine, "--number-of-nodes", str(nodes), "--tasks-per-node", str(tasks)] 
        for item in args_to_add:
           sys.argv.append( item )

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
        
        Copied from test_unit_job_launcher_no_batch should go to utils! """

        with open(filepath, "r") as myfile:
            contents = myfile.read()

        self.assertEqual(expected, contents, msg=msg)

    def test_simple_derecho_args(self):
        """test simple derecho arguments"""
        machine = "derecho"
        nodes = 1
        tasks = 64
        self.add_args(machine, nodes, tasks)
        args = get_parser().parse_args()
        with open(self._jobscript_file, "w", encoding="utf-8") as runfile:
            attribs = write_runscript_part1(nodes, tasks, machine, self._account, runfile)
            self.assertEqual( {"mpilib": "default"}, attribs, msg="attribs not as expected" )

        self.assertFileContentsEqual( self._output_compare, self._jobscript_file )

    def test_bad_machine(self):
        """test bad machine name"""
        machine = "zztop"
        nodes = 1
        tasks = 64
        self.add_args(machine, nodes, tasks)
        with self.assertRaises( SystemExit ):
            args = get_parser().parse_args()

if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
