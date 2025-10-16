#!/usr/bin/env python3

"""
Parent class for some unittest modules relating to gen_mksurfdata_jobscript_single.py
"""

import unittest
import os
import sys
import shutil
from pathlib import Path

import tempfile

from ctsm.path_utils import path_to_ctsm_root


# pylint: disable=too-many-instance-attributes
class TestFGenMkSurfJobscriptSingleParent(unittest.TestCase):
    """Parent class for some unittest modules relating to gen_mksurfdata_jobscript_single.py"""

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
        Path.touch(self._nlfile)
        self.assertTrue(os.path.exists(self._nlfile))
        self._mksurf_exe = os.path.join(self._bld_path, "mksurfdata")
        Path.touch(self._mksurf_exe)
        self.assertTrue(os.path.exists(self._mksurf_exe))
        self._env_mach = os.path.join(self._bld_path, ".env_mach_specific.sh")
        Path.touch(self._env_mach)
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
