#!/usr/bin/env python3

"""
Unit tests for fsurdat_modifier subroutines: read_subgrid, read_varlist
"""

import unittest
import os
import sys
import shutil

import tempfile
from configparser import ConfigParser

from ctsm import unit_testing
from ctsm.path_utils import path_to_ctsm_root
from ctsm.modify_input_files.fsurdat_modifier import fsurdat_modifier_arg_process
from ctsm.modify_input_files.fsurdat_modifier import read_subgrid
from ctsm.modify_input_files.fsurdat_modifier import read_var_list

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


class TestFSurdatModifier(unittest.TestCase):
    """Tests the read_subgrid and read_var_list methods"""

    def setUp(self):
        """Setup for trying out the methods"""
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._cfg_file_path = os.path.join(testinputs_path, "modify_fsurdat_opt_sections.cfg")
        self._testinputs_path = testinputs_path
        self._fsurdat_in = os.path.join(
            testinputs_path,
            "surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self._fsurdat_in = os.path.join(
            testinputs_path,
            "surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc",
        )
        self._fsurdat_out = os.path.join(self._tempdir, "fsurdat_out.nc")
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
            "-i",
            self._fsurdat_in,
            "-o",
            self._fsurdat_out,
        ]
        parser = fsurdat_modifier_arg_process()
        self.cfg_path = str(parser.cfg_path)
        self.config = ConfigParser()
        self.config.read(self.cfg_path)

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_read_subgrid(self):
        """test a simple read of subgrid"""
        read_subgrid(self.config, self.cfg_path)

    def test_read_var_list(self):
        """test a simple read of var_list"""
        read_var_list(self.config, self.cfg_path)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
