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
from ctsm.modify_input_files.fsurdat_modifier import check_no_subgrid_section
from ctsm.modify_input_files.fsurdat_modifier import check_no_varlist_section

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

    def test_subgrid_outofrange(self):
        """test a read of subgrid that's out of range"""
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "101.")
        with self.assertRaisesRegex(SystemExit, "is out of range of 0 to 100 ="):
            read_subgrid(self.config, self.cfg_path)

    def test_subgrid_notsumtohundred(self):
        """test a read of subgrid that's doesn't sum to a hundred"""
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "0.")
        self.config.set(section, "pct_lake", "0.")
        self.config.set(section, "pct_wetland", "0.")
        self.config.set(section, "pct_glacier", "0.")
        self.config.set(section, "pct_natveg", "0.")
        self.config.set(section, "pct_crop", "0.")
        with self.assertRaisesRegex(
            SystemExit, "PCT fractions in subgrid section do NOT sum to a hundred as they should"
        ):
            read_subgrid(self.config, self.cfg_path)

    def test_subgrid_badvar(self):
        """test a read of subgrid for a variable thats not in the list"""
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "badvariable", "100.")
        with self.assertRaisesRegex(SystemExit, "is not a valid variable name. Valid vars ="):
            read_subgrid(self.config, self.cfg_path)

    def test_varlist_varinidealized(self):
        """test a read of varlist for a variable thats in the idealized list"""
        section = "modify_fsurdat_variable_list"
        self.config.set(section, "PCT_SAND", "100.")
        with self.assertRaisesRegex(
            SystemExit,
            "is a special variable handled in the idealized section."
            + " This should NOT be handled in the variiable list section. Special idealized vars =",
        ):
            read_var_list(self.config, self.cfg_path)

    def test_varlist_varinsubgrid(self):
        """test a read of varlist for a variable thats in the subgrid list"""
        section = "modify_fsurdat_variable_list"
        self.config.set(section, "PCT_GLACIER", "100.")
        with self.assertRaisesRegex(
            SystemExit,
            "is a variable handled in the subgrid section."
            + " This should NOT be handled in the variiable list section. Subgrid vars =",
        ):
            read_var_list(self.config, self.cfg_path)

    def test_varlist_monthlyvar(self):
        """test a read of varlist for a variable thats one of the monthly
        variables handled in the dom_pft section"""
        section = "modify_fsurdat_variable_list"
        self.config.set(section, "MONTHLY_LAI", "100.")
        with self.assertRaisesRegex(
            SystemExit,
            "is a variable handled as part of the dom_pft handling."
            + " This should NOT be handled in the variiable list section."
            + " Monthly vars handled this way =",
        ):
            read_var_list(self.config, self.cfg_path)

    def test_subgrid_remove(self):
        """test a read of subgrid when it's section has been removed"""
        section = "modify_fsurdat_subgrid_fractions"
        self.config.remove_section(section)
        with self.assertRaisesRegex(SystemExit, "Config file does not have the expected section"):
            read_subgrid(self.config, self.cfg_path)

    def test_subgrid_not_thereifoff(self):
        """test that a graceful error happens if subgrid section is off,
        but it appears in the file"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "process_subgrid_section", "False")
        with self.assertRaisesRegex(SystemExit, "Config file does have a section"):
            check_no_subgrid_section(self.config)

    def test_varlist_not_thereifoff(self):
        """test that a graceful error happens if varlist section is off,
        but it appears in the file"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "process_var_list_section", "False")
        with self.assertRaisesRegex(SystemExit, "Config file does have a section"):
            check_no_varlist_section(self.config)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
