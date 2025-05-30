#!/usr/bin/env python3

"""
Unit tests for fsurdat_modifier subroutines:
"""

import unittest
import os
import sys
import shutil

import tempfile
from configparser import ConfigParser
import xarray as xr

from ctsm import unit_testing
from ctsm.path_utils import path_to_ctsm_root
from ctsm.longitude import Longitude
from ctsm.modify_input_files.fsurdat_modifier import fsurdat_modifier_arg_process
from ctsm.modify_input_files.fsurdat_modifier import read_cfg_subgrid
from ctsm.modify_input_files.fsurdat_modifier import read_cfg_option_control
from ctsm.modify_input_files.fsurdat_modifier import read_cfg_var_list
from ctsm.modify_input_files.fsurdat_modifier import check_no_subgrid_section
from ctsm.modify_input_files.fsurdat_modifier import check_no_varlist_section
from ctsm.modify_input_files.modify_fsurdat import ModifyFsurdat

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


# Allow as many public methods as needed...
# pylint: disable=too-many-public-methods
# Allow all the instance attributes that we need
# pylint: disable=too-many-instance-attributes
class TestFSurdatModifier(unittest.TestCase):
    """Tests the fsurdat_modifier subroutines"""

    def setUp(self):
        """Setup for trying out the methods"""
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._cfg_file_path = os.path.join(testinputs_path, "modify_fsurdat_opt_sections.cfg")
        self._testinputs_path = testinputs_path
        self._fsurdat_in = os.path.join(
            testinputs_path,
            "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031.nc",
        )
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
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
        my_data = xr.open_dataset(self._fsurdat_in)
        lon_type = 360
        self.modify_fsurdat = ModifyFsurdat(
            my_data=my_data,
            lon_1=Longitude(0.0, lon_type),
            lon_2=Longitude(360.0, lon_type),
            lat_1=90.0,
            lat_2=90.0,
            landmask_file=None,
            lat_dimname=None,
            lon_dimname=None,
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_subgrid_and_idealized_fails(self):
        """test that subgrid and idealized fails gracefully"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "True")
        self.config.set(section, "include_nonveg", "False")
        self.config.set(section, "process_subgrid_section", "True")
        self.config.set(section, "dom_pft", "UNSET")
        with self.assertRaisesRegex(
            SystemExit,
            "idealized AND process_subgrid_section can NOT both be on, pick one or the other",
        ):
            read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)

    def test_dompft_and_splitcropland_fails(self):
        """test that setting dompft crop with evenly_split_cropland True fails gracefully"""
        section = "modify_fsurdat_basic_options"
        crop_pft = max(self.modify_fsurdat.file.natpft.values) + 1
        self.config.set(section, "dom_pft", str(crop_pft))
        self.config.set(section, "evenly_split_cropland", "True")
        with self.assertRaisesRegex(
            SystemExit,
            "dom_pft must not be set to a crop PFT when evenly_split_cropland is True",
        ):
            read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)

    def test_optional_only_true_and_false(self):
        """test that optional settings can only be true or false"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "dom_pft", "1")
        varlist = (
            "idealized",
            "include_nonveg",
            "process_subgrid_section",
            "process_var_list_section",
        )
        for var in varlist:
            self.config.set(section, var, "True")
        self.config.set(section, "idealized", "False")
        read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)
        for var in varlist:
            self.config.set(section, var, "False")
        read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)
        self.config.set(section, "dom_pft", "UNSET")
        read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)
        varlist = (
            "idealized",
            "evenly_split_cropland",
        )
        for var in varlist:
            orig_value = self.config.get(section, var)
            self.config.set(section, var, "Thing")
            with self.assertRaisesRegex(
                SystemExit, "Non-boolean value found for .cfg file variable: " + var
            ):
                read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)
            self.config.set(section, var, orig_value)

    def test_read_subgrid(self):
        """test a simple read of subgrid"""
        read_cfg_subgrid(self.config, self.cfg_path)

    def test_read_subgrid_allglacier(self):
        """test a read of subgrid that's for all glacier"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "0. 0. 0.")
        self.config.set(section, "pct_lake", "0.")
        self.config.set(section, "pct_wetland", "0.")
        self.config.set(section, "pct_ocean", "0.")
        self.config.set(section, "pct_glacier", "100.")
        self.config.set(section, "pct_natveg", "0.")
        self.config.set(section, "pct_crop", "0.")
        read_cfg_subgrid(self.config, self.cfg_path)

    def test_read_subgrid_allspecial(self):
        """test a read of subgrid that's all special landunits"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "0. 0. 0.")
        self.config.set(section, "pct_lake", "25.")
        self.config.set(section, "pct_wetland", "35.")
        self.config.set(section, "pct_ocean", "0.")
        self.config.set(section, "pct_glacier", "40.")
        self.config.set(section, "pct_natveg", "0.")
        self.config.set(section, "pct_crop", "0.")
        read_cfg_subgrid(self.config, self.cfg_path)

    def test_read_subgrid_allurban(self):
        """test a read of subgrid that's all urban"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "100.0 0.0 0.0")
        self.config.set(section, "pct_lake", "0.")
        self.config.set(section, "pct_wetland", "0.")
        self.config.set(section, "pct_ocean", "0.")
        self.config.set(section, "pct_glacier", "0.")
        self.config.set(section, "pct_natveg", "0.")
        self.config.set(section, "pct_crop", "0.")
        read_cfg_subgrid(self.config, self.cfg_path)

    def test_read_subgrid_split_cropland(self):
        """
        test a read of subgrid that's 50/50 natural and
        cropland, with cropland split evenly among
        crop types
        """
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        self.config.set(section, "evenly_split_cropland", "True")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "0.0 0.0 0.0")
        self.config.set(section, "pct_lake", "0.")
        self.config.set(section, "pct_wetland", "0.")
        self.config.set(section, "pct_glacier", "0.")
        self.config.set(section, "pct_natveg", "50.")
        self.config.set(section, "pct_crop", "50.")
        read_cfg_subgrid(self.config, self.cfg_path)

    def test_read_var_list(self):
        """test a simple read of var_list"""
        read_cfg_var_list(self.config, idealized=True)

    def test_subgrid_outofrange(self):
        """test a read of subgrid that's out of range"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "101. 0. 0.")
        with self.assertRaisesRegex(SystemExit, "is out of range of 0 to 100 ="):
            read_cfg_subgrid(self.config, self.cfg_path)

    def test_subgrid_pct_urban_toosmall(self):
        """test a read of subgrid for PCT_URBAN that's an array too small"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "100. 0.")
        with self.assertRaisesRegex(
            SystemExit, "PCT_URBAN is not a list of the expected size of 3"
        ):
            read_cfg_subgrid(self.config, self.cfg_path)

    def test_subgrid_pct_urban_toobig(self):
        """test a read of subgrid for PCT_URBAN that's an array too big"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "100. 0. 0. 0.")
        with self.assertRaisesRegex(
            SystemExit, "PCT_URBAN is not a list of the expected size of 3"
        ):
            read_cfg_subgrid(self.config, self.cfg_path)

    def test_subgrid_pct_urban_singlevalue(self):
        """test a read of subgrid for PCT_URBAN that's a single value"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "100.")
        with self.assertRaisesRegex(
            SystemExit, "PCT_URBAN is not a list of the expected size of 3"
        ):
            read_cfg_subgrid(self.config, self.cfg_path)

    def test_subgrid_notsumtohundred(self):
        """test a read of subgrid that's doesn't sum to a hundred"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "pct_urban", "0. 0. 0.")
        self.config.set(section, "pct_lake", "0.")
        self.config.set(section, "pct_wetland", "0.")
        self.config.set(section, "pct_ocean", "0.")
        self.config.set(section, "pct_glacier", "0.")
        self.config.set(section, "pct_natveg", "0.")
        self.config.set(section, "pct_crop", "0.")
        with self.assertRaisesRegex(
            SystemExit, "PCT fractions in subgrid section do NOT sum to a hundred as they should"
        ):
            read_cfg_subgrid(self.config, self.cfg_path)

    def test_subgrid_badvar(self):
        """test a read of subgrid for a variable thats not in the list"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_subgrid_fractions"
        self.config.set(section, "badvariable", "100.")
        with self.assertRaisesRegex(SystemExit, "is not a valid variable name. Valid vars ="):
            read_cfg_subgrid(self.config, self.cfg_path)

    def test_varlist_varinidealized(self):
        """test a read of varlist for a variable thats in the idealized list,
        when idealized is on"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "True")
        section = "modify_fsurdat_variable_list"
        self.config.set(section, "PCT_SAND", "100.")
        with self.assertRaisesRegex(
            SystemExit,
            "is a special variable handled in the idealized section."
            + " This should NOT be handled in the variable list section. Special idealized vars =",
        ):
            read_cfg_var_list(self.config, idealized=True)

    def test_varlist_varinsubgrid(self):
        """test a read of varlist for a variable thats in the subgrid list"""
        section = "modify_fsurdat_basic_options"
        self.config.set(section, "idealized", "False")
        section = "modify_fsurdat_variable_list"
        self.config.set(section, "PCT_GLACIER", "100.")
        with self.assertRaisesRegex(
            SystemExit,
            "is a variable handled in the subgrid section."
            + " This should NOT be handled in the variable list section. Subgrid vars =",
        ):
            read_cfg_var_list(self.config, idealized=False)

    def test_varlist_monthlyvar(self):
        """test a read of varlist for a variable thats one of the monthly
        variables handled in the dom_pft section"""
        section = "modify_fsurdat_variable_list"
        self.config.set(section, "MONTHLY_LAI", "100.")
        with self.assertRaisesRegex(
            SystemExit,
            "is a variable handled as part of the dom_pft handling."
            + " This should NOT be handled in the variable list section."
            + " Monthly vars handled this way =",
        ):
            read_cfg_var_list(self.config, idealized=False)

    def test_subgrid_remove(self):
        """test a read of subgrid when it's section has been removed"""
        section = "modify_fsurdat_subgrid_fractions"
        self.config.remove_section(section)
        with self.assertRaisesRegex(SystemExit, "Config file does not have the expected section"):
            read_cfg_subgrid(self.config, self.cfg_path)

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
