#!/usr/bin/env python3

"""
Unit tests for regrid_ggcmi_shdates subroutines:
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
from ctsm.crop_calendars.regrid_ggcmi_shdates import import_coord_1d
from ctsm.crop_calendars.regrid_ggcmi_shdates import import_coord_2d

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

# pylint: disable=protected-access


# Allow as many public methods as needed...
# pylint: disable=too-many-public-methods
# Allow all the instance attributes that we need
# pylint: disable=too-many-instance-attributes
class TestRegridGgcmiShdates(unittest.TestCase):
    # Tests the regrid_ggcmi_shdates subroutines

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
        my_data = xr.open_dataset(self._fsurdat_in)
        self.modify_fsurdat = ModifyFsurdat(
            my_data=my_data,
            lon_1=0.0,
            lon_2=360.0,
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
        shutil.rmtree(self._tempdir, ignore_errors=True)

    # def test_subgrid_and_idealized_fails(self):
    #     # Test that 
    #     section = "modify_fsurdat_basic_options"
    #     self.config.set(section, "idealized", "True")
    #     self.config.set(section, "include_nonveg", "False")
    #     self.config.set(section, "process_subgrid_section", "True")
    #     self.config.set(section, "dom_pft", "UNSET")
    #     with self.assertRaisesRegex(
    #         SystemExit,
    #         "idealized AND process_subgrid_section can NOT both be on, pick one or the other",
    #     ):
    #         read_cfg_option_control(self.modify_fsurdat, self.config, section, self.cfg_path)



if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
