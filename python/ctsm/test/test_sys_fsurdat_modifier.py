#!/usr/bin/env python3

"""System tests for fsurdat_modifier

"""

import os

import unittest
import tempfile
import shutil
from subprocess import call
from configparser import ConfigParser

import xarray as xr

import ctsm.unit_testing
from ctsm.utils import get_config_value, fill_template_file
from ctsm.path_utils import path_to_ctsm_root


class TestSysFsurdatModifier(unittest.TestCase):
    """System tests for fsurdat_modifier"""

    def setUp(self):
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_minimalInfo(self):
        """
        This test specifies a minimal amount of information
        """

        cfg_path, fsurdat_in, fsurdat_out = self.prepare_cfg_file(
            'modify_minimal')

        # run the fsurdat_modifier tool
        tool_path = os.path.join(path_to_ctsm_root(),
            'tools/modify_fsurdat/fsurdat_modifier')
        tool_command = tool_path + ' ' + cfg_path
        call(tool_command, shell=True)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        fsurdat_in_data = xr.open_dataset(fsurdat_in)
        fsurdat_out_data = xr.open_dataset(fsurdat_out)
        # assert that fsurdat_out equals fsurdat_in
        self.assertTrue(fsurdat_out_data.equals(fsurdat_in_data))


    def test_allInfo(self):
        """
        This version specifies all possible information
        """

        cfg_path, fsurdat_in, fsurdat_out = self.prepare_cfg_file(
            'modify_all')

        # run the fsurdat_modifier tool
        tool_path = os.path.join(path_to_ctsm_root(),
            'tools/modify_fsurdat/fsurdat_modifier')
        tool_command = tool_path + ' --debug ' + cfg_path
        call(tool_command, shell=True)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # compare fsurdat_out to fsurdat_in
        fsurdat_in_data = xr.open_dataset(fsurdat_in)
        fsurdat_out_data = xr.open_dataset(fsurdat_out)
        # assert that fsurdat_out does not equal fsurdat_in
        self.assertFalse(fsurdat_out_data.equals(fsurdat_in_data))

        # compare fsurdat_out to fsurdat_out_baseline
        fsurdat_out_baseline = fsurdat_in[:-3] + '_modified' + fsurdat_in[-3:]
        fsurdat_out_base_data = xr.open_dataset(fsurdat_out_baseline)
        # assert that fsurdat_out equals fsurdat_out_baseline
        self.assertTrue(fsurdat_out_data.equals(fsurdat_out_base_data))


    def prepare_cfg_file(self, case):
        """
        Enter fsurdat_in and fsurdat_out paths into the .cfg file
        using the fill_template_file function
        """
        testinputs_path = os.path.join(path_to_ctsm_root(),
                                       'python/ctsm/test/testinputs')
        fsurdat_in = os.path.join(testinputs_path,
            'surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc')
        file_name = case + '.nc'
        fsurdat_out = os.path.join(self._tempdir, file_name)
        substitutions = {'FSURDAT_IN': fsurdat_in,
                         'FSURDAT_OUT': fsurdat_out}
        file_name = case + '.cfg'
        cfg_template_path = os.path.join(testinputs_path, file_name)
        cfg_path = os.path.join(self._tempdir, file_name)
        fill_template_file(cfg_template_path, cfg_path, substitutions)

        return cfg_path, fsurdat_in, fsurdat_out


if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
