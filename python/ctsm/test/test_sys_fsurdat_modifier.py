#!/usr/bin/env python3

"""System tests for fsurdat_modifier

"""

import os
import re

import unittest
import tempfile
import shutil

import xarray as xr

from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.modify_fsurdat.fsurdat_modifier import fsurdat_modifier

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name

class TestSysFsurdatModifier(unittest.TestCase):
    """System tests for fsurdat_modifier"""

    def setUp(self):
        """
        Obtain path to the existing:
        - modify_template.cfg file
        - /testinputs directory and fsurdat_in, located in /testinputs
        Make /_tempdir for use by these tests.
        Obtain path and names for the files being created in /_tempdir:
        - modify_fsurdat.cfg
        - fsurdat_out.nc
        """
        self._cfg_template_path = os.path.join(path_to_ctsm_root(),
            'tools/modify_fsurdat/modify_template.cfg')
        testinputs_path = os.path.join(path_to_ctsm_root(),
            'python/ctsm/test/testinputs')
        self._fsurdat_in = os.path.join(testinputs_path,
            'surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc')
        self._tempdir = tempfile.mkdtemp()
        self._cfg_file_path = os.path.join(self._tempdir, 'modify_fsurdat.cfg')
        self._fsurdat_out = os.path.join(self._tempdir, 'fsurdat_out.nc')

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_minimalInfo(self):
        """
        This test specifies a minimal amount of information
        Create .cfg file, run the tool, compare fsurdat_in to fsurdat_out
        """

        self._create_config_file_minimal()

        # run the fsurdat_modifier tool
        fsurdat_modifier(self._cfg_file_path)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        # assert that fsurdat_out equals fsurdat_in
        self.assertTrue(fsurdat_out_data.equals(fsurdat_in_data))


    def test_crop(self):
        """
        This version replaces the vegetation with a crop
        Create .cfg file, run the tool, compare fsurdat_in to fsurdat_out
        """

        self._create_config_file_crop()

        # run the fsurdat_modifier tool
        fsurdat_modifier(self._cfg_file_path)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # compare fsurdat_out to fsurdat_in
        fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        # assert that fsurdat_out does not equal fsurdat_in
        self.assertFalse(fsurdat_out_data.equals(fsurdat_in_data))

        # compare fsurdat_out to fsurdat_out_baseline located in /testinputs
        fsurdat_out_baseline = self._fsurdat_in[:-3] + '_modified_with_crop' + \
                               self._fsurdat_in[-3:]
        fsurdat_out_base_data = xr.open_dataset(fsurdat_out_baseline)
        # assert that fsurdat_out equals fsurdat_out_baseline
        self.assertTrue(fsurdat_out_data.equals(fsurdat_out_base_data))


    def test_allInfo(self):
        """
        This version specifies all possible information
        Create .cfg file, run the tool, compare fsurdat_in to fsurdat_out
        """

        self._create_config_file_complete()

        # run the fsurdat_modifier tool
        fsurdat_modifier(self._cfg_file_path)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # compare fsurdat_out to fsurdat_in
        fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        # assert that fsurdat_out does not equal fsurdat_in
        self.assertFalse(fsurdat_out_data.equals(fsurdat_in_data))

        # compare fsurdat_out to fsurdat_out_baseline located in /testinputs
        fsurdat_out_baseline = self._fsurdat_in[:-3] + '_modified' + \
                               self._fsurdat_in[-3:]
        fsurdat_out_base_data = xr.open_dataset(fsurdat_out_baseline)
        # assert that fsurdat_out equals fsurdat_out_baseline
        self.assertTrue(fsurdat_out_data.equals(fsurdat_out_base_data))


    def _create_config_file_minimal(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open (self._cfg_file_path, 'w', encoding='utf-8') as cfg_out:
            with open (self._cfg_template_path, 'r', encoding='utf-8') as cfg_in:
                for line in cfg_in:
                    if re.match(r' *fsurdat_in *=', line):
                        line = f'fsurdat_in = {self._fsurdat_in}'
                    elif re.match(r' *fsurdat_out *=', line):
                        line = f'fsurdat_out = {self._fsurdat_out}'
                    cfg_out.write(line)


    def _create_config_file_crop(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open (self._cfg_file_path, 'w', encoding='utf-8') as cfg_out:
            with open (self._cfg_template_path, 'r', encoding='utf-8') as cfg_in:
                for line in cfg_in:
                    if re.match(r' *fsurdat_in *=', line):
                        line = f'fsurdat_in = {self._fsurdat_in}'
                    elif re.match(r' *fsurdat_out *=', line):
                        line = f'fsurdat_out = {self._fsurdat_out}'
                    elif re.match(r' *lnd_lat_1 *=', line):
                        line = 'lnd_lat_1 = -10\n'
                    elif re.match(r' *lnd_lat_2 *=', line):
                        line = 'lnd_lat_2 = -7\n'
                    elif re.match(r' *lnd_lon_1 *=', line):
                        line = 'lnd_lon_1 = 295\n'
                    elif re.match(r' *lnd_lon_2 *=', line):
                        line = 'lnd_lon_2 = 300\n'
                    elif re.match(r' *dom_plant *=', line):
                        line = 'dom_plant = 15'
                    elif re.match(r' *lai *=', line):
                        line = 'lai = 0 1 2 3 4 5 5 4 3 2 1 0\n'
                    elif re.match(r' *sai *=', line):
                        line = 'sai = 1 1 1 1 1 1 1 1 1 1 1 1\n'
                    elif re.match(r' *hgt_top *=', line):
                        line = 'hgt_top = 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n'
                    elif re.match(r' *hgt_bot *=', line):
                        line = 'hgt_bot = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1\n'
                    cfg_out.write(line)


    def _create_config_file_complete(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open (self._cfg_file_path, 'w', encoding='utf-8') as cfg_out:
            with open (self._cfg_template_path, 'r', encoding='utf-8') as cfg_in:
                for line in cfg_in:
                    if re.match(r' *fsurdat_in *=', line):
                        line = f'fsurdat_in = {self._fsurdat_in}'
                    elif re.match(r' *fsurdat_out *=', line):
                        line = f'fsurdat_out = {self._fsurdat_out}'
                    elif re.match(r' *idealized *=', line):
                        line = 'idealized = True'
                    elif re.match(r' *lnd_lat_1 *=', line):
                        line = 'lnd_lat_1 = -10\n'
                    elif re.match(r' *lnd_lat_2 *=', line):
                        line = 'lnd_lat_2 = -7\n'
                    elif re.match(r' *lnd_lon_1 *=', line):
                        line = 'lnd_lon_1 = 295\n'
                    elif re.match(r' *lnd_lon_2 *=', line):
                        line = 'lnd_lon_2 = 300\n'
                    elif re.match(r' *dom_plant *=', line):
                        line = 'dom_plant = 1'
                    elif re.match(r' *lai *=', line):
                        line = 'lai = 0 1 2 3 4 5 5 4 3 2 1 0\n'
                    elif re.match(r' *sai *=', line):
                        line = 'sai = 1 1 1 1 1 1 1 1 1 1 1 1\n'
                    elif re.match(r' *hgt_top *=', line):
                        line = 'hgt_top = 5 5 5 5 5 5 5 5 5 5 5 5\n'
                    elif re.match(r' *hgt_bot *=', line):
                        line = 'hgt_bot = 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n'
                    elif re.match(r' *soil_color *=', line):
                        line = 'soil_color = 5'
                    elif re.match(r' *std_elev *=', line):
                        line = 'std_elev = 0.1'
                    elif re.match(r' *max_sat_area *=', line):
                        line = 'max_sat_area = 0.2'
                    cfg_out.write(line)


if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
