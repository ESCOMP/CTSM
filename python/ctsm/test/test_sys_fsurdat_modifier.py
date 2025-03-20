#!/usr/bin/env python3

"""System tests for fsurdat_modifier"""

import os
import re

import unittest
import tempfile
import shutil
import sys

import xarray as xr
import numpy as np

from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.modify_input_files.fsurdat_modifier import fsurdat_modifier
from ctsm.modify_input_files.fsurdat_modifier import fsurdat_modifier_arg_process

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysFsurdatModifier(unittest.TestCase):
    """System tests for fsurdat_modifier"""

    def setUp(self):
        """
        Obtain path to the existing:
        - modify_fsurdat_template.cfg file
        - /testinputs directory and fsurdat_in, located in /testinputs
        Make /_tempdir for use by these tests.
        Obtain path and names for the files being created in /_tempdir:
        - modify_fsurdat.cfg
        - fsurdat_out.nc
        """
        self._previous_dir = os.getcwd()
        self._cfg_template_path = os.path.join(
            path_to_ctsm_root(), "tools/modify_input_files/modify_fsurdat_template.cfg"
        )
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._testinputs_path = testinputs_path
        self._fsurdat_in = os.path.join(
            testinputs_path,
            "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self._cfg_file_path = os.path.join(self._tempdir, "modify_fsurdat.cfg")
        self._fsurdat_out = os.path.join(self._tempdir, "fsurdat_out.nc")

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_no_files_given_fail(self):
        """
        Test that if no input or output files are given that it will gracefully fail
        """
        self._cfg_file_path = os.path.join(
            self._testinputs_path, "modify_fsurdat_short_nofiles.cfg"
        )
        sys.argv = ["fsurdat_modifier", self._cfg_file_path]
        parser = fsurdat_modifier_arg_process()
        with self.assertRaisesRegex(SystemExit, "must contain item 'fsurdat_in'"):
            fsurdat_modifier(parser)

    def test_short_config(self):
        """
        Test that a short config file works
        """
        self._cfg_file_path = os.path.join(self._testinputs_path, "modify_fsurdat_short.cfg")
        sys.argv = ["fsurdat_modifier", self._cfg_file_path]
        parser = fsurdat_modifier_arg_process()
        fsurdat_out = os.path.join(
            self._testinputs_path, "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031.out.nc"
        )
        if os.path.exists(fsurdat_out):
            os.remove(fsurdat_out)
        fsurdat_modifier(parser)
        # Run it again with the overwrite option so that it will overwrite the file just created
        sys.argv = ["fsurdat_modifier", self._cfg_file_path, "--overwrite"]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)
        # Cleanup
        os.remove(fsurdat_out)

    def test_short_infile_both_cmdline_and_cfg(self):
        """
        Test that a graceful fail happens when the infile
        is given both in the command line and the config file
        """
        self._cfg_file_path = os.path.join(self._testinputs_path, "modify_fsurdat_short.cfg")
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
            "-i",
            "specify_fsurdat_in_on_cmd_line.nc",
        ]
        parser = fsurdat_modifier_arg_process()
        with self.assertRaisesRegex(
            SystemExit,
            "fsurdat_in is specified in both the command line and the config file, pick one",
        ):
            fsurdat_modifier(parser)

    def test_short_outfile_both_cmdline_and_cfg(self):
        """
        Test that a graceful fail happens when the outfile is given
        both in the command line and the config file
        """
        self._cfg_file_path = os.path.join(self._testinputs_path, "modify_fsurdat_short.cfg")
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
            "-o",
            "specify_fsurdat_out_on_cmd_line.nc",
        ]
        parser = fsurdat_modifier_arg_process()
        with self.assertRaisesRegex(
            SystemExit,
            "fsurdat_out is specified in both the command line and the config file, pick one",
        ):
            fsurdat_modifier(parser)

    def test_opt_sections(self):
        """
        Test that a simple file with the optional sections works
        """
        self._cfg_file_path = os.path.join(self._testinputs_path, "modify_fsurdat_opt_sections.cfg")
        outfile = os.path.join(
            self._tempdir,
            "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031_output_urban.nc",
        )
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
            "-i",
            os.path.join(
                self._testinputs_path, "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031.nc"
            ),
            "-o",
            outfile,
        ]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)
        # Read the resultant output file and make sure the fields are changed as expected
        fsurdat_out_data = xr.open_dataset(outfile)
        zero0d = np.zeros((5, 5))
        one0d = np.ones((5, 5))
        pct_urban = np.array(
            [
                [
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                ],
                [
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                ],
                [
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0],
                ],
            ]
        )
        lev2_two = np.empty((2, 3, 5, 5))
        lev2_two[0, :, :, :] = 200.0
        lev2_two[1, :, :, :] = 100.0
        lev2_ten = np.empty((10, 3, 5, 5))
        for x in range(10):
            lev2_ten[x, :, :, :] = float(x + 1)
        lev1 = np.array(
            [
                [
                    [200.0, 200.0, 200.0, 200.0, 200.0],
                    [200.0, 200.0, 200.0, 200.0, 200.0],
                    [200.0, 200.0, 200.0, 200.0, 200.0],
                    [200.0, 200.0, 200.0, 200.0, 200.0],
                    [200.0, 200.0, 200.0, 200.0, 200.0],
                ],
                [
                    [150.0, 150.0, 150.0, 150.0, 150.0],
                    [150.0, 150.0, 150.0, 150.0, 150.0],
                    [150.0, 150.0, 150.0, 150.0, 150.0],
                    [150.0, 150.0, 150.0, 150.0, 150.0],
                    [150.0, 150.0, 150.0, 150.0, 150.0],
                ],
                [
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                    [100.0, 100.0, 100.0, 100.0, 100.0],
                ],
            ]
        )
        np.testing.assert_array_equal(fsurdat_out_data.PCT_NATVEG, zero0d)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_CROP, zero0d)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_LAKE, zero0d)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_WETLAND, zero0d)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_OCEAN, zero0d)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_GLACIER, zero0d)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_URBAN, pct_urban)
        np.testing.assert_array_equal(fsurdat_out_data.LAKEDEPTH, one0d * 200.0)
        np.testing.assert_array_equal(fsurdat_out_data.T_BUILDING_MIN, lev1)
        np.testing.assert_array_equal(fsurdat_out_data.ALB_ROOF_DIR, lev2_two)
        np.testing.assert_array_equal(fsurdat_out_data.TK_ROOF, lev2_ten)

    def test_evenly_split_cropland(self):
        """
        Test that evenly splitting cropland works
        """
        self._create_config_file_evenlysplitcrop()
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
        ]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)
        # Read the resultant output file and make sure the fields are changed as expected
        fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        Ncrops = fsurdat_out_data.dims["cft"]
        pct_cft = np.full_like(fsurdat_out_data.PCT_CFT, 100 / Ncrops)
        np.testing.assert_array_equal(fsurdat_in_data.PCT_NATVEG, fsurdat_out_data.PCT_NATVEG)
        np.testing.assert_array_equal(fsurdat_in_data.PCT_CROP, fsurdat_out_data.PCT_CROP)
        np.testing.assert_array_equal(fsurdat_in_data.PCT_LAKE, fsurdat_out_data.PCT_LAKE)
        np.testing.assert_array_equal(fsurdat_in_data.PCT_WETLAND, fsurdat_out_data.PCT_WETLAND)
        np.testing.assert_array_equal(fsurdat_in_data.PCT_GLACIER, fsurdat_out_data.PCT_GLACIER)
        np.testing.assert_array_equal(fsurdat_in_data.PCT_URBAN, fsurdat_out_data.PCT_URBAN)
        np.testing.assert_array_equal(fsurdat_out_data.PCT_CFT, pct_cft)

    def test_1x1_mexicocity(self):
        """
        Test that the mexicocity file is handled correctly
        """
        self._cfg_file_path = os.path.join(
            self._testinputs_path, "modify_fsurdat_1x1mexicocity.cfg"
        )
        expectfile = os.path.join(
            self._testinputs_path,
            "surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103_modified.nc",
        )
        outfile = os.path.join(
            self._tempdir,
            "surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103_modified.nc",
        )
        infile = os.path.join(
            self._testinputs_path,
            "surfdata_1x1_mexicocityMEX_hist_16pfts_CMIP6_2000_c231103.nc",
        )
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
            "-i",
            infile,
            "-o",
            outfile,
        ]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)

        # Read the resultant output file and make sure the fields are changed as expected
        fsurdat_out_data = xr.open_dataset(outfile)
        fsurdat_inp_data = xr.open_dataset(infile)
        fsurdat_exp_data = xr.open_dataset(expectfile)

        self.assertFalse(fsurdat_out_data.equals(fsurdat_inp_data))
        # assert that fsurdat_out equals fsurdat_out_baseline
        self.assertTrue(fsurdat_out_data.equals(fsurdat_exp_data))

    def test_cfg_file_DNE_fail(self):
        """
        Test that if the config file does not exist that it gracefully fails
        """
        self._cfg_file_path = os.path.join(self._tempdir, "FILE_DOES_NOT_EXIST.cfg")
        sys.argv = ["fsurdat_modifier", self._cfg_file_path]
        with self.assertRaisesRegex(SystemExit, "Config file does NOT exist"):
            fsurdat_modifier_arg_process()

    def test_input_fsurdat_DNE_fail(self):
        """
        Test that if the input fsurdat  file does not exist that it gracefully fails
        """
        self._cfg_file_path = os.path.join(
            self._testinputs_path, "modify_fsurdat_short_nofiles.cfg"
        )
        sys.argv = ["fsurdat_modifier", self._cfg_file_path, "-i", "FILE_DOES_NOT_EXIST.nc"]
        parser = fsurdat_modifier_arg_process()
        with self.assertRaisesRegex(SystemExit, "Input fsurdat_in file does NOT exist"):
            fsurdat_modifier(parser)

    def test_output_fsurdat_EXISTS_fail(self):
        """
        Test that if the output fsurdat file does exist that it gracefully fails
        without --overwrite option
        """
        self._cfg_file_path = os.path.join(
            self._testinputs_path, "modify_fsurdat_short_nofiles.cfg"
        )
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
            "-i",
            self._cfg_file_path,
            "-o",
            self._cfg_file_path,
        ]
        parser = fsurdat_modifier_arg_process()
        with self.assertRaisesRegex(SystemExit, "Output file already exists"):
            fsurdat_modifier(parser)

    def test_cfg_file_empty_fail(self):
        """
        Test that if the config file is empty it gracefully fails
        """
        self._cfg_file_path = os.path.join(self._tempdir, "EMPTY_FILE.cfg")
        fil = open(self._cfg_file_path, "w")
        fil.close()
        sys.argv = ["fsurdat_modifier", self._cfg_file_path]
        parser = fsurdat_modifier_arg_process()
        with self.assertRaisesRegex(SystemExit, "Config file does not have the expected section"):
            fsurdat_modifier(parser)

    def test_minimalInfo(self):
        """
        This test specifies a minimal amount of information
        Create .cfg file, run the tool, compare fsurdat_in to fsurdat_out
        """

        self._create_config_file_minimal()

        # run the fsurdat_modifier tool
        sys.argv = ["fsurdat_modifier", self._cfg_file_path]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)
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
        sys.argv = ["fsurdat_modifier", self._cfg_file_path]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # compare fsurdat_out to fsurdat_in
        fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        # assert that fsurdat_out does not equal fsurdat_in
        self.assertFalse(fsurdat_out_data.equals(fsurdat_in_data))

        # compare fsurdat_out to fsurdat_out_baseline located in /testinputs
        fsurdat_out_baseline = self._fsurdat_in[:-3] + "_modified_with_crop" + self._fsurdat_in[-3:]
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
        sys.argv = [
            "fsurdat_modifier",
            self._cfg_file_path,
        ]
        parser = fsurdat_modifier_arg_process()
        fsurdat_modifier(parser)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # compare fsurdat_out to fsurdat_in
        fsurdat_in_data = xr.open_dataset(self._fsurdat_in)
        fsurdat_out_data = xr.open_dataset(self._fsurdat_out)
        # assert that fsurdat_out does not equal fsurdat_in
        self.assertFalse(fsurdat_out_data.equals(fsurdat_in_data))

        # compare fsurdat_out to fsurdat_out_baseline located in /testinputs
        fsurdat_out_baseline = self._fsurdat_in[:-3] + "_modified" + self._fsurdat_in[-3:]
        fsurdat_out_base_data = xr.open_dataset(fsurdat_out_baseline)
        # assert that fsurdat_out equals fsurdat_out_baseline
        self.assertTrue(fsurdat_out_data.equals(fsurdat_out_base_data))

    def _create_config_file_minimal(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *fsurdat_in *=", line):
                        line = f"fsurdat_in = {self._fsurdat_in}"
                    elif re.match(r" *fsurdat_out *=", line):
                        line = f"fsurdat_out = {self._fsurdat_out}"
                    cfg_out.write(line)

    def _create_config_file_evenlysplitcrop(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *evenly_split_cropland *=", line):
                        line = "evenly_split_cropland = True"
                    elif re.match(r" *fsurdat_in *=", line):
                        line = f"fsurdat_in = {self._fsurdat_in}"
                    elif re.match(r" *fsurdat_out *=", line):
                        line = f"fsurdat_out = {self._fsurdat_out}"
                    cfg_out.write(line)

    def _create_config_file_crop(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *fsurdat_in *=", line):
                        line = f"fsurdat_in = {self._fsurdat_in}"
                    elif re.match(r" *fsurdat_out *=", line):
                        line = f"fsurdat_out = {self._fsurdat_out}"
                    elif re.match(r" *lnd_lat_1 *=", line):
                        line = "lnd_lat_1 = -10\n"
                    elif re.match(r" *lnd_lat_2 *=", line):
                        line = "lnd_lat_2 = -7\n"
                    elif re.match(r" *lnd_lon_1 *=", line):
                        line = "lnd_lon_1 = 295\n"
                    elif re.match(r" *lnd_lon_2 *=", line):
                        line = "lnd_lon_2 = 300\n"
                    elif re.match(r" *dom_pft *=", line):
                        line = "dom_pft = 15"
                    elif re.match(r" *evenly_split_cropland *=", line):
                        line = "evenly_split_cropland = False"
                    elif re.match(r" *lai *=", line):
                        line = "lai = 0 1 2 3 4 5 5 4 3 2 1 0\n"
                    elif re.match(r" *sai *=", line):
                        line = "sai = 1 1 1 1 1 1 1 1 1 1 1 1\n"
                    elif re.match(r" *hgt_top *=", line):
                        line = "hgt_top = 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n"
                    elif re.match(r" *hgt_bot *=", line):
                        line = "hgt_bot = 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1\n"
                    cfg_out.write(line)

    def _create_config_file_complete(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *fsurdat_in *=", line):
                        line = f"fsurdat_in = {self._fsurdat_in}"
                    elif re.match(r" *fsurdat_out *=", line):
                        line = f"fsurdat_out = {self._fsurdat_out}"
                    elif re.match(r" *idealized *=", line):
                        line = "idealized = True"
                    elif re.match(r" *lnd_lat_1 *=", line):
                        line = "lnd_lat_1 = -10\n"
                    elif re.match(r" *lnd_lat_2 *=", line):
                        line = "lnd_lat_2 = -7\n"
                    elif re.match(r" *lnd_lon_1 *=", line):
                        line = "lnd_lon_1 = 295\n"
                    elif re.match(r" *lnd_lon_2 *=", line):
                        line = "lnd_lon_2 = 300\n"
                    elif re.match(r" *dom_pft *=", line):
                        line = "dom_pft = 1"
                    elif re.match(r" *evenly_split_cropland *=", line):
                        line = "evenly_split_cropland = False"
                    elif re.match(r" *lai *=", line):
                        line = "lai = 0 1 2 3 4 5 5 4 3 2 1 0\n"
                    elif re.match(r" *sai *=", line):
                        line = "sai = 1 1 1 1 1 1 1 1 1 1 1 1\n"
                    elif re.match(r" *hgt_top *=", line):
                        line = "hgt_top = 5 5 5 5 5 5 5 5 5 5 5 5\n"
                    elif re.match(r" *hgt_bot *=", line):
                        line = "hgt_bot = 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5\n"
                    elif re.match(r" *soil_color *=", line):
                        line = "soil_color = 5"
                    elif re.match(r" *std_elev *=", line):
                        line = "std_elev = 0.1"
                    elif re.match(r" *max_sat_area *=", line):
                        line = "max_sat_area = 0.2"
                    cfg_out.write(line)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
