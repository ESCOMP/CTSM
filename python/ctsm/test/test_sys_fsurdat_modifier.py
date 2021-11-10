#!/usr/bin/env python3

"""System tests for fsurdat_modifier

"""

import os

import unittest
import tempfile
import shutil
from subprocess import call

import ctsm.unit_testing

import xarray as xr


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

        # fsurdat_out that we expect to create
        prefix = 'test_'
        case = 'modify_minimal'
        suffix = '.nc'
        file_name_out = prefix + case + suffix

        if os.path.exists(file_name_out):
            os.remove(file_name_out)  # remove old file to write the new one

        # run the fsurdat_modifier tool
        prefix_1 = '../tools/modify_fsurdat/'  # location of tool
        tool = 'fsurdat_modifier'
        prefix_2 = 'ctsm/test/inputdata/'  # location of .cfg and fsurdat_in
        suffix = '.cfg'
        command = prefix_1 + tool + ' ' + prefix_2 + case + suffix
        call(command, shell=True)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # move the generated file to _tempdir
        file_name_final = os.path.join(self._tempdir, file_name_out)
        shutil.move(file_name_out, file_name_final)

        # TODO read fsurdat_in and fsurdat_out as xarray datasets and compare
        # TODO get filename_in from .cfg using argparse rather than hardwiring
        file_name_in = prefix_2 + \
            'surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc'
        fsurdat_in = xr.open_dataset(file_name_in)
        fsurdat_out = xr.open_dataset(file_name_final)

#   def test_allInfo(self):
#       """
#       This version specifies all possible information
#       """

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
