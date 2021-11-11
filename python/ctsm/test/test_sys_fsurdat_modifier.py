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
from ctsm.utils import get_config_value


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

        # get fsurdat_in and fsurdat_out from the .cfg (config) file
        cfg_path = 'ctsm/test/testinputs/modify_minimal.cfg'
        config = ConfigParser()
        config.read(cfg_path)
        section = config.sections()[0]  # name of the first section
        fsurdat_in = get_config_value(config=config,
            section=section, item='fsurdat_in',
            file_path=cfg_path)
        fsurdat_out = get_config_value(config=config,
            section=section, item='fsurdat_out',
            file_path=cfg_path)

        if os.path.exists(fsurdat_out):
            os.remove(fsurdat_out)  # remove old file to write the new one

        # run the fsurdat_modifier tool
        tool = '../tools/modify_fsurdat/fsurdat_modifier'
        command = tool + ' ' + cfg_path
        call(command, shell=True)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # move the generated file to _tempdir
        fsurdat_out_final = os.path.join(self._tempdir, fsurdat_out)
        shutil.move(fsurdat_out, fsurdat_out_final)

        fsurdat_in_data = xr.open_dataset(fsurdat_in)
        fsurdat_out_data = xr.open_dataset(fsurdat_out_final)

        self.assertTrue(fsurdat_out_data.equals(fsurdat_in_data))


    def test_allInfo(self):
        """
        This version specifies all possible information
        """

        # get fsurdat_in and fsurdat_out from the .cfg (config) file
        cfg_path = 'ctsm/test/testinputs/modify_all.cfg'
        config = ConfigParser()
        config.read(cfg_path)
        section = config.sections()[0]  # name of the first section
        fsurdat_in = get_config_value(config=config,
            section=section, item='fsurdat_in',
            file_path=cfg_path)
        fsurdat_out = get_config_value(config=config,
            section=section, item='fsurdat_out',
            file_path=cfg_path)

        if os.path.exists(fsurdat_out):
            os.remove(fsurdat_out)  # remove old file to write the new one

        # run the fsurdat_modifier tool
        tool = '../tools/modify_fsurdat/fsurdat_modifier'
        command = tool + ' ' + cfg_path
        call(command, shell=True)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        # move the generated file to _tempdir
        fsurdat_out_final = os.path.join(self._tempdir, fsurdat_out)
        shutil.move(fsurdat_out, fsurdat_out_final)

        fsurdat_in_data = xr.open_dataset(fsurdat_in)
        fsurdat_out_data = xr.open_dataset(fsurdat_out_final)

        self.assertFalse(fsurdat_out_data.equals(fsurdat_in_data))
        # TODO read template file and self.assertTrue is
        # easier than needing to specify every variable...


if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
