#!/usr/bin/env python3
"""
Unit tests for modify_singlept_site_neon

You can run this by:
    python -m unittest test_unit_modify_singlept_site_neon.py
"""

import unittest
import tempfile
import shutil
import configparser
import argparse
import os
import sys
from getpass import getuser
from datetime import date

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.modify_singlept_site_neon import (
    get_parser,
    get_neon,
    find_surffile,
    find_soil_structure,
    update_metadata,
    update_time_tag,
    sort_print_soil_layers,
    check_neon_time,
    download_file,
    fill_interpolate,
)
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name


class TestModifySingleptSiteNeon(unittest.TestCase):
    """
    Basic class for testing modify_singlept_site_neon.py.
    """

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        Check tempdir for history files
        """
        self._tempdir = tempfile.mkdtemp()

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_get_neon(self):
        """
        Test to see if neon data for valid site name is found
        """
        site_name = "ABBY"
        neon_dir = self._tempdir
        file = get_neon(neon_dir, site_name)
        self.assertEqual(file.split("/")[-1][:4], "ABBY", "CSV file did not download as expected")

    def test_get_neon_false_site(self):
        """
        Test to see if neon data for invalid site name is found
        """
        site_name = "ABY"
        neon_dir = self._tempdir
        with self.assertRaises(SystemExit):
            get_neon(neon_dir, site_name)

    def test_find_surffile(self):
        """
        Test that surface file does not exist in tempdir and raises system exit error
        """
        surf_dir = self._tempdir
        site_name = "BART"
        pft_16 = True
        with self.assertRaises(SystemExit):
            find_surffile(surf_dir, site_name, pft_16)
        # TODO: it would be useful to also include an example where the surface file is actually found or correctly searched for
        # Is there a default surface data file directory we can point to?

    # def test_find_soil_structure(self):
    #    """
    #    Test
    #    """
    #    find_soil_structure()
    #    TODO: need file to use in test

    # def test_update_metadata(self):
    #    TODO: do we have a surface data file that we can put in tempdir to update?
    #    update_metadata(nc, surf_file, neon_file, zb_flag)

    def test_update_time_tag(self):
        """
        Test that file ending is updated
        """
        self.assertEqual(
            update_time_tag("test_YYMMDD.nc")[-9:-3],
            date.today().strftime("%y%m%d"),
            "File ending not as expected",
        )

    def test_check_neon_time(self):
        """
        Test that dictionary containing last modified information is correctly downloaded
        """
        last_abby_download = check_neon_time()[
            "https://storage.neonscience.org/neon-ncar/NEON/surf_files/v1/ABBY_surfaceData.csv"
        ]
        self.assertEqual(
            len(last_abby_download),
            19,
            "last ABBY download has unexpected date format or does not exist",
        )
        # TODO: this checks that data is not pulled from before 2021; we may want to update this occassionally or find another way to check data is actually newest?
        self.assertGreater(
            int(last_abby_download[:4]), 2021, "ABBY download is older than expected"
        )

    # TODO: test sort_print_soil_layers, fill_interpolate, find_soil_structure, and update_metadata


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
