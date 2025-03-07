#!/usr/bin/env python3

"""
System tests for modify_singlept_site_neon.py
"""

import os
import unittest
import tempfile
import shutil
import sys

from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.site_and_regional.modify_singlept_site_neon import main

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysModifySingleptSiteNeon(unittest.TestCase):
    """System tests for modify_singlept_site_neon"""

    def setUp(self):
        """
        Make /_tempdir for use by these tests.
        Check tempdir for history files
        """
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._cfg_file_path = os.path.join(
            testinputs_path, "modify_singlept_site_neon_opt_sections.cfg"
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_modify_site(self):
        """
        Test modifying a singple point site.
        This test currently checks that the run fails due to dir structure

        TODO: The primary items to test here are the following:
        1) Fields are overwritten with site-specific data for neon sites
        2) Downloaded data is used in surface dataset
        3) Check specific fields listed in update_metadata for correct output
        4) Check that a netcdf with correct formatting is created
        """
        sys.argv = [
            "modify_singlept_site_neon",
            "--neon_site",
            path_to_ctsm_root() + "/ctsm/cime_config/usermods_dirs/clm/NEON/ABBY",
        ]
        # TODO: the above requires a full path instead of site name
        #       because of how run_tower is configured.
        # This needs to be fixed in run_tower.
        with self.assertRaises(SystemExit):
            print(
                """This should currently fail due to directory structure in run_tower
                   and the directory structure listed in sys.argv"""
            )
            main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
