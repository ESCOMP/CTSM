#!/usr/bin/env python3

"""System tests for run_tower

"""

import os
import unittest
import tempfile
import shutil
import sys
import filecmp

from ctsm import unit_testing
from ctsm.remove_duplicate_tests import main
from ctsm.path_utils import path_to_ctsm_root

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysRemoveDuplicateTests(unittest.TestCase):
    """
    System tests for remove_duplicate_tests
    """

    def setUp(self):
        """
        Set up for these tests.
        """
        self._tempdir = tempfile.mkdtemp()
        self.output_file = os.path.join(self._tempdir, "test.xml")
        self.testdata_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "testinputs",
        )
        self.target_file = os.path.join(
                self.testdata_dir,
                "testlist_clm_target.xml"
            )
        self.script = os.path.join(
            path_to_ctsm_root(),
            "python",
            "ctsm",
            "remove_duplicate_tests.py",
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_extra_whitespace(self):
        """
        This test checks that whitespace is correctly removed
        """
        sys.argv = [
            self.script,
            "-i",
            os.path.join(
                self.testdata_dir,
                "testlist_clm_extra_whitespace.xml"
            ),
            "-o",
            self.output_file,
        ]
        main()
        assert(filecmp.cmp(
            self.output_file,
            self.target_file,
            shallow=True,
        ))

    def test_duplicate_machine(self):
        """
        This test checks that duplicate machines are removed
        """
        sys.argv = [
            self.script,
            "-i",
            os.path.join(
                self.testdata_dir,
                "testlist_clm_duplicate_machine.xml"
            ),
            "-o",
            self.output_file,
        ]
        main()
        assert(filecmp.cmp(
            self.output_file,
            self.target_file,
            shallow=True,
        ))

    def test_duplicate_test(self):
        """
        This test checks that duplicate tests are removed
        """
        sys.argv = [
            self.script,
            "-i",
            os.path.join(
                self.testdata_dir,
                "testlist_clm_duplicate_test.xml"
            ),
            "-o",
            self.output_file,
        ]
        main()
        print(self.output_file)
        print(self.target_file)
        assert(filecmp.cmp(
            self.output_file,
            self.target_file,
            shallow=False,
        ))

if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
