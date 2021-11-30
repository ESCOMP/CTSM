#!/usr/bin/env python3

"""Unit tests for utils
"""

import tempfile
import shutil
import unittest
import os

from ctsm import unit_testing
from ctsm.utils import fill_template_file, lon_range_0_to_360

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

class TestUtilsFillTemplateFile(unittest.TestCase):
    """Tests of utils: fill_template_file"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_fillTemplateFile_basic(self):
        """Basic test of fill_template_file"""
        template_path = os.path.join(self._testdir, 'template.txt')
        final_path = os.path.join(self._testdir, 'final.txt')
        template_contents = """\
Hello
$foo
Goodbye
$bar
"""
        with open(template_path, 'w') as f:
            f.write(template_contents)

        fillins = {'foo':'aardvark',
                   'bar':'zyzzyva'}
        fill_template_file(template_path, final_path, fillins)

        expected_final_text = """\
Hello
aardvark
Goodbye
zyzzyva
"""
        with open(final_path) as f:
            final_contents = f.read()

        self.assertEqual(final_contents, expected_final_text)


class TestUtilsLonRange0to360(unittest.TestCase):
    """Test of utils: lon_range_0_to_360"""

    def test_lonRange0To360_lonNegative(self):
        """
        Tests that negative inputs to lon_range_0_to_360 in the range [-180, 0)
        get 360 added to them
        """
        inputs = [-180, -0.001]
        for input in inputs:
            result = lon_range_0_to_360(input)
            self.assertEqual(result, input + 360)

    def test_lonRange0To360_lonNotNegative(self):
        """
        Tests that inputs to lon_range_0_to_360 of 0 to 360 remain unchanged
        """
        inputs = [0, 360]
        for input in inputs:
            result = lon_range_0_to_360(input)
            self.assertEqual(result, input)

    def test_lonRange0To360_outOfBounds(self):
        with self.assertRaisesRegex(SystemExit,
                                    "lon_in needs to be in the range 0 to 360"):
            _ = lon_range_0_to_360(361)


if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
