#!/usr/bin/env python

"""Unit tests for utils
"""

import tempfile
import shutil
import unittest
import os

from ctsm import unit_testing
from ctsm.utils import fill_template_file

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

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
