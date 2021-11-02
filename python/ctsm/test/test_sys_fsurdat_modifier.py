#!/usr/bin/env python3

"""System tests for fsurdat_modifier

"""

import os

import unittest
import tempfile
import shutil
from subprocess import call

import ctsm.unit_testing


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

        prefix_1 = '../tools/modify_fsurdat/'
        prefix_2 = 'test_'
        case = 'modify_minimal'
        suffix = '.nc'
        file_name = prefix_1 + prefix_2 + case + suffix
        file_exists = os.path.exists(file_name)

        if file_exists:
            call('rm ' + file_name)

        command = 'fsurdat_modifier'
        suffix = '.cfg'
        call(prefix_1 + command + ' ' + prefix_1 + case + suffix)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions
        # below

        # ensure that ncdump of file created matches archived ncdump
        # of original version of that file
        call('ncdump test_modify_minimal.nc > test_modify_minimal.asc')
        call('diff test_modify_minimal.asc test_modify_minimal_orig.asc > dif.out')
        number_of_lines,file_name = call('wc -l dif.out')
        self.assertTrue(number_of_lines < 11)

    def test_allInfo(self):
        """
        This version specifies all possible information
        """
        prefix_1 = '../tools/modify_fsurdat/'
        prefix_2 = 'test_'
        case = 'modify'  # the only difference from previous test
        suffix = '.nc'
        file_name = prefix_1 + prefix_2 + case + suffix
        file_exists = os.path.exists(file_name)

        if file_exists:
            call('rm ' + file_name)

        command = 'fsurdat_modifier'
        suffix = '.cfg'
        call(prefix_1 + command + ' ' + prefix_1 + case + suffix)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions
        # below

        # ensure that ncdump of file created matches archived ncdump
        # of original version of that file
        call('ncdump test_modify.nc > test_modify.asc')
        call('diff test_modify.asc test_modify_orig.asc > dif.out')
        number_of_lines,file_name = call('wc -l dif.out')
        self.assertTrue(number_of_lines < 11)

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
