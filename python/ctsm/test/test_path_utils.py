#!/usr/bin/env python

"""Unit tests for path_utils
"""

import unittest
import tempfile
import shutil
import os

import six
from six_additions import mock
from ctsm import unit_testing
from ctsm import path_utils

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

class TestPathUtils(unittest.TestCase):
    """Tests of path_utils"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def _ctsm_path_in_cesm(self):
        """Returns the path to a ctsm directory nested inside a typical cesm
        directory structure, where self._testdir is the root of the cesm
        checkout
        """
        return os.path.join(self._testdir, 'components', 'clm')

    def _cime_path_in_cesm(self):
        """Returns the path to a cime directory nested inside a typical
        cesm directory structure, where self._testdir is the root of the
        cesm checkout
        """
        return os.path.join(self._testdir, 'cime')

    def _make_cesm_dirs(self):
        """Makes a directory structure for a typical CESM layout, where
        self._testdir is the root of the CESM checkout.

        This makes the ctsm directory and the cime directory.

        Returns a tuple, (ctsm_path, cime_path)
        """
        ctsm_path = self._ctsm_path_in_cesm()
        cime_path = self._cime_path_in_cesm()
        os.makedirs(ctsm_path)
        os.makedirs(cime_path)
        return (ctsm_path, cime_path)

    def test_pathToCime_standaloneOnlyWithCime(self):
        """Test path_to_cime with standalone_only, where cime is in the location
        it should be with a standalone checkout
        """
        ctsm_path = os.path.join(self._testdir, 'ctsm')
        actual_path_to_cime = os.path.join(ctsm_path, 'cime')
        os.makedirs(actual_path_to_cime)

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            path_to_cime = path_utils.path_to_cime(standalone_only=True)

        self.assertEqual(path_to_cime, actual_path_to_cime)

    def test_pathToCime_standaloneOnlyWithoutCime(self):
        """Test path_to_cime with standalone_only, where cime is missing"""
        ctsm_path = os.path.join(self._testdir, 'ctsm')
        os.makedirs(ctsm_path)

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            with six.assertRaisesRegex(self, RuntimeError, "Cannot find cime"):
                _ = path_utils.path_to_cime(standalone_only=True)

    def test_pathToCime_standaloneOnlyWithCimeInCesm(self):
        """Test path_to_cime with standalone_only, where cime is missing from
        the standalone structure, but cime is present in the CESM
        directory structure: should raise an exception rather than
        finding that cime
        """
        ctsm_path, _ = self._make_cesm_dirs()

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            with six.assertRaisesRegex(self, RuntimeError, "Cannot find cime"):
                _ = path_utils.path_to_cime(standalone_only=True)

    def test_pathToCime_cimeInCesm(self):
        """Test path_to_cime, where cime is not in the standalone directory but
        is present in the CESM directory structure
        """
        ctsm_path, actual_path_to_cime = self._make_cesm_dirs()

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            path_to_cime = path_utils.path_to_cime()

        self.assertEqual(path_to_cime, actual_path_to_cime)

    def test_pathToCime_notInCesmCheckout(self):
        """Test path_to_cime, where cime is not in the standalone directory, and
        we don't appear to be in a CESM checkout
        """
        ctsm_path = os.path.join(self._testdir, 'components', 'foo')
        os.makedirs(ctsm_path)
        os.makedirs(self._cime_path_in_cesm())

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            with six.assertRaisesRegex(self, RuntimeError,
                                       "Cannot find cime.*don't seem to be within a CESM checkout"):
                _ = path_utils.path_to_cime()

    def test_pathToCime_noCimeInCesm(self):
        """Test path_to_cime, where we appear to be within a CESM checkout, but
        there is no cime directory"""
        ctsm_path = self._ctsm_path_in_cesm()
        os.makedirs(ctsm_path)

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            with six.assertRaisesRegex(self, RuntimeError,
                                       "Cannot find cime.*or within CESM checkout"):
                _ = path_utils.path_to_cime()

    def test_pathToCime_cimeInStandaloneAndCesm(self):
        """Test path_to_cime, where there is a cime directory both in the
        standalone checkout and in the enclosing CESM checkout. Should
        give us the cime in the standalone checkout.
        """
        ctsm_path, _ = self._make_cesm_dirs()
        actual_path_to_cime = os.path.join(ctsm_path, 'cime')
        os.makedirs(actual_path_to_cime)

        with mock.patch('ctsm.path_utils.path_to_ctsm_root',
                        return_value=ctsm_path):
            path_to_cime = path_utils.path_to_cime()

        self.assertEqual(path_to_cime, actual_path_to_cime)

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
