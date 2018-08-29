#!/usr/bin/env python

"""Unit tests for run_sys_tests
"""

from __future__ import print_function
import unittest
import tempfile
import shutil
import os
from datetime import datetime
from six_additions import mock

from ctsm import add_cime_to_path #pylint:disable=unused-import
from ctsm import unit_testing
from ctsm.run_sys_tests import run_sys_tests
from ctsm.machine_defaults import MACHINE_DEFAULTS
from ctsm.machine import create_machine
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_FAKE

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
#pylint:disable=invalid-name

class TestRunSysTests(unittest.TestCase):
    """Tests of run_sys_tests"""

    _MACHINE_NAME = 'fake_machine'

    def setUp(self):
        self._original_wd = os.getcwd()
        self._curdir = tempfile.mkdtemp()
        os.chdir(self._curdir)
        self._scratch = os.path.join(self._curdir, 'scratch')
        os.makedirs(self._scratch)

    def tearDown(self):
        os.chdir(self._original_wd)
        shutil.rmtree(self._curdir, ignore_errors=True)

    def _make_machine(self):
        machine = create_machine(machine_name=self._MACHINE_NAME,
                                 defaults=MACHINE_DEFAULTS,
                                 job_launcher_type=JOB_LAUNCHER_FAKE,
                                 scratch_dir=self._scratch)
        return machine

    @staticmethod
    def _fake_now():
        return datetime(year=2001,
                        month=2,
                        day=3,
                        hour=4,
                        minute=5,
                        second=6)

    @staticmethod
    def _expected_testroot():
        """Returns an expected testroot based on values set in _fake_now and _MACHINE_NAME"""
        return 'tests_0203-0405fa'

    def test_testroot_setup(self):
        """Ensure that the appropriate test root directory is created and populated"""
        machine = self._make_machine()
        with mock.patch('ctsm.run_sys_tests.datetime') as mock_date:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(machine=machine, testlist=['foo'])

        expected_dir = os.path.join(self._scratch,
                                    self._expected_testroot())
        self.assertTrue(os.path.isdir(expected_dir))
        expected_link = os.path.join(self._curdir,
                                     self._expected_testroot())
        self.assertTrue(os.path.islink(expected_link))
        self.assertEqual(os.readlink(expected_link), expected_dir)
        # FIXME(wjs, 2018-08-28) Make sure it's populated with a cs.status file

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
