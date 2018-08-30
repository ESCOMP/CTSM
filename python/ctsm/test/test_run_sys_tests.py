#!/usr/bin/env python

"""Unit tests for run_sys_tests
"""

from __future__ import print_function
import unittest
import tempfile
import shutil
import os
import re
from datetime import datetime

import six
from six_additions import mock, assertNotRegex

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

    def _cime_path(self):
        # For the sake of paths to scripts used in run_sys_tests: Pretend that cime exists
        # under the current directory, even though it doesn't
        return os.path.join(self._curdir, 'cime')

    @staticmethod
    def _expected_testid():
        """Returns an expected testid based on values set in _fake_now and _MACHINE_NAME"""
        return '0203-0405fa'

    def _expected_testroot(self):
        """Returns an expected testroot based on values set in _fake_now and _MACHINE_NAME

        This just returns the name of the testroot directory, not the full path"""
        return 'tests_{}'.format(self._expected_testid())

    def test_testroot_setup(self):
        """Ensure that the appropriate test root directory is created and populated"""
        machine = self._make_machine()
        with mock.patch('ctsm.run_sys_tests.datetime') as mock_date:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(machine=machine, cime_path=self._cime_path(),
                          testlist=['foo'])

        expected_dir = os.path.join(self._scratch,
                                    self._expected_testroot())
        self.assertTrue(os.path.isdir(expected_dir))
        expected_link = os.path.join(self._curdir,
                                     self._expected_testroot())
        self.assertTrue(os.path.islink(expected_link))
        self.assertEqual(os.readlink(expected_link), expected_dir)

    def test_createTestCommand_testnames(self):
        """The correct create_test command should be run when providing a list of test names

        This test does basic checking of the standard arguments to create_test
        """
        machine = self._make_machine()
        with mock.patch('ctsm.run_sys_tests.datetime') as mock_date:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(machine=machine, cime_path=self._cime_path(),
                          testlist=['test1', 'test2'])

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        command = all_commands[0]
        expected_create_test = os.path.join(self._cime_path(), 'scripts', 'create_test')
        six.assertRegex(self, command, r'^ *{}'.format(re.escape(expected_create_test)))
        six.assertRegex(self, command, r'--test-id +{}'.format(self._expected_testid()))
        expected_testroot_path = os.path.join(self._scratch, self._expected_testroot())
        six.assertRegex(self, command, r'--test-root +{}'.format(expected_testroot_path))
        six.assertRegex(self, command, r'test1 +test2 *$')
        assertNotRegex(self, command, r'--compare')
        assertNotRegex(self, command, r'--generate')

        # FIXME(wjs, 2018-08-29) Similar to the above test, but with testid_base, testroot_base specified and some optional args specified (compare_name, generate_name, baseline_root, walltime, queue, extra_create_test_args; also account?)

                          # compare_name='mycompare',
                          # generate_name='mygenerate',
                          # baseline_root='myblroot',
                          # walltime='3:45:67',
                          # queue='runqueue',
                          # extra_create_test_args='--some extra --createtest args'

        # cime/scripts/create_test --generate $newtag --compare $oldtag --baseline-root $baselineroot --test-id ${testid}_${compiler:0:1} --project ${account} --walltime ${walltime} --queue ${queue} ${extra_create_test_args} --test-root ${testroot} testname1 testname2


    # FIXME(wjs, 2018-08-29) Test with a test suite

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
