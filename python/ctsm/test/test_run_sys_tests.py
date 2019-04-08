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

from ctsm import add_cime_to_path # pylint: disable=unused-import
from ctsm import unit_testing
from ctsm.run_sys_tests import run_sys_tests
from ctsm.machine_defaults import MACHINE_DEFAULTS
from ctsm.machine import create_machine
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_FAKE

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

# Replace the slow _record_git_status with a fake that does nothing
@mock.patch('ctsm.run_sys_tests._record_git_status', mock.MagicMock(return_value=None))
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

    def _make_machine(self, account=None):
        machine = create_machine(machine_name=self._MACHINE_NAME,
                                 defaults=MACHINE_DEFAULTS,
                                 job_launcher_type=JOB_LAUNCHER_FAKE,
                                 scratch_dir=self._scratch,
                                 account=account)
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
        return '0203-040506fa'

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

        This test covers three things:

        (1) The use of a testlist argument

        (2) The standard arguments to create_test (the path to create_test, the arguments
        --test-id and --output-root, and the absence of --compare and --generate)

        (3) That a cs.status.fails file was created
        """
        machine = self._make_machine()
        with mock.patch('ctsm.run_sys_tests.datetime') as mock_date:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(machine=machine, cime_path=self._cime_path(),
                          testlist=['test1', 'test2'])

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        command = all_commands[0].cmd
        expected_create_test = os.path.join(self._cime_path(), 'scripts', 'create_test')
        six.assertRegex(self, command, r'^ *{}\s'.format(re.escape(expected_create_test)))
        six.assertRegex(self, command, r'--test-id +{}\s'.format(self._expected_testid()))
        expected_testroot_path = os.path.join(self._scratch, self._expected_testroot())
        six.assertRegex(self, command, r'--output-root +{}\s'.format(expected_testroot_path))
        six.assertRegex(self, command, r'test1 +test2(\s|$)')
        assertNotRegex(self, command, r'--compare\s')
        assertNotRegex(self, command, r'--generate\s')

        expected_cs_status = os.path.join(self._scratch,
                                          self._expected_testroot(),
                                          'cs.status.fails')
        self.assertTrue(os.path.isfile(expected_cs_status))

    def test_createTestCommand_testfileAndExtraArgs(self):
        """The correct create_test command should be run with a testfile and extra arguments

        This test covers three things:

        (1) The use of a testfile argument

        (2) The use of a bunch of optional arguments that are passed along to create_test

        (3) That a cs.status.fails file was created
        """
        machine = self._make_machine(account='myaccount')
        testroot_base = os.path.join(self._scratch, 'my', 'testroot')
        run_sys_tests(machine=machine, cime_path=self._cime_path(),
                      testfile='/path/to/testfile',
                      testid_base='mytestid',
                      testroot_base=testroot_base,
                      compare_name='mycompare',
                      generate_name='mygenerate',
                      baseline_root='myblroot',
                      walltime='3:45:67',
                      queue='runqueue',
                      extra_create_test_args='--some extra --createtest args')

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        command = all_commands[0].cmd
        six.assertRegex(self, command, r'--test-id +mytestid(\s|$)')
        expected_testroot = os.path.join(testroot_base, 'tests_mytestid')
        six.assertRegex(self, command, r'--output-root +{}(\s|$)'.format(expected_testroot))
        six.assertRegex(self, command, r'--testfile +/path/to/testfile(\s|$)')
        six.assertRegex(self, command, r'--compare +mycompare(\s|$)')
        six.assertRegex(self, command, r'--generate +mygenerate(\s|$)')
        six.assertRegex(self, command, r'--baseline-root +myblroot(\s|$)')
        six.assertRegex(self, command, r'--walltime +3:45:67(\s|$)')
        six.assertRegex(self, command, r'--queue +runqueue(\s|$)')
        six.assertRegex(self, command, r'--project +myaccount(\s|$)')
        six.assertRegex(self, command, r'--some +extra +--createtest +args(\s|$)')

        expected_cs_status = os.path.join(expected_testroot,
                                          'cs.status.fails')
        self.assertTrue(os.path.isfile(expected_cs_status))

    def test_createTestCommands_testsuite(self):
        """The correct create_test commands should be run with a test suite

        This tests that multiple create_test commands are run, one with each compiler in
        the given test suite for the given machine

        This test also checks the stdout and stderr files used for each command

        It also ensures that the cs.status.fails and cs.status files are created
        """
        machine = self._make_machine()
        with mock.patch('ctsm.run_sys_tests.datetime') as mock_date, \
             mock.patch('ctsm.run_sys_tests.get_tests_from_xml') as mock_get_tests:
            mock_date.now.side_effect = self._fake_now
            mock_get_tests.return_value = [{'compiler': 'intel'},
                                           {'compiler': 'pgi'},
                                           {'compiler': 'intel'}]
            run_sys_tests(machine=machine, cime_path=self._cime_path(),
                          suite_name='my_suite')

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 2)
        for command in all_commands:
            six.assertRegex(self, command.cmd,
                            r'--xml-category +{}(\s|$)'.format('my_suite'))
            six.assertRegex(self, command.cmd,
                            r'--xml-machine +{}(\s|$)'.format(self._MACHINE_NAME))

        six.assertRegex(self, all_commands[0].cmd, r'--xml-compiler +intel(\s|$)')
        six.assertRegex(self, all_commands[1].cmd, r'--xml-compiler +pgi(\s|$)')

        expected_testid1 = '{}_int'.format(self._expected_testid())
        expected_testid2 = '{}_pgi'.format(self._expected_testid())
        six.assertRegex(self, all_commands[0].cmd,
                        r'--test-id +{}(\s|$)'.format(expected_testid1))
        six.assertRegex(self, all_commands[1].cmd,
                        r'--test-id +{}(\s|$)'.format(expected_testid2))

        expected_testroot_path = os.path.join(self._scratch, self._expected_testroot())
        self.assertEqual(all_commands[0].out, os.path.join(expected_testroot_path,
                                                           'STDOUT.'+expected_testid1))
        self.assertEqual(all_commands[0].err, os.path.join(expected_testroot_path,
                                                           'STDERR.'+expected_testid1))
        self.assertEqual(all_commands[1].out, os.path.join(expected_testroot_path,
                                                           'STDOUT.'+expected_testid2))
        self.assertEqual(all_commands[1].err, os.path.join(expected_testroot_path,
                                                           'STDERR.'+expected_testid2))

        expected_cs_status = os.path.join(self._scratch,
                                          self._expected_testroot(),
                                          'cs.status')
        expected_cs_status = os.path.join(self._scratch,
                                          self._expected_testroot(),
                                          'cs.status.fails')
        self.assertTrue(os.path.isfile(expected_cs_status))

    def test_createTestCommands_testsuiteSpecifiedCompilers(self):
        """The correct commands should be run with a test suite where compilers are specified"""
        machine = self._make_machine()
        with mock.patch('ctsm.run_sys_tests.get_tests_from_xml') as mock_get_tests:
            # This value should be ignored; we just set it to make sure it's different
            # from the passed-in compiler list
            mock_get_tests.return_value = [{'compiler': 'intel'},
                                           {'compiler': 'pgi'},
                                           {'compiler': 'gnu'}]
            run_sys_tests(machine=machine, cime_path=self._cime_path(),
                          suite_name='my_suite',
                          suite_compilers=['comp1a', 'comp2b'])

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 2)
        six.assertRegex(self, all_commands[0].cmd, r'--xml-compiler +comp1a(\s|$)')
        six.assertRegex(self, all_commands[1].cmd, r'--xml-compiler +comp2b(\s|$)')

    def test_withDryRun_nothingDone(self):
        """With dry_run=True, no directories should be created, and no commands should be run"""
        machine = self._make_machine()
        run_sys_tests(machine=machine, cime_path=self._cime_path(), testlist=['foo'],
                      dry_run=True)
        self.assertEqual(os.listdir(self._scratch), [])
        self.assertEqual(machine.job_launcher.get_commands(), [])

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
