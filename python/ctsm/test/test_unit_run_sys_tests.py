#!/usr/bin/env python3

"""Unit tests for run_sys_tests"""

from __future__ import print_function
import argparse
import unittest
import tempfile
import shutil
import os
import re
from datetime import datetime
from unittest import mock

from ctsm import add_cime_to_path  # pylint: disable=unused-import
from ctsm import unit_testing
from ctsm.run_sys_tests import run_sys_tests, _get_testmod_list, _check_arg_validity
from ctsm.machine_defaults import MACHINE_DEFAULTS
from ctsm.machine import create_machine
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_FAKE, JOB_LAUNCHER_QSUB

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


# Replace the slow _record_git_status with a fake that does nothing
@mock.patch("ctsm.run_sys_tests._record_git_status", mock.MagicMock(return_value=None))
class TestRunSysTests(unittest.TestCase):
    """Tests of run_sys_tests"""

    _MACHINE_NAME = "fake_machine"

    def setUp(self):
        self._original_wd = os.getcwd()
        self._curdir = tempfile.mkdtemp()
        os.chdir(self._curdir)
        self._scratch = os.path.join(self._curdir, "scratch")
        os.makedirs(self._scratch)

    def tearDown(self):
        os.chdir(self._original_wd)
        shutil.rmtree(self._curdir, ignore_errors=True)

    def _make_machine(self, account=None):
        machine = create_machine(
            machine_name=self._MACHINE_NAME,
            defaults=MACHINE_DEFAULTS,
            job_launcher_type=JOB_LAUNCHER_FAKE,
            scratch_dir=self._scratch,
            account=account,
        )
        return machine

    @staticmethod
    def _fake_now():
        return datetime(year=2001, month=2, day=3, hour=4, minute=5, second=6)

    def _cime_path(self):
        # For the sake of paths to scripts used in run_sys_tests: Pretend that cime exists
        # under the current directory, even though it doesn't
        return os.path.join(self._curdir, "cime")

    @staticmethod
    def _expected_testid():
        """Returns an expected testid based on values set in _fake_now and _MACHINE_NAME"""
        return "0203-040506fa"

    def _expected_testroot(self):
        """Returns an expected testroot based on values set in _fake_now and _MACHINE_NAME

        This just returns the name of the testroot directory, not the full path"""
        return "tests_{}".format(self._expected_testid())

    def test_testroot_setup(self):
        """Ensure that the appropriate test root directory is created and populated"""
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.datetime") as mock_date:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(machine=machine, cime_path=self._cime_path(), testlist=["foo"])

        expected_dir = os.path.join(self._scratch, self._expected_testroot())
        self.assertTrue(os.path.isdir(expected_dir))
        expected_link = os.path.join(self._curdir, self._expected_testroot())
        self.assertTrue(os.path.islink(expected_link))
        self.assertEqual(os.readlink(expected_link), expected_dir)

    def test_createTestCommand_testnames(self):
        """The correct create_test command should be run when providing a list of test names

        This test covers three things:

        (1) The use of a testlist argument

        (2) The standard arguments to create_test (the path to create_test, the arguments
        --test-id, --output-root and --retry, the absence of --compare and --generate, and
        (on this unknown machine) the absence of --baseline-root)

        (3) That a cs.status.fails file was created
        """
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.datetime") as mock_date:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                testlist=["test1", "test2"],
            )

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        command = all_commands[0].cmd
        expected_create_test = os.path.join(self._cime_path(), "scripts", "create_test")
        self.assertRegex(command, r"^ *{}\s".format(re.escape(expected_create_test)))
        self.assertRegex(command, r"--test-id +{}\s".format(self._expected_testid()))
        expected_testroot_path = os.path.join(self._scratch, self._expected_testroot())
        self.assertRegex(command, r"--output-root +{}\s".format(expected_testroot_path))
        self.assertRegex(command, r"--retry +0(\s|$)")
        self.assertRegex(command, r"test1 +test2(\s|$)")
        self.assertNotRegex(command, r"--compare\s")
        self.assertNotRegex(command, r"--generate\s")
        self.assertNotRegex(command, r"--baseline-root\s")
        # In the machine object for this test, create_test_queue will be 'unspecified';
        # verify that this results in there being no '--queue' argument:
        self.assertNotRegex(command, r"--queue\s")

        expected_cs_status = os.path.join(
            self._scratch, self._expected_testroot(), "cs.status.fails"
        )
        self.assertTrue(os.path.isfile(expected_cs_status))

    def test_createTestCommand_testfileAndExtraArgs(self):
        """The correct create_test command should be run with a testfile and extra arguments

        This test covers three things:

        (1) The use of a testfile argument

        (2) The use of a bunch of optional arguments that are passed along to create_test

        (3) That a cs.status.fails file was created
        """
        machine = self._make_machine(account="myaccount")
        testroot_base = os.path.join(self._scratch, "my", "testroot")
        run_sys_tests(
            machine=machine,
            cime_path=self._cime_path(),
            testfile="/path/to/testfile",
            testid_base="mytestid",
            testroot_base=testroot_base,
            compare_name="mycompare",
            generate_name="mygenerate",
            baseline_root="myblroot",
            walltime="3:45:67",
            queue="runqueue",
            retry=5,
            extra_create_test_args="--some extra --createtest args",
        )

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        command = all_commands[0].cmd
        self.assertRegex(command, r"--test-id +mytestid(\s|$)")
        expected_testroot = os.path.join(testroot_base, "tests_mytestid")
        self.assertRegex(command, r"--output-root +{}(\s|$)".format(expected_testroot))
        self.assertRegex(command, r"--testfile +/path/to/testfile(\s|$)")
        self.assertRegex(command, r"--compare +mycompare(\s|$)")
        self.assertRegex(command, r"--generate +mygenerate(\s|$)")
        self.assertRegex(command, r"--baseline-root +myblroot(\s|$)")
        self.assertRegex(command, r"--walltime +3:45:67(\s|$)")
        self.assertRegex(command, r"--queue +runqueue(\s|$)")
        self.assertRegex(command, r"--project +myaccount(\s|$)")
        self.assertRegex(command, r"--retry +5(\s|$)")
        self.assertRegex(command, r"--some +extra +--createtest +args(\s|$)")

        expected_cs_status = os.path.join(expected_testroot, "cs.status.fails")
        self.assertTrue(os.path.isfile(expected_cs_status))

    def test_createTestCommands_testsuite(self):
        """The correct create_test commands should be run with a test suite

        This tests that multiple create_test commands are run, one with each compiler in
        the given test suite for the given machine

        This test also checks the stdout and stderr files used for each command

        It also ensures that the cs.status.fails and cs.status files are created
        """
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.datetime") as mock_date, mock.patch(
            "ctsm.run_sys_tests.get_tests_from_xml"
        ) as mock_get_tests:
            mock_date.now.side_effect = self._fake_now
            mock_get_tests.return_value = [
                {"compiler": "intel"},
                {"compiler": "pgi"},
                {"compiler": "intel"},
            ]
            run_sys_tests(machine=machine, cime_path=self._cime_path(), suite_name="my_suite")

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 2)
        for command in all_commands:
            self.assertRegex(command.cmd, r"--xml-category +{}(\s|$)".format("my_suite"))
            self.assertRegex(command.cmd, r"--xml-machine +{}(\s|$)".format(self._MACHINE_NAME))

        self.assertRegex(all_commands[0].cmd, r"--xml-compiler +intel(\s|$)")
        self.assertRegex(all_commands[1].cmd, r"--xml-compiler +pgi(\s|$)")

        expected_testid1 = "{}_int".format(self._expected_testid())
        expected_testid2 = "{}_pgi".format(self._expected_testid())
        self.assertRegex(all_commands[0].cmd, r"--test-id +{}(\s|$)".format(expected_testid1))
        self.assertRegex(all_commands[1].cmd, r"--test-id +{}(\s|$)".format(expected_testid2))

        expected_testroot_path = os.path.join(self._scratch, self._expected_testroot())
        self.assertEqual(
            all_commands[0].out,
            os.path.join(expected_testroot_path, "STDOUT." + expected_testid1),
        )
        self.assertEqual(
            all_commands[0].err,
            os.path.join(expected_testroot_path, "STDERR." + expected_testid1),
        )
        self.assertEqual(
            all_commands[1].out,
            os.path.join(expected_testroot_path, "STDOUT." + expected_testid2),
        )
        self.assertEqual(
            all_commands[1].err,
            os.path.join(expected_testroot_path, "STDERR." + expected_testid2),
        )

        expected_cs_status = os.path.join(self._scratch, self._expected_testroot(), "cs.status")
        expected_cs_status = os.path.join(
            self._scratch, self._expected_testroot(), "cs.status.fails"
        )
        self.assertTrue(os.path.isfile(expected_cs_status))

    def test_createTestCommands_testsuiteSpecifiedCompilers(self):
        """The correct commands should be run with a test suite where compilers are specified"""
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.get_tests_from_xml") as mock_get_tests:
            # This value should be ignored; we just set it to make sure it's different
            # from the passed-in compiler list
            mock_get_tests.return_value = [
                {"compiler": "intel"},
                {"compiler": "pgi"},
                {"compiler": "gnu"},
            ]
            run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                suite_name="my_suite",
                suite_compilers=["comp1a", "comp2b"],
            )

        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 2)
        self.assertRegex(all_commands[0].cmd, r"--xml-compiler +comp1a(\s|$)")
        self.assertRegex(all_commands[1].cmd, r"--xml-compiler +comp2b(\s|$)")

    def test_createTestCommands_testsuiteXmlMachine(self):
        """With xml_machine given, the suite is queried for THAT machine

        The machine we are actually running on is unchanged: the testid still
        carries this machine's prefix. This is what lets the container run
        derecho's testlist without claiming to be derecho.
        """
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.datetime") as mock_date, mock.patch(
            "ctsm.run_sys_tests.get_tests_from_xml"
        ) as mock_get_tests:
            mock_date.now.side_effect = self._fake_now
            mock_get_tests.return_value = [{"compiler": "gnu"}]
            run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                suite_name="my_suite",
                xml_machine="derecho",
            )

        mock_get_tests.assert_called_once_with(xml_machine="derecho", xml_category="my_suite")
        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        self.assertRegex(all_commands[0].cmd, r"--xml-machine +derecho(\s|$)")
        self.assertNotRegex(
            all_commands[0].cmd, r"--xml-machine +{}(\s|$)".format(self._MACHINE_NAME)
        )
        # The machine we are running on is unaffected: testid prefix is still "fa"
        self.assertRegex(
            all_commands[0].cmd,
            r"--test-id +{}_gnu(\s|$)".format(self._expected_testid()),
        )

    def test_createTestCommands_testsuiteXmlMachineAndCompilers(self):
        """With --xml-machine and an explicit --suite-compiler together

        This is exactly the combination the container wrapper uses: it already
        knows which compiler to run (so no auto-discovery is needed), but it
        still needs to query a different machine's testlist entries via
        xml_machine. get_tests_from_xml must NOT be called in this case, since
        suite_compilers is already given.
        """
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.datetime") as mock_date, mock.patch(
            "ctsm.run_sys_tests.get_tests_from_xml"
        ) as mock_get_tests:
            mock_date.now.side_effect = self._fake_now
            run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                suite_name="my_suite",
                suite_compilers=["gnu"],
                xml_machine="derecho",
            )

        mock_get_tests.assert_not_called()
        all_commands = machine.job_launcher.get_commands()
        self.assertEqual(len(all_commands), 1)
        command = all_commands[0].cmd
        self.assertRegex(command, r"--xml-machine +derecho(\s|$)")
        self.assertRegex(command, r"--xml-category +my_suite(\s|$)")
        self.assertRegex(command, r"--xml-compiler +gnu(\s|$)")

    def test_withDryRun_nothingDone(self):
        """With dry_run=True, no directories should be created, and no commands should be run"""
        machine = self._make_machine()
        run_sys_tests(machine=machine, cime_path=self._cime_path(), testlist=["foo"], dry_run=True)
        self.assertEqual(os.listdir(self._scratch), [])
        self.assertEqual(machine.job_launcher.get_commands(), [])

    def test_wait_returnsCreateTestStatus(self):
        """With wait=True, run_sys_tests returns create_test's exit status

        A failing test therefore surfaces as a nonzero exit, which is what CI
        needs and what makes the container wrapper report failures.
        """
        machine = self._make_machine()
        machine.job_launcher.set_return_code(3)
        return_code = run_sys_tests(
            machine=machine, cime_path=self._cime_path(), testlist=["foo"], wait=True
        )
        self.assertEqual(return_code, 3)

    def test_noWait_returnsZero(self):
        """Without wait, run_sys_tests returns 0 even if create_test would fail

        This is the pre-existing behavior and must not change: run_sys_tests
        normally returns as soon as the job is launched.
        """
        machine = self._make_machine()
        machine.job_launcher.set_return_code(3)
        return_code = run_sys_tests(machine=machine, cime_path=self._cime_path(), testlist=["foo"])
        self.assertEqual(return_code, 0)

    def test_wait_withDryRun_returnsZero(self):
        """With dry_run there is no process to wait for, so the status is 0"""
        machine = self._make_machine()
        machine.job_launcher.set_return_code(3)
        return_code = run_sys_tests(
            machine=machine,
            cime_path=self._cime_path(),
            testlist=["foo"],
            wait=True,
            dry_run=True,
        )
        self.assertEqual(return_code, 0)

    def test_wait_multiCompiler_firstFails_returnsNonzero(self):
        """With wait=True and a multi-compiler suite, a failure in the FIRST-launched
        compiler's create_test must still produce a nonzero exit status, even though the
        last-launched compiler succeeds.

        This is exactly the scenario --wait exists to protect against: inside a
        container, run_sys_tests exiting tears down the PID namespace and kills any
        create_test processes still running, so we must wait for (and check the status
        of) every launched process, not just the last one. Under the old
        wait-for-last-only behavior, this would incorrectly return 0.
        """
        machine = self._make_machine()
        with mock.patch("ctsm.run_sys_tests.get_tests_from_xml") as mock_get_tests:
            mock_get_tests.return_value = [
                {"compiler": "intel"},
                {"compiler": "pgi"},
            ]
            # Compilers are launched in sorted order (intel, then pgi); the
            # first-launched (intel) fails, the last-launched (pgi) succeeds.
            machine.job_launcher.set_return_codes([5, 0])
            return_code = run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                suite_name="my_suite",
                wait=True,
            )
        self.assertEqual(return_code, 5)
        # The fake launcher's canned return codes are independent of what was actually
        # launched, so assert the launch count too: without this, the test would pass
        # even if the multi-compiler loop stopped launching a create_test per compiler.
        self.assertEqual(len(machine.job_launcher.get_commands()), 2)

    def test_wait_withoutNoBatchLauncher_raisesBeforeLaunching(self):
        """--wait with a launcher that can't wait (e.g. qsub) must fail fast

        This must raise before any create_test is launched, rather than after
        submitting jobs via qsub (which the test environment can't actually
        run) and only then failing.
        """
        machine = create_machine(
            machine_name=self._MACHINE_NAME,
            defaults=MACHINE_DEFAULTS,
            job_launcher_type=JOB_LAUNCHER_QSUB,
            scratch_dir=self._scratch,
            account="fakeaccount",
            job_launcher_queue="somequeue",
            job_launcher_walltime="00:10:00",
        )
        with self.assertRaises(RuntimeError):
            run_sys_tests(
                machine=machine,
                cime_path=self._cime_path(),
                testlist=["foo"],
                wait=True,
            )
        # Nothing should have been launched
        self.assertEqual(os.listdir(self._scratch), [])

    def test_getTestmodList_suite(self):
        """Ensure that _get_testmod_list() works correctly with suite-style input"""
        testmod_list_input = [
            "clm/default",
            "clm/default",
            "clm/crop",
            "clm/cropMonthlyOutput",
        ]
        target = [
            "clm-default",
            "clm-default",
            "clm-crop",
            "clm-cropMonthlyOutput",
        ]
        output = _get_testmod_list(testmod_list_input, unique=False)
        self.assertEqual(output, target)

    def test_getTestmodList_suite_unique(self):
        """Ensure that _get_testmod_list() works correctly with unique=True"""
        testmod_list_input = [
            "clm/default",
            "clm/default",
            "clm/crop",
            "clm/cropMonthlyOutput",
        ]
        target = [
            "clm-default",
            "clm-crop",
            "clm-cropMonthlyOutput",
        ]

        output = _get_testmod_list(testmod_list_input, unique=True)
        self.assertEqual(output, target)

    def test_getTestmodList_testname(self):
        """Ensure that _get_testmod_list() works correctly with full test name(s) specified"""
        testmod_list_input = [
            "ERS_D_Ld15.f45_f45_mg37.I2000Clm50FatesRs.izumi_nag.clm-crop",
            "ERS_D_Ld15.f45_f45_mg37.I2000Clm50FatesRs.izumi_nag.clm-default",
        ]
        target = ["clm-crop", "clm-default"]
        output = _get_testmod_list(testmod_list_input)
        self.assertEqual(output, target)

    def test_checkArgValidity_xmlMachineWithoutSuiteName_raises(self):
        """--xml-machine without --suite-name does nothing, so it must be rejected

        This mirrors the pre-existing check that --suite-compiler requires
        --suite-name.
        """
        args = argparse.Namespace(
            suite_name=None,
            suite_compiler=None,
            xml_machine="derecho",
            rerun_existing_failures=False,
            testid_base=None,
        )
        with self.assertRaises(RuntimeError):
            _check_arg_validity(args)

    def test_getTestmodList_twomods(self):
        """
        Ensure that _get_testmod_list() works correctly with full test name(s) specified and two
        mods in one test
        """
        testmod_list_input = [
            "ERS_D_Ld15.f45_f45_mg37.I2000Clm50FatesRs.izumi_nag.clm-default--clm-crop"
        ]
        target = ["clm-default", "clm-crop"]
        output = _get_testmod_list(testmod_list_input)
        self.assertEqual(output, target)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
