#!/usr/bin/env python3

"""Unit tests for job_launcher_no_batch"""

import unittest
import tempfile
import shutil
import os
from ctsm.joblauncher.job_launcher_factory import (
    create_job_launcher,
    JOB_LAUNCHER_NOBATCH,
)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestJobLauncherNoBatch(unittest.TestCase):
    """Tests of job_launcher_no_batch"""

    def setUp(self):
        self._previous_dir = os.getcwd()
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        os.chdir(self._previous_dir)
        shutil.rmtree(self._testdir, ignore_errors=True)

    def assertFileContentsEqual(self, expected, filepath, msg=None):
        """Asserts that the contents of the file given by 'filepath' are equal to
        the string given by 'expected'. 'msg' gives an optional message to be
        printed if the assertion fails."""

        with open(filepath, "r") as myfile:
            contents = myfile.read()

        self.assertEqual(expected, contents, msg=msg)

    def test_runCommand(self):
        """Test that the given command gets executed"""
        job_launcher = create_job_launcher(job_launcher_type=JOB_LAUNCHER_NOBATCH)
        stdout = os.path.join(self._testdir, "stdout")
        job_launcher.run_command(
            command=["echo", "hello", "world"],
            stdout_path=stdout,
            stderr_path=os.path.join(self._testdir, "stderr"),
        )
        job_launcher.wait_for_processes_to_complete()
        self.assertTrue(os.path.isfile(stdout))
        self.assertFileContentsEqual("hello world\n", stdout)

    def test_waitForProcesses_waitsForAllAndReportsFailure(self):
        """With multiple launched processes, wait for all of them and report a failure

        The first-launched process fails immediately while the second is still
        running. We must (a) still be waiting for that second process when the
        first has already exited, and (b) return the failing status. This is what
        protects a containerized run: run_sys_tests returning early would let
        podman tear down the PID namespace and kill the processes still going.
        """
        job_launcher = create_job_launcher(job_launcher_type=JOB_LAUNCHER_NOBATCH)
        slow_stdout = os.path.join(self._testdir, "stdout_slow")
        job_launcher.run_command(
            command=["sh", "-c", "exit 5"],
            stdout_path=os.path.join(self._testdir, "stdout_fail"),
            stderr_path=os.path.join(self._testdir, "stderr_fail"),
        )
        job_launcher.run_command(
            command=["sh", "-c", "sleep 1; echo slow done"],
            stdout_path=slow_stdout,
            stderr_path=os.path.join(self._testdir, "stderr_slow"),
        )

        return_code = job_launcher.wait_for_processes_to_complete()

        self.assertEqual(return_code, 5)
        # If we had only waited for the last process, or had returned as soon as the
        # first one failed, this file would still be empty.
        self.assertFileContentsEqual("slow done\n", slow_stdout)

    def test_runCommand_dryRun(self):
        """With dry_run, testdir should be empty"""
        job_launcher = create_job_launcher(job_launcher_type=JOB_LAUNCHER_NOBATCH)
        job_launcher.run_command(
            command=["echo", "hello", "world"],
            stdout_path=os.path.join(self._testdir, "stdout"),
            stderr_path=os.path.join(self._testdir, "stderr"),
            dry_run=True,
        )
        # There shouldn't be any launched processes, but in case there are, wait for them
        # to complete so we can be confident that the test isn't passing simply because a
        # process hasn't completed yet. (This relies on there being logic in
        # wait_for_processes_to_complete so that it succeeds even if no process was
        # launched.)
        job_launcher.wait_for_processes_to_complete()
        self.assertEqual(os.listdir(self._testdir), [])
