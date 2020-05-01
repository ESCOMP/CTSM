#!/usr/bin/env python

"""Unit tests for job_launcher_no_batch
"""

import unittest
import tempfile
import shutil
import os
from ctsm.joblauncher.job_launcher_factory import create_job_launcher, JOB_LAUNCHER_NOBATCH

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

class TestJobLauncherNoBatch(unittest.TestCase):
    """Tests of job_launcher_no_batch"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def assertFileContentsEqual(self, expected, filepath, msg=None):
        """Asserts that the contents of the file given by 'filepath' are equal to
        the string given by 'expected'. 'msg' gives an optional message to be
        printed if the assertion fails."""

        with open(filepath, 'r') as myfile:
            contents = myfile.read()

        self.assertEqual(expected, contents, msg=msg)

    def test_runCommand(self):
        """Test that the given command gets executed"""
        job_launcher = create_job_launcher(job_launcher_type=JOB_LAUNCHER_NOBATCH)
        stdout = os.path.join(self._testdir, 'stdout')
        job_launcher.run_command(command=['echo', 'hello', 'world'],
                                 stdout_path=stdout,
                                 stderr_path=os.path.join(self._testdir, 'stderr'))
        job_launcher.wait_for_last_process_to_complete()
        self.assertTrue(os.path.isfile(stdout))
        self.assertFileContentsEqual('hello world\n', stdout)

    def test_runCommand_dryRun(self):
        """With dry_run, testdir should be empty"""
        job_launcher = create_job_launcher(job_launcher_type=JOB_LAUNCHER_NOBATCH)
        job_launcher.run_command(command=['echo', 'hello', 'world'],
                                 stdout_path=os.path.join(self._testdir, 'stdout'),
                                 stderr_path=os.path.join(self._testdir, 'stderr'),
                                 dry_run=True)
        # There shouldn't be a "last process", but in case there is, wait for it to
        # complete so we can be confident that the test isn't passing simply because the
        # process hasn't completed yet. (This relies on there being logic in
        # wait_for_last_process_to_complete so that it succeeds even if there is no "last
        # process".)
        job_launcher.wait_for_last_process_to_complete()
        self.assertEqual(os.listdir(self._testdir), [])
