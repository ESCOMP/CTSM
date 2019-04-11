#!/usr/bin/env python

"""Unit tests for machine
"""

import unittest
import os

from ctsm import add_cime_to_path # pylint: disable=unused-import
from ctsm import unit_testing

from ctsm.machine import create_machine
from ctsm.machine_utils import get_user
from ctsm.machine_defaults import MACHINE_DEFAULTS, MachineDefaults, QsubDefaults
from ctsm.joblauncher.job_launcher_no_batch import JobLauncherNoBatch
from ctsm.joblauncher.job_launcher_qsub import JobLauncherQsub
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_QSUB, JOB_LAUNCHER_NOBATCH

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

class TestCreateMachine(unittest.TestCase):
    """Tests of create_machine"""

    def assertMachineInfo(self, machine, name, scratch_dir, account):
        """Asserts that the basic machine info is as expected.

        This does NOT dive down into the job launcher"""
        self.assertEqual(machine.name, name)
        self.assertEqual(machine.scratch_dir, scratch_dir)
        self.assertEqual(machine.account, account)

    def assertNoBatchInfo(self, machine, nice_level=None):
        """Asserts that the machine's launcher is of type JobLauncherNoBatch"""
        launcher = machine.job_launcher
        self.assertIsInstance(launcher, JobLauncherNoBatch)
        if nice_level is None:
            # the default nice level should be 0
            nice_level = 0
        self.assertEqual(launcher.get_nice_level(), nice_level)

    def assertQsubInfo(self, machine, queue, walltime, account, required_args, extra_args):
        """Asserts that the machine's launcher is of type JobLauncherQsub, and has values as
        expected"""
        launcher = machine.job_launcher
        self.assertIsInstance(launcher, JobLauncherQsub)
        self.assertEqual(launcher.get_queue(), queue)
        self.assertEqual(launcher.get_walltime(), walltime)
        self.assertEqual(launcher.get_account(), account)
        self.assertEqual(launcher.get_required_args(), required_args)
        self.assertEqual(launcher.get_extra_args(), extra_args)

    @staticmethod
    def create_defaults(default_job_launcher=JOB_LAUNCHER_QSUB):
        """Creates test-specific defaults so we don't tie the tests to changes in the real
        defaults"""
        defaults = {
            'cheyenne': MachineDefaults(
                job_launcher_type=default_job_launcher,
                scratch_dir=os.path.join(os.path.sep, 'glade', 'scratch', get_user()),
                account_required=True,
                job_launcher_defaults={
                    JOB_LAUNCHER_QSUB: QsubDefaults(
                        queue='regular',
                        walltime='06:00:00',
                        extra_args='',
                        required_args=
                        '-l select=1:ncpus=36:mpiprocs=1 -r n -l inception=login')
                })
            }
        return defaults

    def test_unknownMachine_defaults(self):
        """Tests a machine not in the defaults structure, with no overriding arguments"""
        machine = create_machine('unknown_test_machine', MACHINE_DEFAULTS,
                                 account='a123')
        self.assertMachineInfo(machine=machine,
                               name='unknown_test_machine',
                               scratch_dir=None,
                               account='a123')
        self.assertNoBatchInfo(machine)

    def test_noBatchMachine_niceLevel(self):
        """Tests a no-batch machine where the nice level is explicit"""
        machine = create_machine('unknown_test_machine', MACHINE_DEFAULTS,
                                 job_launcher_type=JOB_LAUNCHER_NOBATCH,
                                 scratch_dir='/path/to/scratch',
                                 account='a123',
                                 job_launcher_nice_level=13)
        self.assertMachineInfo(machine=machine,
                               name='unknown_test_machine',
                               scratch_dir='/path/to/scratch',
                               account='a123')
        self.assertNoBatchInfo(machine, nice_level=13)

    def test_unknownMachine_argsExplicit(self):
        """Tests a machine not in the defaults structure, with explicit arguments"""
        machine = create_machine('unknown_test_machine', MACHINE_DEFAULTS,
                                 scratch_dir='/path/to/scratch',
                                 job_launcher_type=JOB_LAUNCHER_QSUB,
                                 account='a123',
                                 job_launcher_queue='my_queue',
                                 job_launcher_walltime='1:23:45',
                                 job_launcher_extra_args='--some args')
        self.assertMachineInfo(machine=machine,
                               name='unknown_test_machine',
                               scratch_dir='/path/to/scratch',
                               account='a123')
        self.assertQsubInfo(machine=machine,
                            queue='my_queue',
                            walltime='1:23:45',
                            account='a123',
                            required_args='',
                            extra_args='--some args')

    def test_knownMachine_defaults(self):
        """Tests a machine known in the defaults structure, with no overriding arguments"""
        defaults = self.create_defaults()
        machine = create_machine('cheyenne', defaults, account='a123')
        self.assertMachineInfo(machine=machine,
                               name='cheyenne',
                               scratch_dir=os.path.join(os.path.sep,
                                                        'glade',
                                                        'scratch',
                                                        get_user()),
                               account='a123')
        self.assertQsubInfo(machine=machine,
                            queue='regular',
                            walltime='06:00:00',
                            account='a123',
                            required_args='-l select=1:ncpus=36:mpiprocs=1 -r n -l inception=login',
                            extra_args='')

    def test_knownMachine_argsExplicit(self):
        """Tests a machine known in the defaults structure, with explicit arguments"""
        defaults = self.create_defaults(default_job_launcher=JOB_LAUNCHER_NOBATCH)
        machine = create_machine('cheyenne', defaults,
                                 job_launcher_type=JOB_LAUNCHER_QSUB,
                                 scratch_dir='/custom/path/to/scratch',
                                 account='a123',
                                 job_launcher_queue='custom_queue',
                                 job_launcher_walltime='9:87:65',
                                 job_launcher_extra_args='--custom args')
        self.assertMachineInfo(machine=machine,
                               name='cheyenne',
                               scratch_dir='/custom/path/to/scratch',
                               account='a123')
        self.assertQsubInfo(machine=machine,
                            queue='custom_queue',
                            walltime='9:87:65',
                            account='a123',
                            required_args='-l select=1:ncpus=36:mpiprocs=1 -r n -l inception=login',
                            extra_args='--custom args')

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
