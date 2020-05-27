"""Factory function for creating an appropriate job launcher class for a given machine.
"""

from __future__ import print_function

import logging
from ctsm.joblauncher.job_launcher_no_batch import JobLauncherNoBatch
from ctsm.joblauncher.job_launcher_qsub import JobLauncherQsub
from ctsm.joblauncher.job_launcher_fake import JobLauncherFake

logger = logging.getLogger(__name__)

# Possible values of job_launcher_type
JOB_LAUNCHER_NOBATCH = "no_batch"
JOB_LAUNCHER_QSUB = "qsub"
JOB_LAUNCHER_FAKE = "fake"

def create_job_launcher(job_launcher_type,
                        account=None, queue=None, walltime=None,
                        nice_level=None,
                        required_args=None, extra_args=None,
                        allow_missing_entries=False):
    """
    Creates and returns a job launcher object of the specified type

    Args:
    job_launcher_type: One of the JOB_LAUNCHER constants defined at the top of this module
    account (str or None): Account number for launching jobs. Not applicable for
        JOB_LAUNCHER_NOBATCH and JOB_LAUNCHER_FAKE. For other job launchers, this must be
        provided if an account number is needed to launch jobs on this machine.
    queue (str or None): Queue to use. Not applicable for JOB_LAUNCHER_NOBATCH and
        JOB_LAUNCHER_FAKE; must be provided for other job launcher types.
    walltime (str or None): Walltime to use. Not applicable for JOB_LAUNCHER_NOBATCH and
        JOB_LAUNCHER_FAKE; much be provided for other job launcher types.
    nice_level (int or None): Level used for the nice command; only applicable for
        JOB_LAUNCHER_NOBATCH
    required_args (str or None): Arguments to the job launcher that cannot be overridden by the
        user. Not applicable for JOB_LAUNCHER_NOBATCH and JOB_LAUNCHER_FAKE.
    extra_args (str or None): Arguments to the job launcher that can be overridden by the
        user. Not applicable for JOB_LAUNCHER_NOBATCH and JOB_LAUNCHER_FAKE.
    allow_missing_entries (bool): For a job launcher type that generally requires certain
        entries (e.g., queue): If allow_missing_entries is True, then we proceed even if
        these entries are missing. This is intended for when create_job_launcher is just
        called for the sake of getting default values.
    """

    if job_launcher_type == JOB_LAUNCHER_FAKE:
        return JobLauncherFake()

    if job_launcher_type == JOB_LAUNCHER_NOBATCH:
        return JobLauncherNoBatch(nice_level=nice_level)

    if job_launcher_type == JOB_LAUNCHER_QSUB:
        if not allow_missing_entries:
            _assert_not_none(queue, 'queue', job_launcher_type)
            _assert_not_none(walltime, 'walltime', job_launcher_type)
        return JobLauncherQsub(queue=queue,
                               walltime=walltime,
                               account=account,
                               required_args=required_args,
                               extra_args=extra_args)

    raise RuntimeError("Unexpected job launcher type: {}".format(job_launcher_type))

def _assert_not_none(arg, arg_name, job_launcher_type):
    """Raises an exception if the given argument has value None"""
    if arg is None:
        raise TypeError("{} cannot be None for job launcher of type {}".format(
            arg_name, job_launcher_type))
