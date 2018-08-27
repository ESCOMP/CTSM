"""Factory function for creating an appropriate job launcher class for a given machine.
"""

from __future__ import print_function

import logging
from ctsm.joblauncher.job_launcher_no_batch import JobLauncherNoBatch
from ctsm.joblauncher.job_launcher_qsub import JobLauncherQsub

logger = logging.getLogger(__name__)

# Possible values of job_launcher_type
JOB_LAUNCHER_NOBATCH = "no_batch"
JOB_LAUNCHER_QSUB = "qsub"

def create_job_launcher(job_launcher_type, account=None, queue=None, walltime=None,
                        required_args=None, extra_args=None,
                        allow_missing_entries=False):
    """

    Args:
    account (str): Account number for launching jobs. This must be provided (not None) if
        we're using a job launcher other than JOB_LAUNCHER_NOBATCH (unless
        allow_missing_entries is True)
    allow_missing_entries (bool): For a job launcher type that generally requires certain
        entries (e.g., account): If allow_missing_entries is True, then we proceed even if
        these entries are missing. This is intended for when create_job_launcher is just
        called for the sake of getting default values.
    """

    if job_launcher_type is JOB_LAUNCHER_NOBATCH:
        return JobLauncherNoBatch()

    if account is None and not allow_missing_entries:
        raise RuntimeError("Could not find an account code")

    if job_launcher_type is JOB_LAUNCHER_QSUB:
        if not allow_missing_entries:
            _assert_not_none(queue, 'queue', job_launcher_type)
            _assert_not_none(walltime, 'walltime', job_launcher_type)
            _assert_not_none(account, 'account', job_launcher_type)
            # required_args and extra_args can be the empty string, but not None
            _assert_not_none(required_args, 'required_args', job_launcher_type)
            _assert_not_none(extra_args, 'extra_args', job_launcher_type)
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
