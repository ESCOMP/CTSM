"""Factory function and class for creating and storing machine-related information.

Note that this module contains machine-specific information. To allow running
out-of-the-box on other machines, add code here.
"""
# TODO(wjs, 2018-08-26) Consider moving these machine-specific settings into a config file

import logging
import os
from collections import namedtuple
from ctsm.joblauncher.job_launcher_factory import create_job_launcher, JOB_LAUNCHER_NOBATCH, JOB_LAUNCHER_QSUB
from ctsm.machine_utils import get_user
from CIME.utils import get_charge_account

logger = logging.getLogger(__name__)

Machine = namedtuple('Machine', ['name',           # str
                                 'scratch_dir',    # str
                                 'account',        # str or None
                                 'job_launcher'])  # subclass of JobLauncherBase

def create_machine(machine_name, scratch_dir=None, account=None,
                   job_launcher_queue=None, job_launcher_walltime=None,
                   job_launcher_extra_args=None,
                   allow_missing_entries=False):
    """

    Args:
    allow_missing_entries (bool): For a machine that generally requires certain entries
        (e.g., account or scratch_dir): If allow_missing_entries is True, then we proceed
        even if these entries are missing. This is intended for when create_machine is
        just called for the sake of getting default values.
    """

    job_launcher_required_args = None
    if machine_name == 'hobart':
        job_launcher_type = JOB_LAUNCHER_NOBATCH
        if scratch_dir is None:
            scratch_dir = os.path.join(os.path.sep, 'scratch', 'cluster', get_user())
    elif machine_name == 'cheyenne':
        job_launcher_type = JOB_LAUNCHER_QSUB
        if scratch_dir is None:
            scratch_dir = os.path.join(os.path.sep, 'glade', 'scratch', get_user())
        if job_launcher_queue is None:
            job_launcher_queue = 'regular'
        if job_launcher_walltime is None:
            job_launcher_walltime = '06:00:00'
        if job_launcher_extra_args is None:
            job_launcher_extra_args = ''
        job_launcher_required_args = '-l select=1:ncpus=36:mpiprocs=1 -r n -l inception=login'
    else:
        if not allow_missing_entries:
            # This isn't exactly a missing entry, but the times we don't care about this
            # warning tend to be the same as the times when allow_missing_entries is true
            logger.warning("machine {} not recognized; using generic no-batch settings".format(
                machine_name))
        job_launcher_type = JOB_LAUNCHER_NOBATCH
        if scratch_dir is None and not allow_missing_entries:
            raise RuntimeError("For generic machine, must specify scratch directory")

    if account is None:
        account = _get_account()

    job_launcher = create_job_launcher(job_launcher_type=job_launcher_type,
                                       account=account,
                                       queue=job_launcher_queue,
                                       walltime=job_launcher_walltime,
                                       required_args=job_launcher_required_args,
                                       extra_args=job_launcher_extra_args,
                                       allow_missing_entries=allow_missing_entries)

    return Machine(name=machine_name,
                   scratch_dir=scratch_dir,
                   account=account,
                   job_launcher=job_launcher)

def _get_account():
    account = get_charge_account()
    return account
