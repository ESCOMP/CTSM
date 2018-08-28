"""Factory function and class for creating and storing machine-related information.
"""

import logging
from collections import namedtuple
from ctsm.joblauncher.job_launcher_factory import \
    create_job_launcher, JOB_LAUNCHER_NOBATCH
from CIME.utils import get_charge_account  # pylint:disable=import-error

logger = logging.getLogger(__name__)

Machine = namedtuple('Machine', ['name',           # str
                                 'scratch_dir',    # str
                                 'account',        # str or None
                                 'job_launcher'])  # subclass of JobLauncherBase

def create_machine(machine_name, defaults, job_launcher_type=None,
                   scratch_dir=None, account=None,
                   job_launcher_queue=None, job_launcher_walltime=None,
                   job_launcher_extra_args=None,
                   allow_missing_entries=False):
    # FIXME(wjs, 2018-08-27) finish documenting the interface
    """

    Args:
    defaults: dict of MachineDefaults (as defined in machine_defaults)
    job_launcher_type: one of the JOB_LAUNCHER constants defined in job_launcher_factory,
        or None. If None, we pick the default for this machine.
    allow_missing_entries (bool): For a machine that generally requires certain entries
        (e.g., account or scratch_dir): If allow_missing_entries is True, then we proceed
        even if these entries are missing. This is intended for when create_machine is
        just called for the sake of getting default values.
    """

    # ------------------------------------------------------------------------
    # Settings that are independent of both machine and job launcher type
    # ------------------------------------------------------------------------

    if account is None:
        account = _get_account()

    # ------------------------------------------------------------------------
    # Settings that depend on machine
    # ------------------------------------------------------------------------

    mach_defaults = defaults.get(machine_name)
    if mach_defaults is not None:
        if job_launcher_type is None:
            job_launcher_type = mach_defaults.job_launcher_type
        if scratch_dir is None:
            scratch_dir = mach_defaults.scratch_dir
    else:
        if not allow_missing_entries:
            # This isn't exactly a missing entry, but the times we don't care about this
            # warning tend to be the same as the times when allow_missing_entries is true
            logger.warning("machine %s not recognized; using generic no-batch settings",
                           machine_name)
        if job_launcher_type is None:
            job_launcher_type = JOB_LAUNCHER_NOBATCH
        if scratch_dir is None and not allow_missing_entries:
            raise RuntimeError("For generic machine, must specify scratch directory")

    # ------------------------------------------------------------------------
    # Settings that depend on both machine and job launcher type
    # ------------------------------------------------------------------------

    # Required args are ones that cannot be overridden by the user. So it doesn't make
    # sense to have these in the argument list for this function: they should only come
    # from the defaults structure. If the user wants to provide their own arguments, those
    # should be provided via extra_args.
    job_launcher_required_args = ''

    if mach_defaults is not None:
        these_defaults = mach_defaults.job_launcher_defaults.get(job_launcher_type)
        if these_defaults is not None:
            if job_launcher_queue is None:
                job_launcher_queue = these_defaults.queue
            if job_launcher_walltime is None:
                job_launcher_walltime = these_defaults.walltime
            if job_launcher_extra_args is None:
                job_launcher_extra_args = these_defaults.extra_args
            job_launcher_required_args = these_defaults.required_args

    # ------------------------------------------------------------------------
    # Create the job launcher and the full machine object
    # ------------------------------------------------------------------------

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
