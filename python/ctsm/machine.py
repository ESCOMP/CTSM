"""Factory function and class for creating and storing machine-related information.
"""

import logging
from collections import namedtuple
from ctsm.joblauncher.job_launcher_factory import \
    create_job_launcher, JOB_LAUNCHER_NOBATCH
from CIME.utils import get_project  # pylint: disable=import-error

logger = logging.getLogger(__name__)

# TODO(wjs, 2018-08-31) Turn this into a real class, with getter methods.
#
# An attempt to call get_scratch_dir() should raise a RuntimeError if self._scratch_dir is
# None. (That is, we allow a machine to be created with a scratch_dir of None - because
# for some applications we don't need a scratch_dir - but we raise an exception if the
# application tries to use the scratch_dir when it wasn't properly set, in order to print
# a meaningful error message rather than whatever message you'd get when trying to use the
# None value somewhere.)
#
# For now, just keep in mind that it's possible that a machine's scratch_dir may be None,
# if it wasn't set explicitly and there was no default available. For now, it's up to the
# user of the machine object to check for that possibility if need be.
#
# Similar notes apply to baseline_dir.
Machine = namedtuple('Machine', ['name',           # str
                                 'scratch_dir',    # str
                                 'baseline_dir',   # str
                                 'account',        # str or None
                                 'job_launcher'])  # subclass of JobLauncherBase

def create_machine(machine_name, defaults, job_launcher_type=None,
                   scratch_dir=None, account=None,
                   job_launcher_queue=None, job_launcher_walltime=None,
                   job_launcher_nice_level=None,
                   job_launcher_extra_args=None,
                   allow_missing_entries=False):
    """Create a machine object (of type Machine, as given above)

    This uses the provided (non-None) arguments to override any defaults provided via the
    'defaults' argument.

    Args:
    machine_name (str): name of machine; this is used to index into the 'defaults'
        argument, and is used as the machine name in the returned object
    defaults: dict of MachineDefaults (as defined in machine_defaults)
    scratch_dir: path to scratch directory (if not provided, will attempt to get it from
        machine defaults)
    account: account to use for job submission to a queue (if not provided, will attempt
        to determine it automatically using the CIME method)
    job_launcher_type: one of the JOB_LAUNCHER constants defined in job_launcher_factory,
        or None. If None, we pick the default for this machine.
    job_launcher_queue (str or None): Queue to use (not applicable for
        JOB_LAUNCHER_NOBATCH and JOB_LAUNCHER_FAKE)
    job_launcher_walltime (str or None): Walltime to use (not applicable for
        JOB_LAUNCHER_NOBATCH and JOB_LAUNCHER_FAKE)
    job_launcher_nice_level (int or None): Level used for the nice command; only
        applicable for JOB_LAUNCHER_NOBATCH
    job_launcher_extra_args (str or None): Arguments to the job launcher that can be
        overridden by the user. Not applicable for JOB_LAUNCHER_NOBATCH and
        JOB_LAUNCHER_FAKE
    allow_missing_entries (bool): For a machine that generally requires certain entries
        (e.g., account): If allow_missing_entries is True, then we proceed even if these
        entries are missing. This is intended for when create_machine is just called for
        the sake of getting default values.
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
    baseline_dir = None
    if mach_defaults is not None:
        if job_launcher_type is None:
            job_launcher_type = mach_defaults.job_launcher_type
        if scratch_dir is None:
            scratch_dir = mach_defaults.scratch_dir
        # NOTE(wjs, 2019-05-17) Note that we don't provide a way to override the default
        # baseline_dir. The idea is that this should always provide the standard baseline
        # directory for this machine, even if a particular application wants to use
        # something different (in contrast to, say, scratch_dir, for which the value in
        # the machines object can itself be overridden). I think this will be smoother and
        # more intuitive if different baseline directories are needed for different
        # purposes in a single application (e.g., different baseline directories for
        # generation and comparison, or making a link in some temporary location that
        # points to the standard baselines).
        baseline_dir = mach_defaults.baseline_dir
        if account is None and mach_defaults.account_required and not allow_missing_entries:
            raise RuntimeError("Could not find an account code")
    else:
        if not allow_missing_entries:
            # This isn't exactly a missing entry, but the times we don't care about this
            # warning tend to be the same as the times when allow_missing_entries is true
            logger.warning("machine %s not recognized; using generic no-batch settings",
                           machine_name)
        if job_launcher_type is None:
            job_launcher_type = JOB_LAUNCHER_NOBATCH

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
                                       nice_level=job_launcher_nice_level,
                                       required_args=job_launcher_required_args,
                                       extra_args=job_launcher_extra_args,
                                       allow_missing_entries=allow_missing_entries)

    return Machine(name=machine_name,
                   scratch_dir=scratch_dir,
                   baseline_dir=baseline_dir,
                   account=account,
                   job_launcher=job_launcher)

def get_possibly_overridden_baseline_dir(machine, baseline_dir=None):
    """Get the baseline directory to use here, or None

    If baseline_dir is provided (not None), use that. Otherwise use the baseline directory
    from machine (which may be None).

    Args:
    machine (Machine)
    baseline_dir (str or None): gives the overriding baseline directory to use
    """
    if baseline_dir is not None:
        return baseline_dir
    return machine.baseline_dir

def _get_account():
    account = get_project()
    return account
