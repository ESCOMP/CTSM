"""Machine-specific default values.

To allow running out-of-the-box on other machines, add code here."""

from collections import namedtuple
import os
from ctsm.joblauncher.job_launcher_factory import JOB_LAUNCHER_QSUB
from ctsm.machine import CREATE_TEST_QUEUE_UNSPECIFIED
from ctsm.machine_utils import get_user

MachineDefaults = namedtuple(
    "MachineDefaults",
    [
        "job_launcher_type",
        "scratch_dir",
        "baseline_dir",
        "account_required",
        "create_test_retry",
        "create_test_queue",
        "job_launcher_defaults",
    ],
)
# job_launcher_type: one of the JOB_LAUNCHERs defined in job_launcher_factory
# scratch_dir: str
# baseline_dir: str: The standard location for CTSM baselines on this machine
# job_launcher_defaults: dict: keys are the JOB_LAUNCHERs defined in job_launcher_factory,
#     values are types defined here (like _QsubDefaults). A given machine's defaults can
#     have 0, 1 or multiple job_launcher_defaults. (It can be useful to have defaults even
#     for the non-default job launcher for this machine, in case the user chooses a
#     non-default launcher.)
# create_test_retry: int: Default number of times to retry a create_test job on this machine
# create_test_queue: str: Default queue to use for create_test; if this is
#     CREATE_TEST_QUEUE_UNSPECIFIED, then we won't add a '--queue' option to create_test,
#     instead leaving that value unspecified, allowing CIME to pick an appropriate queue
#     for each test using its standard mechanisms.
# account_required: bool: whether an account number is required on this machine (not
#     really a default, but used for error-checking)

# Note that the different job launcher types have different structures defining their
# defaults, because different ones require different elements to be set. For now we only
# have defaults for qsub, because other launchers (like no_batch) don't need any
# arguments.

QsubDefaults = namedtuple("QsubDefaults", ["queue", "walltime", "extra_args", "required_args"])

MACHINE_DEFAULTS = {
    "cheyenne": MachineDefaults(
        job_launcher_type=JOB_LAUNCHER_QSUB,
        scratch_dir=os.path.join(os.path.sep, "glade", "scratch", get_user()),
        baseline_dir=os.path.join(
            os.path.sep, "glade", "p", "cgd", "tss", "To_Be_Safely_Deleted", "ctsm_baselines"
        ),
        account_required=True,
        create_test_retry=0,
        # NOTE(wjs, 2022-02-23) By default, use the regular queue, even for
        # single-processor jobs. This is because the share queue has been really flaky,
        # with lots of job failures or slow-running jobs.
        create_test_queue="regular",
        job_launcher_defaults={
            JOB_LAUNCHER_QSUB: QsubDefaults(
                queue="regular",
                walltime="11:50:00",
                extra_args="",
                # The following assumes a single node, with a single mpi proc; we may want
                # to add more flexibility in the future, making the node / proc counts
                # individually selectable
                required_args="-l select=1:ncpus=36:mpiprocs=1 -V -r n -l inception=login -k oed",
            )
        },
    ),
    "derecho": MachineDefaults(
        job_launcher_type=JOB_LAUNCHER_QSUB,
        scratch_dir=os.path.join(os.path.sep, "glade", "derecho", "scratch", get_user()),
        baseline_dir=os.path.join(os.path.sep, "glade", "campaign", "cgd", "tss", "ctsm_baselines"),
        account_required=True,
        create_test_retry=0,
        create_test_queue=CREATE_TEST_QUEUE_UNSPECIFIED,
        job_launcher_defaults={
            JOB_LAUNCHER_QSUB: QsubDefaults(
                queue="main",
                walltime="03:50:00",
                extra_args="",
                # The following assumes a single node, with a single mpi proc; we may want
                # to add more flexibility in the future, making the node / proc counts
                # individually selectable
                required_args="-l select=1:ncpus=128:mpiprocs=1 -V -r n -k oed",
            )
        },
    ),
    "hobart": MachineDefaults(
        job_launcher_type=JOB_LAUNCHER_QSUB,
        scratch_dir=os.path.join(os.path.sep, "scratch", "cluster", get_user()),
        baseline_dir=os.path.join(os.path.sep, "fs", "cgd", "csm", "ccsm_baselines"),
        account_required=False,
        create_test_retry=0,
        create_test_queue=CREATE_TEST_QUEUE_UNSPECIFIED,
        job_launcher_defaults={
            JOB_LAUNCHER_QSUB: QsubDefaults(
                queue="medium",
                walltime="04:00:00",
                extra_args="",
                required_args="-l nodes=1:ppn=48 -r n",
            )
        },
    ),
    "izumi": MachineDefaults(
        job_launcher_type=JOB_LAUNCHER_QSUB,
        scratch_dir=os.path.join(os.path.sep, "scratch", "cluster", get_user()),
        baseline_dir=os.path.join(os.path.sep, "fs", "cgd", "csm", "ccsm_baselines"),
        account_required=False,
        # jobs on izumi experience a high frequency of failures, often at the very end of
        # the job; so we'll automatically retry a failed job twice before giving up on it
        create_test_retry=2,
        create_test_queue=CREATE_TEST_QUEUE_UNSPECIFIED,
        job_launcher_defaults={
            JOB_LAUNCHER_QSUB: QsubDefaults(
                queue="medium",
                walltime="04:00:00",
                extra_args="",
                required_args="-l nodes=1:ppn=48 -r n",
            )
        },
    ),
}
