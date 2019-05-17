"""Machine-specific default values.

To allow running out-of-the-box on other machines, add code here."""

from collections import namedtuple
import os
from ctsm.joblauncher.job_launcher_factory import \
    JOB_LAUNCHER_QSUB
from ctsm.machine_utils import get_user

MachineDefaults = namedtuple('MachineDefaults', ['job_launcher_type',
                                                 'scratch_dir',
                                                 'account_required',
                                                 'job_launcher_defaults'])
# job_launcher_type: one of the JOB_LAUNCHERs defined in job_launcher_factory
# scratch_dir: str
# job_launcher_defaults: dict: keys are the JOB_LAUNCHERs defined in job_launcher_factory,
#     values are types defined here (like _QsubDefaults). A given machine's defaults can
#     have 0, 1 or multiple job_launcher_defaults. (It can be useful to have defaults even
#     for the non-default job launcher for this machine, in case the user chooses a
#     non-default launcher.)
# account_required: bool: whether an account number is required on this machine (not
#     really a default, but used for error-checking)

# Note that the different job launcher types have different structures defining their
# defaults, because different ones require different elements to be set. For now we only
# have defaults for qsub, because other launchers (like no_batch) don't need any
# arguments.

QsubDefaults = namedtuple('QsubDefaults', ['queue',
                                           'walltime',
                                           'extra_args',
                                           'required_args'])

MACHINE_DEFAULTS = {
    'cheyenne': MachineDefaults(
        job_launcher_type=JOB_LAUNCHER_QSUB,
        scratch_dir=os.path.join(os.path.sep, 'glade', 'scratch', get_user()),
        account_required=True,
        job_launcher_defaults={
            JOB_LAUNCHER_QSUB: QsubDefaults(
                queue='regular',
                walltime='06:00:00',
                extra_args='',
                # The following assumes a single node, with a single mpi proc; we may want
                # to add more flexibility in the future, making the node / proc counts
                # individually selectable
                required_args=
                '-l select=1:ncpus=36:mpiprocs=1 -r n -l inception=login')
            }),
    'hobart': MachineDefaults(
        job_launcher_type=JOB_LAUNCHER_QSUB,
        scratch_dir=os.path.join(os.path.sep, 'scratch', 'cluster', get_user()),
        account_required=False,
        job_launcher_defaults={
            JOB_LAUNCHER_QSUB: QsubDefaults(
                queue='medium',
                walltime='04:00:00',
                extra_args='',
                required_args='-l nodes=1:ppn=48 -r n')
        })
}
