"""A JobLauncher for systems that use qsub"""

from ctsm.joblauncher.job_launcher_base import JobLauncherBase

class JobLauncherQsub(JobLauncherBase):
    """Job launcher for systems where we run big jobs via qsub"""

    def __init__(self, queue, walltime, account, required_args, extra_args):
        JobLauncherBase.__init__(self,
                                 queue=queue,
                                 walltime=walltime,
                                 account=account,
                                 required_args=required_args,
                                 extra_args=extra_args)

    def run_command_impl(self, command):
        # FIXME(wjs, 2018-08-25) Implement this
        pass
