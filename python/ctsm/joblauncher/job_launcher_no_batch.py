"""A JobLauncher for systems where we can run a big job directly on the login node"""

from ctsm.joblauncher.job_launcher_base import JobLauncherBase

class JobLauncherNoBatch(JobLauncherBase):
    """Job launcher for systems where we can run a big job directly on the login node"""

    def __init__(self):
        JobLauncherBase.__init__(self)

    def run_command(self, command, dry_run=False):
        # FIXME(wjs, 2018-08-25) Implement this
        pass
