"""A JobLauncher for systems where we can run a big job directly on the login node"""

from ctsm.joblauncher.job_launcher_base import JobLauncherBase

class JobLauncherNoBatch(JobLauncherBase):
    """Job launcher for systems where we can run a big job directly on the login node"""

    def __init__(self, nice_level=None):
        """
        Args:
        nice_level: Level used for the nice command (larger = lower priority, up to 19)
        """
        JobLauncherBase.__init__(self)
        if nice_level is None:
            self._nice_level = 0
        else:
            self._nice_level = nice_level

    def get_nice_level(self):
        """Returns the level used for the nice command (larger = lower priority)"""
        return self._nice_level

    def run_command(self, command, dry_run=False):
        # FIXME(wjs, 2018-08-25) Implement this
        pass
