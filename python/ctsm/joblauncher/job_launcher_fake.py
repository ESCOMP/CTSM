"""A fake JobLauncher that just records the commands it is told to run"""

from ctsm.joblauncher.job_launcher_base import JobLauncherBase

class JobLauncherFake(JobLauncherBase):
    """A fake JobLauncher that just records the commands it is told to run"""

    def __init__(self):
        JobLauncherBase.__init__(self)
        self._commands = []

    def run_command(self, command, dry_run=False):
        self._commands.append(' '.join(command))

    def get_commands(self):
        """Return a list of commands that this job launcher has been told to run

        Each element of the list is a space-delimited string giving this command
        """
        return self._commands
