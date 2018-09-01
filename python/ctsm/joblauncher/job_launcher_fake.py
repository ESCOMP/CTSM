"""A fake JobLauncher that just records the commands it is told to run"""

from collections import namedtuple
from ctsm.joblauncher.job_launcher_base import JobLauncherBase

# cmd (str): space-delimited string giving this command
# out (str): path to stdout
# err (str): path to stderr
Command = namedtuple('Command', ['cmd', 'out', 'err'])

class JobLauncherFake(JobLauncherBase):
    """A fake JobLauncher that just records the commands it is told to run"""

    def __init__(self):
        JobLauncherBase.__init__(self)
        self._commands = []

    def run_command_impl(self, command, stdout_path, stderr_path):
        self._commands.append(Command(cmd=' '.join(command),
                                      out=stdout_path,
                                      err=stderr_path))

    def run_command_logger_message(self, command, stdout_path, stderr_path):
        message = 'Appending: <{}> ' \
                  'with stdout = {} ' \
                  'and stderr = {}'.format(' '.join(command),
                                           stdout_path,
                                           stderr_path)
        return message

    def get_commands(self):
        """Return a list of commands that this job launcher has been told to run

        Each element of the list is a Command namedtuple (see above)
        """
        return self._commands
