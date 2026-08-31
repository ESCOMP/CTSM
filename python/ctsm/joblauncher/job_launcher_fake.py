"""A fake JobLauncher that just records the commands it is told to run"""

from collections import namedtuple
from ctsm.joblauncher.job_launcher_base import JobLauncherBase

# cmd (str): space-delimited string giving this command
# out (str): path to stdout
# err (str): path to stderr
Command = namedtuple("Command", ["cmd", "out", "err"])


class JobLauncherFake(JobLauncherBase):
    """A fake JobLauncher that just records the commands it is told to run"""

    def __init__(self):
        JobLauncherBase.__init__(self)
        self._commands = []
        self._return_code = 0

    def run_command_impl(self, command, stdout_path, stderr_path):
        self._commands.append(Command(cmd=" ".join(command), out=stdout_path, err=stderr_path))

    def run_command_logger_message(self, command, stdout_path, stderr_path):
        message = (
            "Appending: <{}> "
            "with stdout = {} "
            "and stderr = {}".format(" ".join(command), stdout_path, stderr_path)
        )
        return message

    def get_commands(self):
        """Return a list of commands that this job launcher has been told to run

        Each element of the list is a Command namedtuple (see above)
        """
        return self._commands

    def set_return_code(self, return_code):
        """Set the value returned by wait_for_last_process_to_complete"""
        self._return_code = return_code

    def wait_for_last_process_to_complete(self):
        """Pretend to wait; returns the canned return code set by set_return_code"""
        return self._return_code
