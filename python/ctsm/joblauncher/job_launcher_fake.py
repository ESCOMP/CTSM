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
        self._return_codes = None

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
        """Set the value returned by wait_for_processes_to_complete

        Intended for tests that only launch a single process; for tests that launch
        multiple processes (e.g. a multi-compiler suite) and need per-launch return
        codes, use set_return_codes instead.
        """
        self._return_code = return_code

    def set_return_codes(self, return_codes):
        """Set a list of per-launch return codes (in launch order), used by
        wait_for_processes_to_complete to simulate multiple launched processes
        finishing with different statuses. Overrides set_return_code.
        """
        self._return_codes = return_codes

    def supports_waiting(self):
        """This launcher can wait for its launched process(es) to complete"""
        return True

    def wait_for_processes_to_complete(self):
        """Pretend to wait for all launched processes.

        If set_return_codes was used, returns the first nonzero code among those (in
        launch order), or 0 if none is nonzero. Otherwise returns the single canned
        return code set via set_return_code (0 by default).
        """
        if self._return_codes is not None:
            for return_code in self._return_codes:
                if return_code != 0:
                    return return_code
            return 0
        return self._return_code
