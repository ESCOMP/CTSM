"""A JobLauncher for systems where we can run a big job directly on the login node"""

import subprocess
import signal
import os
from ctsm.joblauncher.job_launcher_base import JobLauncherBase


class JobLauncherNoBatch(JobLauncherBase):
    """Job launcher for systems where we can run a big job directly on the login node

    This runs the job in the background, allowing the python process to continue

    Some parts of this implementation are POSIX-only, so this won't work on Windows
    """

    def __init__(self, nice_level=None):
        """
        Args:
        nice_level: Level used for the nice command (larger = lower priority, up to 19)
        """
        JobLauncherBase.__init__(self)
        self._processes = []
        if nice_level is None:
            # We have a default value of None rather than 0 in the argument list because
            # objects of this class may be created via a few layers, and we want to allow
            # for the possibility that the passed-down argument has a default of None in
            # one of these creation layers.
            self._nice_level = 0
        else:
            self._nice_level = nice_level

    def get_nice_level(self):
        """Returns the level used for the nice command (larger = lower priority)"""
        return self._nice_level

    def _preexec(self):
        """This function is run in the child process just before the actual command is run"""
        # This is equivalent to running a command with 'nohup'
        signal.signal(signal.SIGHUP, signal.SIG_IGN)

        # This is equivalent to running a command with 'nice'; note that this is POSIX-only
        os.nice(self._nice_level)

    def run_command_impl(self, command, stdout_path, stderr_path):
        with open(stdout_path, "w") as outfile, open(stderr_path, "w") as errfile:
            # Note that preexec_fn is POSIX-only; also, it may be unsafe in the presence
            # of threads (hence the need for disabling the pylint warning)
            # pylint: disable=subprocess-popen-preexec-fn
            process = subprocess.Popen(
                command, stdout=outfile, stderr=errfile, preexec_fn=self._preexec
            )
            self._processes.append(process)

    def run_command_logger_message(self, command, stdout_path, stderr_path):
        message = (
            "Running: <{command_str}> "
            "with stdout = {outfile} "
            "and stderr = {errfile}".format(
                command_str=" ".join(command), outfile=stdout_path, errfile=stderr_path
            )
        )
        return message

    def supports_waiting(self):
        """This launcher can wait for its launched process(es) to complete"""
        return True

    def wait_for_processes_to_complete(self):
        """Waits for every process started by run_command_impl to complete

        We wait for ALL of them (not just the first one found to have failed),
        so that - e.g. when running inside a container - every launched
        create_test gets a chance to finish before we return and the caller
        potentially exits, which inside a container would tear down the PID
        namespace and kill any processes still running.

        Returns the first nonzero return code encountered, in launch order, or
        0 if all of them succeeded (or if none was started - which is the
        case after a dry run).
        """
        first_failure = 0
        for process in self._processes:
            return_code = process.wait()
            if return_code != 0 and first_failure == 0:
                first_failure = return_code
        return first_failure
