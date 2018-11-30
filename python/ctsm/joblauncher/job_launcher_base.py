"""Base JobLauncher class

This should not be instantiated directly - rather, it should be extended by concrete
classes.

Classes that extend this should provide implementations of run_command_impl and
run_command_logger_message.
"""

import logging

logger = logging.getLogger(__name__)

class JobLauncherBase(object):
    """Base class for job launchers. Not meant to be instantiated directly"""

    def __init__(self, queue=None, walltime=None, account=None,
                 required_args=None, extra_args=None):
        """Initialize a job launcher object.

        Note that some of these arguments (e.g., queue and walltime) aren't needed by some
        job launcher types. They are included in the base class so that users of the class
        can call the getter methods without needing to worry about knowing what type they
        have (e.g., this is used to print default values in help messages; if we're using
        a class that doesn't have a particular value, it might simply print a default of
        'None').

        For a job launcher type that expects a particular argument (e.g., extra_args): if
        that argument is unneeded in a given instance, use an empty string rather than
        None. None should be reserved for job launcher types that do not reference that
        argument at all.

        Args:
        queue: str or None
        walltime: str or None
        account: str or None
        required_args: str or None: arguments to the job launcher that cannot be
            overridden by the user
        extra_args: str or None: arguments to the job launcher that can be set or
            overridden by the user

        (required_args and extra_args are separated so that extra_args can be completely
        replaced without affecting required_args)
        """
        self._queue = queue
        self._walltime = walltime
        self._account = account
        self._required_args = required_args
        self._extra_args = extra_args

    def get_queue(self):
        """Return the queue (str or None)"""
        return self._queue

    def get_walltime(self):
        """Return the walltime (str or None)"""
        return self._walltime

    def get_account(self):
        """Return the account (str or None)"""
        return self._account

    def get_required_args(self):
        """Return the launcher's required arguments (str or None)

        These are arguments to the job launcher that cannot be overridden by the user
        """
        return self._required_args

    def get_extra_args(self):
        """Return the launcher's extra arguments (str or None)

        These are arguments to the job launcher that can be set or overridden by the user
        """
        return self._extra_args

    def run_command(self, command, stdout_path, stderr_path, dry_run=False):
        """Run a command with this job launcher.

        Command should be a list (as is typically passed to subprocess calls)

        stdout_path and stderr_path are the paths to the files that will hold stdout and
        stderr from the job

        If dry_run is True, then just print the command to be run without actually running it.
        """
        logger.info("%s", self.run_command_logger_message(command, stdout_path, stderr_path))
        if not dry_run:
            self.run_command_impl(command, stdout_path, stderr_path)

    def run_command_impl(self, command, stdout_path, stderr_path):
        """Actually runs the command

        stdout_path and stderr_path are the paths to the files that will hold stdout and
        stderr from the job

        Command should be a list (as is typically passed to subprocess calls)
        """
        raise NotImplementedError

    def run_command_logger_message(self, command, stdout_path, stderr_path):
        """Return a string for output to the log describing the command that will be run

        stdout_path and stderr_path are the paths to the files that will hold stdout and
        stderr from the job

        Command should be a list (as is typically passed to subprocess calls)
        """
        raise NotImplementedError

    def __repr__(self):
        return (type(self).__name__ +
                "(queue='{queue}', "
                "walltime='{walltime}', "
                "account='{account}', "
                "required_args='{required_args}', "
                "extra_args='{extra_args}')".format(queue=self._queue,
                                                    walltime=self._walltime,
                                                    account=self._account,
                                                    required_args=self._required_args,
                                                    extra_args=self._extra_args))
