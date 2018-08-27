"""Base JobLauncher class

This should not be instantiated directly - rather, it should be extended by concrete
classes
"""

class JobLauncherBase(object):

    def __init__(self, queue=None, walltime=None, account=None,
                 required_args=None, extra_args=None):

        self._queue = queue
        self._walltime = walltime
        self._account = account
        self._required_args = required_args
        self._extra_args = extra_args

    def get_queue(self):
        return self._queue

    def get_walltime(self):
        return self._walltime

    def get_account(self):
        return self._account

    def get_required_args(self):
        return self._required_args

    def get_extra_args(self):
        return self._extra_args

    def run_command(self, command, dry_run):
        raise NotImplementedError

    def __repr__(self):
        return ("(queue='{queue}', "
                "walltime='{walltime}', "
                "account='{account}', "
                "required_args='{required_args}', "
                "extra_args='{extra_args}')".format(queue=self._queue,
                                                  walltime=self._walltime,
                                                  account=self._account,
                                                  required_args=self._required_args,
                                                  extra_args=self._extra_args))
