"""A JobLauncher for systems that use qsub"""

import logging
import subprocess

from ctsm.joblauncher.job_launcher_base import JobLauncherBase

logger = logging.getLogger(__name__)

class JobLauncherQsub(JobLauncherBase):
    """Job launcher for systems where we run big jobs via qsub"""

    def __init__(self, queue, walltime, account, required_args, extra_args):
        """
        account can be None or the empty string on a machine that doesn't require an account
        """
        JobLauncherBase.__init__(self,
                                 queue=queue,
                                 walltime=walltime,
                                 account=account,
                                 required_args=required_args,
                                 extra_args=extra_args)

    def run_command_impl(self, command, stdout_path, stderr_path):
        qsub_process = subprocess.Popen(self._qsub_command(stdout_path, stderr_path),
                                        stdin=subprocess.PIPE,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
        (out, err) = qsub_process.communicate(' '.join(command))
        if err:
            logger.info('qsub ERROR:\n%s', err)
        if out:
            logger.info('qsub OUTPUT:\n%s', out)

    def run_command_logger_message(self, command, stdout_path, stderr_path):
        message = 'Running: <{command_str}> ' \
                  'via <{qsub_str}>'.format(
                      command_str=' '.join(command),
                      qsub_str=' '.join(self._qsub_command(stdout_path,
                                                           stderr_path)))
        return message

    def _qsub_command(self, stdout_path, stderr_path):
        """Returns the qsub command"""
        qsub_cmd = ['qsub',
                    '-q', self._queue,
                    '-l', 'walltime={}'.format(self._walltime),
                    '-o', stdout_path,
                    '-e', stderr_path]
        if self._account:
            qsub_cmd.extend(['-A', self._account])
        if self._required_args:
            qsub_cmd.extend(self._required_args.split())
        if self._extra_args:
            qsub_cmd.extend(self._extra_args.split())
        return qsub_cmd
