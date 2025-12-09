"""Various OS-related utility functions"""

import os
import subprocess
from ctsm.utils import abort


def run_cmd_output_on_error(cmd, errmsg, cwd=None):
    """Run the given command; suppress output but print it if there is an error

    If there is an error running the command, print the output from the command and abort
    with the given errmsg.

    Args:
    cmd: list of strings - command and its arguments
    errmsg: string - error message to print if the command returns an error code
    cwd: string or None - path from which the command should be run
    """
    try:
        _ = subprocess.check_output(cmd, stderr=subprocess.STDOUT, universal_newlines=True, cwd=cwd)
    except subprocess.CalledProcessError as error:
        print("ERROR while running:")
        print(" ".join(cmd))
        if cwd is not None:
            print("From {}".format(cwd))
        print("")
        print(error.output)
        print("")
        abort(errmsg)
    except:
        print("ERROR trying to run:")
        print(" ".join(cmd))
        if cwd is not None:
            print("From {}".format(cwd))
        raise


def make_link(src, dst):
    """Makes a link pointing to src named dst

    Does nothing if link is already set up correctly
    """
    if os.path.islink(dst) and os.readlink(dst) == src:
        # Link is already set up correctly: do nothing (os.symlink raises an exception if
        # you try to replace an existing file)
        pass
    else:
        os.symlink(src, dst)
