"""General-purpose utility functions"""

import logging
import sys
import string
import pdb

logger = logging.getLogger(__name__)


def abort(errmsg):
    """Abort the program with the given error message

    No traceback is given, but if the logging level is DEBUG, then we'll enter pdb
    """
    if logger.isEnabledFor(logging.DEBUG):
        pdb.set_trace()

    sys.exit("ERROR: {}".format(errmsg))


def fill_template_file(path_to_template, path_to_final, substitutions):
    """Given a template file (based on python's template strings), write a copy of the
    file with template values filled in.

    Args:
    path_to_template (str): path to the existing template file
    path_to_final (str): path to where the final version will be written
    substitutions (dict): key-value pairs for the template string substitutions
    """

    with open(path_to_template) as template_file:
        template_file_contents = template_file.read()
    template = string.Template(template_file_contents)
    final_file_contents = template.substitute(substitutions)
    with open(path_to_final, "w") as final_file:
        final_file.write(final_file_contents)


def str2bool(v):
    """
    Function for converting different forms of
    command line boolean strings to boolean value.

    Args:
        v (str): String bool input

    Raises:
        if the argument is not an acceptable boolean string
        (such as yes or no ; true or false ; y or n ; t or f ; 0 or 1).
        argparse.ArgumentTypeError: The string should be one of the mentioned values.

    Returns:
        bool: Boolean value corresponding to the input.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError(
            "Boolean value expected. [true or false] or [y or n]"
        )
