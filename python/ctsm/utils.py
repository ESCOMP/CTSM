"""General-purpose utility functions"""

import logging
import sys

logger = logging.getLogger(__name__)

def abort(errmsg):
    """Abort the program with the given error message

    No traceback is given, but if the logging level is DEBUG, then we'll enter pdb
    """
    if logger.isEnabledFor(logging.DEBUG):
        import pdb
        pdb.set_trace()

    sys.exit('ERROR: {}'.format(errmsg))
