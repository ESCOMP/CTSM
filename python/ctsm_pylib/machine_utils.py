"""Various machine-related utility functions
"""

from __future__ import print_function

import getpass
import socket
import re

# ========================================================================
# Public functions
# ========================================================================

def get_user():
    """Return the current user name (string)"""
    return getpass.getuser()

def get_machine_name():
    """Return the current machine name (string)"""
    full_hostname = socket.gethostname()
    hostname = full_hostname.split('.')[0]
    return _machine_from_hostname(hostname)

# ========================================================================
# Private functions
# ========================================================================

def _machine_from_hostname(hostname):
    """Given a hostname (string), return the machine name (string)

    This uses some machine-specific rules. The logic here will need to
    be extended if there are other machines with special translation
    rules from hostname to machine name.
    """
    if re.match(r'cheyenne\d+', hostname):
        machine = 'cheyenne'
    else:
        machine = hostname

    return machine
