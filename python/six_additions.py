"""Additions to the six library needed for python 2/3 compatibility"""

import six

try:
    # Only available in python 3.3+
    from unittest import mock
except ImportError:
    # Requires that mock be installed
    import mock

if six.PY3:
    # This is only available in python3.2 and later, so this code won't
    # run with python3 versions prior to 3.2
    _assertNotRegex = "assertNotRegex"
else:
    _assertNotRegex = "assertNotRegexpMatches"

def assertNotRegex(self, *args, **kwargs):
    return getattr(self, _assertNotRegex)(*args, **kwargs)
