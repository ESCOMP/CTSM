"""Additions to the six library needed for python 2/3 compatibility"""

import six

try:
    # Only available in python 3.3+
    from unittest import mock
except ImportError:
    # Requires that mock be installed
    import mock
