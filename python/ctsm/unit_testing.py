"""Functions to aid unit tests"""

from ctsm.ctsm_logging import setup_logging_for_tests

def setup_for_tests(enable_critical_logs=False):
    """Call this at the beginning of unit testing

    Does various setup that would normally be done by the top-level application/script

    Args:
    enable_critical_logs (bool): If True, then critical logging messages will be output;
        if False, then even critical messages will be suppressed
    """
    setup_logging_for_tests(enable_critical_logs)
