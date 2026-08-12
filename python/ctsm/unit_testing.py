"""Functions to aid unit tests"""

from pathlib import Path
import sys
from ctsm.ctsm_logging import setup_logging_for_tests


def add_machine_node_args(machine, nodes, tasks):
    """add arguments to sys.argv"""
    args_to_add = [
        "--machine",
        machine,
        "--number-of-nodes",
        str(nodes),
        "--tasks-per-node",
        str(tasks),
    ]
    sys.argv += args_to_add


def setup_for_tests(enable_critical_logs=False):
    """Call this at the beginning of unit testing

    Does various setup that would normally be done by the top-level application/script

    Args:
    enable_critical_logs (bool): If True, then critical logging messages will be output;
        if False, then even critical messages will be suppressed
    """
    setup_logging_for_tests(enable_critical_logs)


def get_test_input_data_dir():
    """
    Get the absolute path to the directory containing Python unit/system test input data
    """
    test_input_data_dir = Path(__file__) / ".." / "test" / "testinputs"
    return str(test_input_data_dir.resolve())
