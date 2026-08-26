#!/usr/bin/env python3
"""
This is a just top-level skeleton script that calls
add_super_testlists.py.
The original code (add_super_testlists.py) is located under the
python/ctsm folder.

For full instructions on how to run the code and different options,
please check the python/ctsm/add_super_testlists.py file, or run this
script with --help.

----------------------------------------------------------------
Instructions for running using conda python environments:

../../py_env_create
conda activate ctsm_py

Utilities to ensure tests in `testlist_clm.xml` are included in
super-testlists (for example, `ctsm_release`).
"""
import os
import sys

# -- add python/ctsm to path
_CTSM_PYTHON = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir, "python"
)
sys.path.insert(1, _CTSM_PYTHON)

from ctsm.add_super_testlists import main  # pylint: disable=wrong-import-position

if __name__ == "__main__":
    main()
