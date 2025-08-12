#!/usr/bin/env python3

"""System tests for query_paramfile"""

import unittest
import os
import sys
import io
from contextlib import redirect_stdout

from ctsm import unit_testing

from ctsm.param_utils import query_paramfile as qp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


class TestSysQueryParamfile(unittest.TestCase):
    """System tests of query_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv

    def tearDown(self):
        sys.argv = self.orig_argv

    def test_query_paramfile_scalar(self):
        """Test that print_values works with scalar parameter"""

        sys.argv = ["get_arguments", "-i", PARAMFILE, "phenology_soil_depth"]

        f = io.StringIO()
        with redirect_stdout(f):
            qp.main()
        out = f.getvalue()
        self.assertEqual("phenology_soil_depth: 0.08\n", out)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
