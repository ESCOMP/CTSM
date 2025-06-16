#!/usr/bin/env python3

"""System tests for hillslope processing"""

import os

import unittest
import tempfile
import shutil
import sys
import glob

import xarray as xr

# -- add python/ctsm  to path (needed if we want to run test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)
# pylint: disable=wrong-import-position
from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.hillslopes.combine_gridcell_files import main as combine_gridcell_files
from ctsm.hillslopes.combine_chunk_files import main as combine_chunk_files

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestHillslopes(unittest.TestCase):
    """System tests for CTSM-specific hillslope scripts"""

    def setUp(self):
        # Where in the /testinputs directory are the input and comparison file(s)?
        testinputs_path = os.path.join(path_to_ctsm_root(), "python", "ctsm", "test", "testinputs")
        testinputs_hills_path = os.path.join(testinputs_path, "hillslopes")
        self._testinputs_hills_path = testinputs_hills_path

        # Make /_tempdir for use by these tests.
        self._previous_dir = os.getcwd()
        self._tempdir = tempfile.mkdtemp()

        # Obtain paths for the files being created in /_tempdir
        self._outdir = os.path.join(self._tempdir, "hillslopes")
        os.makedirs(self._outdir)

        # What is the surface dataset we're using?
        self._fsurdat = os.path.join(
            testinputs_path, "surfdata_5x5_amazon_hist_16pfts_CMIP6_2000_c231031.nc"
        )

    def tearDown(self):
        """
        Remove temporary directory
        """
        os.chdir(self._previous_dir)
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def _combine_gridcell_files(self, verbose, outdir):
        """
        Combines gridcell files into chunk files
        """

        # Set arguments
        sys.argv = [
            "combine_gridcell_files",
            "-i",
            self._fsurdat,
            "-d",
            os.path.join(self._testinputs_hills_path, "gridcell_files"),
            "-o",
            outdir,
        ]
        if verbose:
            sys.argv += ["-v"]

        # Call script
        combine_gridcell_files()

    def test_combine_gridcell_files(self):
        """
        Tests combine_gridcell_files
        """

        # Call script
        outdir_chunkfiles = os.path.join(self._outdir, "chunk_files_out")
        self._combine_gridcell_files(verbose=True, outdir=outdir_chunkfiles)

        # Make sure every output file matches its counterpart in the comparison directory
        comparison_dir = os.path.join(self._testinputs_hills_path, "chunk_files")
        output_files = glob.glob(os.path.join(outdir_chunkfiles, "*.nc"))
        output_files.sort()
        for output_file in output_files:
            comparison_file = os.path.join(
                comparison_dir,
                os.path.basename(output_file),
            )
            output_ds = xr.open_dataset(output_file)
            comparison_ds = xr.open_dataset(comparison_file)
            xr.testing.assert_equal(output_ds, comparison_ds)

        # Make sure the output directory contains every expected file
        comparison_files = glob.glob(os.path.join(comparison_dir, "*.nc"))
        comparison_files.sort()
        for comparison_file in comparison_files:
            output_file = os.path.join(
                outdir_chunkfiles,
                os.path.basename(comparison_file),
            )
            if not os.path.exists(output_file):
                raise RuntimeError(f"Expected output file not found: '{output_file}'")

    def test_combine_chunk_files(self):
        """
        Tests combine_chunk_files
        """

        # Call script
        output_file = os.path.join(self._outdir, "hillslope_data.nc")
        sys.argv = [
            "combine_chunk_files",
            "-i",
            self._fsurdat,
            "-d",
            os.path.join(self._testinputs_hills_path, "chunk_files"),
            "-o",
            output_file,
            "-v",
        ]
        combine_chunk_files()

        # Open output file
        output_ds = xr.open_dataset(output_file)

        # Open comparison file
        comparison_file = os.path.join(
            self._testinputs_hills_path,
            "hilldata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc",
        )
        comparison_ds = xr.open_dataset(comparison_file)
        comparison_ds = comparison_ds.drop_vars("AREA")

        # Compare
        xr.testing.assert_equal(output_ds, comparison_ds)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
