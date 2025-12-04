#!/usr/bin/env python3

"""System tests for compare_paramfiles"""

import unittest
import unittest.mock
import os
import sys
import shutil
import tempfile
from io import StringIO
import numpy as np
import xarray as xr

from ctsm.param_utils import compare_paramfiles as cp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

NETCDF_TYPE = "NETCDF4_CLASSIC"


class TestCompareParamfilesMain(unittest.TestCase):
    """System tests for compare_paramfiles.main()"""

    # pylint: disable=too-many-instance-attributes

    def setUp(self):
        """Set up test fixtures - create temporary directory and test paramfiles"""
        self.tempdir = tempfile.mkdtemp()

        # Create PFT names
        pft_names = np.array([b"tree1", b"grass1", b"crop1"])

        # Create first paramfile with various types of variables
        self.file0 = os.path.join(self.tempdir, "params0.nc")
        ds0 = xr.Dataset(
            {
                # Scalar variables
                "scalar_int": xr.DataArray(
                    42, attrs={"units": "unitless", "long_name": "A scalar integer"}
                ),
                "scalar_float": xr.DataArray(
                    3.14, attrs={"units": "m", "long_name": "A scalar float"}
                ),
                # 1D PFT-specific variables
                "pft_param1": xr.DataArray(
                    [1.0, 2.0, 3.0],
                    dims=["pft"],
                    attrs={
                        "units": "kg/m2",
                        "long_name": "PFT parameter 1",
                        "coordinates": "pftname",
                    },
                ),
                "pft_param2": xr.DataArray(
                    [10, 20, 30],
                    dims=["pft"],
                    attrs={
                        "units": "days",
                        "long_name": "PFT parameter 2",
                        "coordinates": "pftname",
                    },
                ),
                # PFT names coordinate
                "pftname": xr.DataArray(
                    pft_names, dims=["pft"], attrs={"long_name": "PFT names", "units": "unitless"}
                ),
            }
        )
        ds0.to_netcdf(self.file0, format=NETCDF_TYPE)

        # Create second paramfile (identical to first)
        self.file1 = os.path.join(self.tempdir, "params1.nc")
        ds1 = ds0.copy(deep=True)
        ds1.to_netcdf(self.file1, format=NETCDF_TYPE)

        # Create third paramfile with some differences
        self.file2 = os.path.join(self.tempdir, "params2.nc")
        ds2 = ds0.copy(deep=True)
        # Change a scalar value
        ds2["scalar_int"] = xr.DataArray(
            99, attrs={"units": "unitless", "long_name": "A scalar integer"}
        )
        # Change a PFT-specific value
        ds2["pft_param1"].values[1] = 5.0
        # Add a new variable
        ds2["new_var"] = xr.DataArray(123, attrs={"units": "m/s", "long_name": "New variable"})
        ds2.to_netcdf(self.file2, format=NETCDF_TYPE)

        # Create fourth paramfile with attribute differences
        self.file3 = os.path.join(self.tempdir, "params3.nc")
        ds3 = ds0.copy(deep=True)
        ds3["scalar_float"].attrs["units"] = "km"  # Different units
        ds3["pft_param1"].attrs["long_name"] = "Different description"  # Different long_name
        ds3.to_netcdf(self.file3, format=NETCDF_TYPE)

        # Create fifth paramfile with missing variable
        self.file4 = os.path.join(self.tempdir, "params4.nc")
        ds4 = ds0.copy(deep=True)
        ds4 = ds4.drop_vars("pft_param2")
        ds4.to_netcdf(self.file4, format=NETCDF_TYPE)

        # Create sixth paramfile with both attribute AND value differences in same variable
        self.file5 = os.path.join(self.tempdir, "params5.nc")
        ds5 = ds0.copy(deep=True)
        # Change both attribute and value for pft_param1
        ds5["pft_param1"].attrs["units"] = "g/m2"  # Different units
        ds5["pft_param1"].values[0] = 9.9  # Different value
        ds5.to_netcdf(self.file5, format=NETCDF_TYPE)

        # Create seventh paramfile (identical to first except file type)
        self.file7 = os.path.join(self.tempdir, "params6.nc")
        ds7 = ds0.copy(deep=True)
        different_format = "NETCDF3_64BIT"
        assert different_format != NETCDF_TYPE
        ds7.to_netcdf(self.file7, format=different_format)

        # Create 8th paramfile (like file2 but with another new variable)
        self.file8 = os.path.join(self.tempdir, "params7.nc")
        ds8 = ds2.copy(deep=True)
        # Add a new variable
        ds8["new_var2"] = xr.DataArray(123, attrs={"units": "m/s", "long_name": "New variable 2"})
        ds8.to_netcdf(self.file8, format=NETCDF_TYPE)

        # Create 9th paramfile (like file4 but with an extra variable)
        self.file9 = os.path.join(self.tempdir, "params8.nc")
        ds9 = ds4.copy(deep=True)
        # Add a new variable
        ds9["new_var2"] = xr.DataArray(123, attrs={"units": "m/s", "long_name": "New variable 2"})
        ds9.to_netcdf(self.file9, format=NETCDF_TYPE)

    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.tempdir)

    def test_identical_files_no_output(self):
        """Test that comparing identical files produces no difference output"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file1]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should print file names and "Identical" message
        self.assertIn(f"File 0: {self.file0}", output)
        self.assertIn(f"File 1: {self.file1}", output)
        self.assertIn("Files are identical.", output)
        self.assertEqual(output.count("\n"), 4)

        # Should not print any variable names (no differences)
        self.assertNotIn("scalar_int:", output)
        self.assertNotIn("pft_param1:", output)
        self.assertNotIn("Values differ:", output)

    def test_file_type_difference(self):
        """Test that comparing files that differ in netCDF type prints message"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file7]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()
        self.assertIn("File types differ: NETCDF4_CLASSIC → NETCDF3_64BIT_OFFSET", output)

    def test_scalar_value_difference(self):
        """Test that scalar value differences are reported"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file2]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should report scalar_int difference (no indices for scalar)
        self.assertIn("Values differ:", output)
        self.assertRegex(output, r"scalar_int:\n\s+Values differ:\n\s+42 → 99")

    def test_pft_value_difference_with_index_and_pftname(self):
        """Test that PFT-specific value differences show indices and PFT name"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file2]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should report pft_param1 difference with index and dimension name
        self.assertRegex(
            output,
            r"pft_param1:\n\s+Values differ:\n\s+\[pft 1 \(grass1\)\] 2\.0 → 5\.0",
        )

    def test_variable_missing_from_first_file(self):
        """Test that variables missing from first file are reported"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file2]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should report new_var only in file 1
        self.assertRegex(
            output,
            r"Variable\(s\) present in File 1 but not File 0:\n\s+new_var",
        )

    def test_variable_missing_from_second_file(self):
        """Test that variables missing from second file are reported"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file4]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should report pft_param2 only in file 0
        self.assertRegex(
            output,
            r"Variable\(s\) present in File 0 but not File 1:\n\s+pft_param2",
        )

    def test_variables_missing_from_both_files(self):
        """Test that variables missing from first AND second file are reported"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file9]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should report pft_param2 only in file 0
        self.assertRegex(
            output,
            r"Variable\(s\) present in File 0 but not File 1:\n\s+pft_param2",
        )

        # Should report new_var2 only in file 1
        self.assertRegex(
            output,
            r"Variable\(s\) present in File 1 but not File 0:\n\s+new_var2",
        )

    def test_attribute_differences(self):
        """Test that attribute differences are reported"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file3]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Should report scalar_float attribute difference
        self.assertRegex(
            output,
            r"scalar_float:"
            r"\s+Attribute\(s\) with different values:"
            r"\s+units, file 0: m"
            r"\s+units, file 1: km",
        )

        # Should report pft_param1 attribute difference
        self.assertRegex(
            output,
            r"pft_param1:"
            r"\s+Attribute\(s\) with different values:"
            r"\s+long_name, file 0: PFT parameter 1"
            r"\s+long_name, file 1: Different description",
        )

    def test_multiple_differences_in_one_variable(self):
        """Test that multiple types of differences in one variable are all reported"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file5]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # pft_param1 should have both attribute AND value differences
        self.assertRegex(
            output,
            r"pft_param1:"
            r"\s+Attribute\(s\) with different values:"
            r"\s+units, file 0: kg/m2"
            r"\s+units, file 1: g/m2"
            r"\s+Values differ:"
            r"\s+\[pft 0 \(tree1\)\] 1\.0 → 9\.9",
        )

    def test_file_paths_printed(self):
        """Test that file paths are printed at the start"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file1]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        lines = output.split("\n")
        # File paths should be in first few lines
        self.assertIn(f"File 0: {self.file0}", lines[0])
        self.assertIn(f"File 1: {self.file1}", lines[1])

    def test_same_file_returns_early(self):
        """Test that comparing a file to itself returns early with message"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file0]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()
            lines = output.split("\n")
            self.assertEqual(len(lines), 2)
            self.assertIn("These are the same file.", output)

    def test_nonexistent_file_raises_error(self):
        """Test that nonexistent file raises FileNotFoundError"""
        nonexistent = os.path.join(self.tempdir, "nonexistent.nc")
        sys.argv = ["compare_paramfiles.py", self.file0, nonexistent]

        with self.assertRaises(FileNotFoundError):
            cp.main()

    def test_variables_reported_in_order(self):
        """Test that variables with differences are reported in sorted order"""
        sys.argv = ["compare_paramfiles.py", self.file0, self.file2]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        # Find positions of variable names in output
        pos_pft_param1 = output.find("pft_param1:")
        pos_scalar_int = output.find("scalar_int:")

        # Both should be present
        self.assertNotEqual(pos_pft_param1, -1)
        self.assertNotEqual(pos_scalar_int, -1)

        # pft_param1 should come before scalar_int (alphabetical order)
        self.assertLess(pos_pft_param1, pos_scalar_int)

    def test_params_arg_in0not1_complete(self):
        """Test that --params arg works with a variable in first file and not second"""

        # Test with all option names
        opt_name_list = ["params", "param", "parameters"]
        for opt_name in opt_name_list:
            sys.argv = [
                "compare_paramfiles.py",
                f"--{opt_name}",
                "pft_param1,new_var",
                self.file8,
                self.file0,
            ]

            with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                cp.main()
                output = mock_stdout.getvalue()

            self.assertRegex(output, r"Variable\(s\) present in File 0 but not File 1:\n\s+new_var")
            self.assertNotIn("new_var2", output)  # Even though it's only in one file
            self.assertIn("pft_param1", output)
            self.assertNotIn("scalar_int", output)  # Even though it does differ
            self.assertNotIn("scalar_float", output)
            self.assertNotIn("pft_param2", output)

    def test_params_arg_in1not0(self):
        """Test that --params arg works with a variable in first file and not second"""
        sys.argv = [
            "compare_paramfiles.py",
            "--params",
            "pft_param1,new_var",
            self.file0,
            self.file8,
        ]

        with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
            cp.main()
            output = mock_stdout.getvalue()

        self.assertRegex(output, r"Variable\(s\) present in File 1 but not File 0:\n\s+new_var")

    def test_params_arg_keyerror(self):
        """Test that --params arg throws KeyError if nonexistent variables requested"""
        sys.argv = [
            "compare_paramfiles.py",
            "--params",
            "missing_var1,missing_var2",
            self.file0,
            self.file8,
        ]

        with self.assertRaisesRegex(KeyError, ": missing_var1, missing_var2"):
            cp.main()


if __name__ == "__main__":
    unittest.main()
