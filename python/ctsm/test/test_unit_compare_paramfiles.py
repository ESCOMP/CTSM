#!/usr/bin/env python3

"""Unit tests for compare_paramfiles"""

import unittest
import unittest.mock
import os
import sys
import shutil
import tempfile
import argparse
from io import StringIO
import xarray as xr
import numpy as np

from ctsm.param_utils import compare_paramfiles as cp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

# pylint: disable=protected-access


class TestGetArguments(unittest.TestCase):
    """Unit tests for get_arguments"""

    def test_get_arguments(self):
        """Test get_arguments correctly parses command line arguments"""
        file0 = "/path/to/first/file.nc"
        file1 = "/path/to/second/file.nc"
        sys.argv = ["compare_paramfiles.py", file0, file1]
        args = cp.get_arguments()

        self.assertEqual(args.file0, file0)
        self.assertEqual(args.file1, file1)
        self.assertIsInstance(args, argparse.Namespace)


class TestCheckArguments(unittest.TestCase):
    """Unit tests for check_arguments"""

    def setUp(self):
        """Set up test fixtures"""
        # Create temporary directory for test files
        self.tempdir = tempfile.mkdtemp()

        # Create two simple empty test files
        self.file0 = os.path.join(self.tempdir, "test_file0.nc")
        self.file1 = os.path.join(self.tempdir, "test_file1.nc")

        # Create empty files. Using "with" ensures release of the allocated resources even in the
        # case of an exception.
        with open(self.file0, "w", encoding="utf-8"):
            pass
        with open(self.file1, "w", encoding="utf-8"):
            pass

    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.tempdir)

    def test_check_arguments_valid(self):
        """Test check_arguments with valid files"""
        args = argparse.Namespace(file0=self.file0, file1=self.file1)
        # Should not raise any exception
        cp.check_arguments(args)

    def test_check_arguments_file0_missing(self):
        """Test check_arguments raises FileNotFoundError when file0 doesn't exist"""
        nonexistent = os.path.join(self.tempdir, "nonexistent.nc")
        args = argparse.Namespace(file0=nonexistent, file1=self.file1)

        with self.assertRaises(FileNotFoundError):
            cp.check_arguments(args)

    def test_check_arguments_file1_missing(self):
        """Test check_arguments raises FileNotFoundError when file1 doesn't exist"""
        nonexistent = os.path.join(self.tempdir, "nonexistent.nc")
        args = argparse.Namespace(file0=self.file0, file1=nonexistent)

        with self.assertRaises(FileNotFoundError):
            cp.check_arguments(args)

    def test_check_arguments_same_file(self):
        """Test check_arguments exits when files are the same"""
        args = argparse.Namespace(file0=self.file0, file1=self.file0)

        with self.assertRaises(SystemExit):
            with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                cp.check_arguments(args)
                self.assertIn("These are the same file.", mock_stdout.getvalue())

    def test_check_arguments_same_file_realpath(self):
        """Test check_arguments exits when files are the same via realpath"""
        # Create a symlink to file0
        symlink = os.path.join(self.tempdir, "symlink.nc")
        os.symlink(self.file0, symlink)

        args = argparse.Namespace(file0=self.file0, file1=symlink)

        with self.assertRaises(SystemExit):
            with unittest.mock.patch("sys.stdout", new_callable=StringIO) as mock_stdout:
                cp.check_arguments(args)
                self.assertIn("These are the same file.", mock_stdout.getvalue())


class TestGetVariablesInOnlyOneDs(unittest.TestCase):
    """Unit tests for _get_variables_in_only_one_ds"""

    def test_variables_in_only_first(self):
        """Test finding variables only in first dataset"""
        ds_a = xr.Dataset({"var1": xr.DataArray([1, 2, 3]), "var2": xr.DataArray([4, 5, 6])})
        ds_b = xr.Dataset({"var2": xr.DataArray([7, 8, 9])})

        result = cp._get_variables_in_only_one_ds(ds_a, ds_b)

        self.assertEqual(result, ["var1"])

    def test_no_unique_variables(self):
        """Test when all variables are shared"""
        ds_a = xr.Dataset({"var1": xr.DataArray([1, 2, 3])})
        ds_b = xr.Dataset({"var1": xr.DataArray([4, 5, 6])})

        result = cp._get_variables_in_only_one_ds(ds_a, ds_b)

        self.assertEqual(result, [])

    def test_all_variables_unique(self):
        """Test when all variables are unique to first dataset"""
        ds_a = xr.Dataset({"var1": xr.DataArray([1, 2, 3]), "var2": xr.DataArray([4, 5, 6])})
        ds_b = xr.Dataset({"var3": xr.DataArray([7, 8, 9])})

        result = cp._get_variables_in_only_one_ds(ds_a, ds_b)

        self.assertEqual(result, ["var1", "var2"])

    def test_sorted_output(self):
        """Test that output is sorted alphabetically"""
        ds_a = xr.Dataset(
            {
                "zebra": xr.DataArray([1]),
                "apple": xr.DataArray([2]),
                "middle": xr.DataArray([3]),
            }
        )
        ds_b = xr.Dataset({"other": xr.DataArray([4])})

        result = cp._get_variables_in_only_one_ds(ds_a, ds_b)

        self.assertEqual(result, ["apple", "middle", "zebra"])


class TestGetVariablesInBothDs(unittest.TestCase):
    """Unit tests for _get_variables_in_both_ds"""

    def test_variables_in_both(self):
        """Test finding variables in both datasets"""
        ds_a = xr.Dataset({"var1": xr.DataArray([1, 2, 3]), "var2": xr.DataArray([4, 5, 6])})
        ds_b = xr.Dataset({"var2": xr.DataArray([7, 8, 9]), "var3": xr.DataArray([10, 11, 12])})

        result = cp._get_variables_in_both_ds(ds_a, ds_b)

        self.assertEqual(result, ["var2"])

    def test_all_variables_shared(self):
        """Test when all variables are shared"""
        ds_a = xr.Dataset({"var1": xr.DataArray([1, 2, 3]), "var2": xr.DataArray([4, 5, 6])})
        ds_b = xr.Dataset({"var1": xr.DataArray([7, 8, 9]), "var2": xr.DataArray([10, 11, 12])})

        result = cp._get_variables_in_both_ds(ds_a, ds_b)

        self.assertEqual(result, ["var1", "var2"])

    def test_no_shared_variables(self):
        """Test when no variables are shared"""
        ds_a = xr.Dataset({"var1": xr.DataArray([1, 2, 3])})
        ds_b = xr.Dataset({"var2": xr.DataArray([4, 5, 6])})

        result = cp._get_variables_in_both_ds(ds_a, ds_b)

        self.assertEqual(result, [])

    def test_sorted_output(self):
        """Test that output is sorted alphabetically"""
        ds_a = xr.Dataset(
            {
                "zebra": xr.DataArray([1]),
                "apple": xr.DataArray([2]),
                "middle": xr.DataArray([3]),
            }
        )
        ds_b = xr.Dataset(
            {
                "zebra": xr.DataArray([4]),
                "apple": xr.DataArray([5]),
                "middle": xr.DataArray([6]),
            }
        )

        result = cp._get_variables_in_both_ds(ds_a, ds_b)

        self.assertEqual(result, ["apple", "middle", "zebra"])


class TestGetAttributesInOnlyOneDa(unittest.TestCase):
    """Unit tests for _get_attributes_in_only_one_da"""

    def test_attributes_in_only_first(self):
        """Test finding attributes only in first DataArray"""
        da_a = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2"})
        da_b = xr.DataArray([4, 5, 6], attrs={"attr2": "value2"})

        result = cp._get_attributes_in_only_one_da(da_a, da_b)

        self.assertEqual(result, ["attr1"])

    def test_no_unique_attributes(self):
        """Test when all attributes are shared"""
        da_a = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da_b = xr.DataArray([4, 5, 6], attrs={"attr1": "value1"})

        result = cp._get_attributes_in_only_one_da(da_a, da_b)

        self.assertEqual(result, [])

    def test_all_attributes_unique(self):
        """Test when all attributes are unique to first DataArray"""
        da_a = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2"})
        da_b = xr.DataArray([4, 5, 6], attrs={"attr3": "value3"})

        result = cp._get_attributes_in_only_one_da(da_a, da_b)

        self.assertEqual(result, ["attr1", "attr2"])

    def test_sorted_output(self):
        """Test that output is sorted alphabetically"""
        da_a = xr.DataArray([1], attrs={"zebra": "z", "apple": "a", "middle": "m"})
        da_b = xr.DataArray([2], attrs={"other": "o"})

        result = cp._get_attributes_in_only_one_da(da_a, da_b)

        self.assertEqual(result, ["apple", "middle", "zebra"])


class TestGetAttributesInBothDa(unittest.TestCase):
    """Unit tests for _get_attributes_in_both_da"""

    def test_attributes_in_both(self):
        """Test finding attributes in both DataArrays"""
        da_a = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2"})
        da_b = xr.DataArray([4, 5, 6], attrs={"attr2": "value2", "attr3": "value3"})

        result = cp._get_attributes_in_both_da(da_a, da_b)

        self.assertEqual(result, ["attr2"])

    def test_all_attributes_shared(self):
        """Test when all attributes are shared"""
        da_a = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2"})
        da_b = xr.DataArray([4, 5, 6], attrs={"attr1": "value1", "attr2": "value2"})

        result = cp._get_attributes_in_both_da(da_a, da_b)

        self.assertEqual(result, ["attr1", "attr2"])

    def test_no_shared_attributes(self):
        """Test when no attributes are shared"""
        da_a = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da_b = xr.DataArray([4, 5, 6], attrs={"attr2": "value2"})

        result = cp._get_attributes_in_both_da(da_a, da_b)

        self.assertEqual(result, [])

    def test_sorted_output(self):
        """Test that output is sorted alphabetically"""
        da_a = xr.DataArray([1], attrs={"zebra": "z", "apple": "a", "middle": "m"})
        da_b = xr.DataArray([2], attrs={"zebra": "z", "apple": "a", "middle": "m"})

        result = cp._get_attributes_in_both_da(da_a, da_b)

        self.assertEqual(result, ["apple", "middle", "zebra"])


class TestCompareAttrs(unittest.TestCase):
    """Unit tests for _compare_attrs"""

    def test_identical_attributes(self):
        """Test when DataArrays have identical attributes"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value1", "attr2": "value2"})

        result = cp._compare_attrs(da0, da1, "")

        self.assertEqual(result, "")

    def test_no_attributes(self):
        """Test when DataArrays have no attributes"""
        da0 = xr.DataArray([1, 2, 3])
        da1 = xr.DataArray([4, 5, 6])

        result = cp._compare_attrs(da0, da1, "")

        self.assertEqual(result, "")

    def test_attribute_only_in_first(self):
        """Test when attribute exists only in first DataArray"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr2": "value2"})

        result = cp._compare_attrs(da0, da1, "")

        self.assertIn("Attribute(s) present in File 0 but not File 1:", result)
        # Check that attr1 appears on the line after the header
        lines = result.split("\n")
        header_idx = next(i for i, line in enumerate(lines) if "File 0 but not File 1:" in line)
        self.assertIn("attr1", lines[header_idx + 1])

    def test_attribute_only_in_second(self):
        """Test when attribute exists only in second DataArray"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value1", "attr2": "value2"})

        result = cp._compare_attrs(da0, da1, "")

        self.assertIn("Attribute(s) present in File 1 but not File 0:", result)
        # Check that attr2 appears on the line after the header
        lines = result.split("\n")
        header_idx = next(i for i, line in enumerate(lines) if "File 1 but not File 0:" in line)
        self.assertIn("attr2", lines[header_idx + 1])

    def test_different_attribute_values(self):
        """Test when attributes have different values"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value2"})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("attr1, file 0: value1", lines[header_idx + 1])
        self.assertIn("attr1, file 1: value2", lines[header_idx + 2])

    def test_numeric_attribute_values(self):
        """Test when attributes have numeric values that differ"""
        da0 = xr.DataArray([1, 2, 3], attrs={"scale_factor": 1.0, "add_offset": 0.0})
        da1 = xr.DataArray([4, 5, 6], attrs={"scale_factor": 2.0, "add_offset": 0.0})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        # scale_factor should be reported (differs)
        self.assertIn("scale_factor", lines[header_idx + 1])
        # add_offset should not be in the result (same in both)
        self.assertNotIn("add_offset", result)

    def test_array_attribute_values_match(self):
        """Test when attributes are matching numpy arrays"""

        da0 = xr.DataArray([1, 2, 3], attrs={"valid_range": np.array([0, 100])})
        da1 = xr.DataArray([4, 5, 6], attrs={"valid_range": np.array([0, 100])})

        result = cp._compare_attrs(da0, da1, "")

        self.assertEqual(result, "")

    def test_array_attribute_values_differ(self):
        """Test when array attributes differ"""

        da0 = xr.DataArray([1, 2, 3], attrs={"valid_range": np.array([0, 100])})
        da1 = xr.DataArray([4, 5, 6], attrs={"valid_range": np.array([0, 200])})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("valid_range", lines[header_idx + 1])

    # TODO: Test arrays and lists of different sizes

    def test_multiple_differences(self):
        """Test when there are multiple types of differences"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2", "attr3": "same"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr2": "different", "attr3": "same", "attr4": "new"})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        # Should report attr1 only in file 0
        self.assertIn("Attribute(s) present in File 0 but not File 1:", result)
        header_idx = next(i for i, line in enumerate(lines) if "File 0 but not File 1:" in line)
        self.assertIn("attr1", lines[header_idx + 1])

        # Should report attr4 only in file 1
        self.assertIn("Attribute(s) present in File 1 but not File 0:", result)
        header_idx = next(i for i, line in enumerate(lines) if "File 1 but not File 0:" in line)
        self.assertIn("attr4", lines[header_idx + 1])

        # Should report attr2 with different values
        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        # attr2 should appear on the next line (file 0 value)
        self.assertIn("attr2, file 0: value2", lines[header_idx + 1])
        # and the line after that (file 1 value)
        self.assertIn("attr2, file 1: different", lines[header_idx + 2])

        # Should not report attr3 (same in both)
        self.assertNotIn("attr3", result)

    def test_appends_to_existing_message(self):
        """Test that function appends to existing message"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value2"})

        existing_msg = "Previous message\n"
        result = cp._compare_attrs(da0, da1, existing_msg)

        self.assertTrue(result.startswith("Previous message\n"))
        self.assertIn("Attribute(s) with different values:", result)

    def test_empty_string_attributes(self):
        """Test when attributes are empty strings"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": ""})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": ""})

        result = cp._compare_attrs(da0, da1, "")

        self.assertEqual(result, "")

    def test_empty_vs_nonempty_string(self):
        """Test when one attribute is empty string and other is not"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": ""})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value"})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("attr1, file 0:", lines[header_idx + 1])
        self.assertIn("attr1, file 1: value", lines[header_idx + 2])

    def test_none_attribute_values(self):
        """Test when attributes have None values"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": None})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": None})

        result = cp._compare_attrs(da0, da1, "")

        self.assertEqual(result, "")

    def test_list_attribute_values(self):
        """Test when attributes are lists"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": [1, 2, 3]})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": [1, 2, 3]})

        result = cp._compare_attrs(da0, da1, "")

        self.assertEqual(result, "")

    def test_list_attribute_values_differ(self):
        """Test when list attributes differ"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": [1, 2, 3]})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": [1, 2, 4]})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("attr1", lines[header_idx + 1])


if __name__ == "__main__":
    unittest.main()
