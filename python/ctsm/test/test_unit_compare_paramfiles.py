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
# pylint: disable=too-many-lines


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

        # Should have 2 lines: header + 1 attribute
        self.assertEqual(result.count("\n"), 2)

        self.assertIn("Attribute(s) present in File 0 but not File 1:", result)
        # Check that attr1 appears on the line after the header
        lines = result.split("\n")
        header_idx = next(i for i, line in enumerate(lines) if "File 0 but not File 1:" in line)
        self.assertIn("attr1: value1", lines[header_idx + 1])

    def test_attribute_only_in_second(self):
        """Test when attribute exists only in second DataArray"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value1", "attr2": "value2"})

        result = cp._compare_attrs(da0, da1, "")

        # Should have 2 lines: header + 1 attribute
        self.assertEqual(result.count("\n"), 2)

        self.assertIn("Attribute(s) present in File 1 but not File 0:", result)
        # Check that attr2 appears on the line after the header
        lines = result.split("\n")
        header_idx = next(i for i, line in enumerate(lines) if "File 1 but not File 0:" in line)
        self.assertIn("attr2: value2", lines[header_idx + 1])

    def test_different_attribute_values(self):
        """Test when attributes have different values"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": "value2"})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        # Should have 3 lines: header + 2 attribute value lines
        self.assertEqual(result.count("\n"), 3)

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

        # Should have 3 lines: header + 2 attribute value lines
        self.assertEqual(result.count("\n"), 3)

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
        print(result)
        lines = result.split("\n")

        # Should have 3 lines: header + 2 attribute value lines
        self.assertEqual(result.count("\n"), 3)

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("valid_range", lines[header_idx + 1])

    def test_array_attribute_values_different_sizes(self):
        """Test when array attributes have different sizes"""

        da0 = xr.DataArray([1, 2, 3], attrs={"valid_range": np.array([0, 100])})
        da1 = xr.DataArray([4, 5, 6], attrs={"valid_range": np.array([0, 100, 200])})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        # Should have 3 lines: header + 2 attribute value lines
        self.assertEqual(result.count("\n"), 3)

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("valid_range", lines[header_idx + 1])

    def test_multiple_differences(self):
        """Test when there are multiple types of differences"""
        da0 = xr.DataArray([1, 2, 3], attrs={"attr1": "value1", "attr2": "value2", "attr3": "same"})
        da1 = xr.DataArray([4, 5, 6], attrs={"attr2": "different", "attr3": "same", "attr4": "new"})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        # Should have 7 lines: header + 6 attribute lines
        self.assertEqual(result.count("\n"), 7)

        # Should report attr1 only in file 0
        self.assertIn("Attribute(s) present in File 0 but not File 1:", result)
        header_idx = next(i for i, line in enumerate(lines) if "File 0 but not File 1:" in line)
        self.assertIn("attr1: value1", lines[header_idx + 1])

        # Should report attr4 only in file 1
        self.assertIn("Attribute(s) present in File 1 but not File 0:", result)
        header_idx = next(i for i, line in enumerate(lines) if "File 1 but not File 0:" in line)
        self.assertIn("attr4: new", lines[header_idx + 1])

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

        # Should have 3 lines: header + 2 attribute value lines
        self.assertEqual(result.count("\n"), 3)

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
        da1 = xr.DataArray([4, 5, 6], attrs={"attr1": [1, 2]})

        result = cp._compare_attrs(da0, da1, "")
        lines = result.split("\n")

        # Should have 3 lines: header + 2 attribute value lines
        self.assertEqual(result.count("\n"), 3)

        self.assertIn("Attribute(s) with different values:", result)
        header_idx = next(i for i, line in enumerate(lines) if "different values:" in line)
        self.assertIn("attr1", lines[header_idx + 1])


class TestOneUnequalValueMsg(unittest.TestCase):
    """Unit tests for _one_unequal_value_msg"""

    def test_raw_differ_ms_same(self):
        """Test when raw values differ but masked/scaled values are the same"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([1.0])  # After scaling
        np1_ms = np.array([1.0])  # After scaling
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("raw: 100 → 200 (but both 1.0 after masking/scaling)", result)

    def test_raw_same_ms_differ(self):
        """Test when raw values are the same but masked/scaled values differ"""
        np0 = np.array([100])
        np1 = np.array([100])
        np0_ms = np.array([1.0])
        np1_ms = np.array([2.0])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("masked/scaled (raw both 100): 1.0 → 2.0", result)

    def test_both_differ_no_scaling(self):
        """Test when both raw and m/s differ, but no scaling applied (values identical)"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([100])  # No scaling
        np1_ms = np.array([200])  # No scaling
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("raw and masked/scaled: 100 → 200", result)

    def test_both_differ_with_scaling(self):
        """Test when both raw and m/s differ, with different scaling"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([1.0])
        np1_ms = np.array([3.0])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 2 lines
        self.assertEqual(result.count("\n"), 2)

        lines = result.split("\n")
        # First line should have raw values
        self.assertIn("raw:           100 → 200", lines[0])
        # Second line should have masked/scaled values
        self.assertIn("masked/scaled: 1.0 → 3.0", lines[1])

    def test_single_element_array_no_indices(self):
        """Test with single element array (no indices in output)"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([100])
        np1_ms = np.array([200])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        # Should not have indices list for single element
        self.assertNotIn("[0]", result)
        self.assertIn("raw and masked/scaled: 100 → 200", result)

    def test_multi_element_array_with_indices(self):
        """Test with multi-element array (indices shown in output)"""
        np0 = np.array([100, 200, 300])
        np1 = np.array([100, 250, 300])
        np0_ms = np.array([100, 200, 300])
        np1_ms = np.array([100, 250, 300])
        indices = (1,)

        result = cp._one_unequal_value_msg(
            np0=np0,
            np1=np1,
            np0_ms=np0_ms,
            np1_ms=np1_ms,
            indices=indices,
            msg="",
            dimnames=["level"],
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        # Should have indices list for multi-element array
        # Check that [level 1] appears at the beginning (after indentation) immediately before the
        # message
        self.assertIn("[level 1] raw and masked/scaled: 200 → 250", result)

    def test_multidimensional_array_indices(self):
        """Test with multidimensional array"""
        np0 = np.array([[1, 2], [3, 4]])
        np1 = np.array([[1, 2], [3, 5]])
        np0_ms = np.array([[1, 2], [3, 4]])
        np1_ms = np.array([[1, 2], [3, 5]])
        indices = (1, 1)
        dimnames = ["dim0", "dim1"]

        result = cp._one_unequal_value_msg(
            np0=np0,
            np1=np1,
            np0_ms=np0_ms,
            np1_ms=np1_ms,
            indices=indices,
            msg="",
            dimnames=dimnames,
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        # Should show both indices immediately before the message
        self.assertIn("[dim0 1, dim1 1] raw and masked/scaled: 4 → 5", result)

    def test_nan_values_both_nan(self):
        """Test when both masked/scaled values are NaN"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([np.nan])
        np1_ms = np.array([np.nan])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        # Both NaN should be considered equal for m/s
        self.assertIn("raw: 100 → 200 (but both nan after masking/scaling)", result)

    def test_nan_value_one_side(self):
        """Test when only one masked/scaled value is NaN"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([np.nan])
        np1_ms = np.array([2.0])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 2 lines
        self.assertEqual(result.count("\n"), 2)

        lines = result.split("\n")
        # NaN != 2.0, so both should differ with separate lines
        self.assertIn("raw:", lines[0])
        self.assertIn("100 → 200", lines[0])
        self.assertIn("masked/scaled:", lines[1])
        self.assertIn("nan → 2.0", lines[1])

    def test_negative_values(self):
        """Test with negative values"""
        np0 = np.array([-100])
        np1 = np.array([-200])
        np0_ms = np.array([-100])
        np1_ms = np.array([-200])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("raw and masked/scaled: -100 → -200", result)

    def test_float_values(self):
        """Test with floating point values"""
        np0 = np.array([1.5])
        np1 = np.array([2.5])
        np0_ms = np.array([1.5])
        np1_ms = np.array([2.5])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("raw and masked/scaled: 1.5 → 2.5", result)

    def test_zero_values(self):
        """Test with zero values"""
        np0 = np.array([0])
        np1 = np.array([1])
        np0_ms = np.array([0])
        np1_ms = np.array([1])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("raw and masked/scaled: 0 → 1", result)

    def test_very_large_values(self):
        """Test with very large values"""
        np0 = np.array([1e10])
        np1 = np.array([2e10])
        np0_ms = np.array([1e10])
        np1_ms = np.array([2e10])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("10000000000.0 → 20000000000.0", result)

    def test_very_small_values(self):
        """Test with very small values"""
        np0 = np.array([1e-10])
        np1 = np.array([2e-10])
        np0_ms = np.array([1e-10])
        np1_ms = np.array([2e-10])
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        # Check that both parts are on the same line
        self.assertIn("raw and masked/scaled: 1e-10 → 2e-10", result)

    def test_appends_to_existing_message(self):
        """Test that function appends to existing message"""
        np0 = np.array([100])
        np1 = np.array([200])
        np0_ms = np.array([100])
        np1_ms = np.array([200])
        indices = (0,)
        existing_msg = "Previous content\n"

        result = cp._one_unequal_value_msg(
            np0=np0,
            np1=np1,
            np0_ms=np0_ms,
            np1_ms=np1_ms,
            indices=indices,
            msg=existing_msg,
            dimnames=None,
        )

        # Should have 2 lines
        self.assertEqual(result.count("\n"), 2)

        self.assertTrue(result.startswith("Previous content\n"))
        self.assertIn("raw and masked/scaled: 100 → 200", result)

    def test_integer_dtype(self):
        """Test with integer data type"""
        np0 = np.array([100], dtype=np.int32)
        np1 = np.array([200], dtype=np.int32)
        np0_ms = np.array([100], dtype=np.int32)
        np1_ms = np.array([200], dtype=np.int32)
        indices = (0,)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        self.assertIn("raw and masked/scaled: 100 → 200", result)

    def test_raises_error_when_values_equal(self):
        """Test that RuntimeError is raised when values are actually equal"""
        np0 = np.array([100])
        np1 = np.array([100])
        np0_ms = np.array([100])
        np1_ms = np.array([100])
        indices = (0,)

        with self.assertRaises(RuntimeError):
            cp._one_unequal_value_msg(
                np0=np0,
                np1=np1,
                np0_ms=np0_ms,
                np1_ms=np1_ms,
                indices=indices,
                msg="",
                dimnames=None,
            )

    def test_indices_alignment_in_two_line_output(self):
        """Test that indices are properly aligned in two-line output"""
        np0 = np.array([100, 200, 300])
        np1 = np.array([100, 250, 300])
        np0_ms = np.array([1.0, 2.0, 3.0])
        np1_ms = np.array([1.0, 2.5, 3.0])
        indices = (1,)

        result = cp._one_unequal_value_msg(
            np0=np0,
            np1=np1,
            np0_ms=np0_ms,
            np1_ms=np1_ms,
            indices=indices,
            msg="",
            dimnames=["some_dim"],
        )

        # Should have 2 lines
        self.assertEqual(result.count("\n"), 2)

        lines = result.split("\n")
        # First line should have indices and raw values
        self.assertIn("[some_dim 1] raw:           200 → 250", lines[0])
        # Second line should have spaces equal to "[1] " (4 chars) to align
        self.assertIn("             masked/scaled: 2.0 → 2.5", lines[1])

    def test_string_dtype(self):
        """Test with string data type"""
        np0 = np.array(["apple", "banana", "cherry"])
        np1 = np.array(["apple", "orange", "cherry"])
        np0_ms = np.array(["apple", "banana", "cherry"])
        np1_ms = np.array(["apple", "orange", "cherry"])
        indices = (1,)

        # Verify that equal_nan=True raises TypeError for string arrays
        with self.assertRaises(TypeError):
            np.array_equal(np0, np1, equal_nan=True)

        result = cp._one_unequal_value_msg(
            np0=np0, np1=np1, np0_ms=np0_ms, np1_ms=np1_ms, indices=indices, msg="", dimnames=None
        )

        # Should have 1 line
        self.assertEqual(result.count("\n"), 1)

        # Should handle string comparison without NaN issues
        self.assertIn("[1] raw and masked/scaled: banana → orange", result)


class TestCompareDaValues(unittest.TestCase):
    """Unit tests for _compare_da_values"""

    def test_identical_values(self):
        """Test when all values are identical"""
        da0 = xr.DataArray([1, 2, 3])
        da1 = xr.DataArray([1, 2, 3])
        da0_ms = xr.DataArray([1, 2, 3])
        da1_ms = xr.DataArray([1, 2, 3])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        self.assertEqual(result, "")

    def test_single_difference(self):
        """Test when there is a single value difference"""
        da0 = xr.DataArray([1, 2, 3])
        da1 = xr.DataArray([1, 5, 3])
        da0_ms = xr.DataArray([1, 2, 3])
        da1_ms = xr.DataArray([1, 5, 3])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 2 lines including header
        self.assertEqual(result.count("\n"), 2)

        self.assertIn("Values differ:", result)
        self.assertIn("2 → 5", result)

    def test_multiple_differences(self):
        """Test when there are multiple value differences"""
        da0 = xr.DataArray([1, 2, 3, 4])
        da1 = xr.DataArray([1, 5, 3, 8])
        da0_ms = xr.DataArray([1, 2, 3, 4])
        da1_ms = xr.DataArray([1, 5, 3, 8])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 3 lines including header
        self.assertEqual(result.count("\n"), 3)

        self.assertIn("Values differ:", result)
        # Should report both differences
        self.assertIn("2 → 5", result)
        self.assertIn("4 → 8", result)

    def test_scalar_values_identical(self):
        """Test with scalar values that are identical"""
        da0 = xr.DataArray(42)
        da1 = xr.DataArray(42)
        da0_ms = xr.DataArray(42)
        da1_ms = xr.DataArray(42)

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        self.assertEqual(result, "")

    def test_scalar_values_differ(self):
        """Test with scalar values that differ"""
        da0 = xr.DataArray(42)
        da1 = xr.DataArray(99)
        da0_ms = xr.DataArray(42)
        da1_ms = xr.DataArray(99)

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 2 lines including header
        self.assertEqual(result.count("\n"), 2)

        self.assertIn("Values differ:", result)
        self.assertIn("42 → 99", result)

    def test_with_nans_both_nan(self):
        """Test when both masked/scaled arrays have NaN at same position"""
        da0 = xr.DataArray([1.0, 2.0, 3.0])
        da1 = xr.DataArray([1.0, 5.0, 3.0])
        da0_ms = xr.DataArray([1.0, np.nan, 3.0])
        da1_ms = xr.DataArray([1.0, np.nan, 3.0])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 2 lines including header
        self.assertEqual(result.count("\n"), 2)

        # Raw values differ but m/s are same (both NaN)
        self.assertIn("Values differ:", result)
        self.assertIn("2.0 → 5.0", result)

    def test_with_nans_different(self):
        """Test when NaN appears in only one masked/scaled array"""
        da0 = xr.DataArray([1.0, 2.0, 3.0])
        da1 = xr.DataArray([1.0, 5.0, 3.0])  # Raw value differs
        da0_ms = xr.DataArray([1.0, np.nan, 3.0])
        da1_ms = xr.DataArray([1.0, 5.0, 3.0])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 3 lines including header
        self.assertEqual(result.count("\n"), 3)

        # Raw values differ and m/s differ
        self.assertIn("Values differ:", result)
        self.assertIn("2.0 → 5.0", result)
        self.assertIn("nan → 5.0", result)

    def test_string_dtype_no_nan_handling(self):
        """Test with string dtype (can't have NaN, should handle TypeError)"""
        da0 = xr.DataArray(np.array(["a", "b", "c"], dtype=str))
        da1 = xr.DataArray(np.array(["a", "x", "c"], dtype=str))
        # Verify that equal_nan=True raises TypeError for string arrays
        with self.assertRaises(TypeError):
            np.array_equal(da0, da1, equal_nan=True)
        da0_ms = xr.DataArray(np.array(["a", "b", "c"], dtype=str))
        da1_ms = xr.DataArray(np.array(["a", "x", "c"], dtype=str))

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 2 lines including header
        self.assertEqual(result.count("\n"), 2)

        # Should still work despite TypeError in equal_nan
        self.assertIn("Values differ:", result)
        self.assertIn("b → x", result)

    def test_appends_to_existing_message(self):
        """Test that function appends to existing message"""
        da0 = xr.DataArray([1, 2])
        da1 = xr.DataArray([1, 5])
        da0_ms = xr.DataArray([1, 2])
        da1_ms = xr.DataArray([1, 5])
        existing_msg = "Previous content\n"

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, existing_msg)

        # Should have 3 lines including header and previous content
        self.assertEqual(result.count("\n"), 3)

        self.assertTrue(result.startswith("Previous content\n"))
        self.assertIn("Values differ:", result)

    def test_multidimensional_array(self):
        """Test with multidimensional arrays"""
        dimnames = ["level", "cohort"]
        da0 = xr.DataArray([[1, 2], [3, 4]], dims=dimnames)
        da1 = xr.DataArray([[1, 2], [3, 9]], dims=dimnames)
        da0_ms = xr.DataArray([[1, 2], [3, 4]], dims=dimnames)
        da1_ms = xr.DataArray([[1, 2], [3, 9]], dims=dimnames)

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 2 lines including header
        self.assertEqual(result.count("\n"), 2)

        self.assertIn("Values differ:", result)
        self.assertIn("[level 1, cohort 1]", result)
        self.assertIn("4 → 9", result)

    def test_multidimensional_array_dimnames_dont_match(self):
        """Test with multidimensional arrays whose dimension names don't match"""
        da0 = xr.DataArray([[1, 2], [3, 4]], dims=["level", "cohort"])
        da1 = xr.DataArray([[1, 2], [3, 9]], dims=["dim0", "dim1"])
        da0_ms = xr.DataArray([[1, 2], [3, 4]], dims=["level", "cohort"])
        da1_ms = xr.DataArray([[1, 2], [3, 9]], dims=["dim0", "dim1"])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 2 lines including header
        self.assertEqual(result.count("\n"), 2)

        self.assertIn("Values differ:", result)
        self.assertIn("[1, 1]", result)
        self.assertIn("4 → 9", result)

    def test_all_values_differ(self):
        """Test when all values differ"""
        da0 = xr.DataArray([1, 2, 3])
        da1 = xr.DataArray([4, 5, 6])
        da0_ms = xr.DataArray([1, 2, 3])
        da1_ms = xr.DataArray([4, 5, 6])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 4 lines including header
        self.assertEqual(result.count("\n"), 4)

        self.assertIn("Values differ:", result)
        # Should report all three differences
        self.assertIn("1 → 4", result)
        self.assertIn("2 → 5", result)
        self.assertIn("3 → 6", result)

    def test_values_differ_header_only_once(self):
        """Test that 'Values differ:' header appears only once even with multiple differences"""
        da0 = xr.DataArray([1, 2, 3])
        da1 = xr.DataArray([4, 5, 6])
        da0_ms = xr.DataArray([1, 2, 3])
        da1_ms = xr.DataArray([4, 5, 6])

        result = cp._compare_da_values(da0_ms, da1_ms, da0, da1, "")

        # Should have 4 lines including header
        self.assertEqual(result.count("\n"), 4)

        # Count occurrences of the header
        self.assertEqual(result.count("Values differ:"), 1)


if __name__ == "__main__":
    unittest.main()
