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

from ctsm.param_utils import compare_paramfiles as cp

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


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


if __name__ == "__main__":
    unittest.main()
