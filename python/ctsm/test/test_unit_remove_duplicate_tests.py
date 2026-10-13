#!/usr/bin/env python3
"""
Unit tests for subset_data

You can run this by:
    python -m unittest test_unit_subset_data.py
"""

import unittest
import os
import sys

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
import ctsm.remove_duplicate_tests as rdt
from ctsm.path_utils import path_to_ctsm_root

# pylint: disable=invalid-name


class TestRemoveDuplicateTests(unittest.TestCase):
    """
    Basic class for testing remove_duplicate_tests.py.
    """

    def setUp(self):
        self.testdata_dir = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "testinputs",
        )
        if not os.path.exists(self.testdata_dir):
            raise FileNotFoundError(self.testdata_dir)
        self.script = os.path.join(
            path_to_ctsm_root(),
            "python",
            "ctsm",
            "remove_duplicate_tests.py",
        )

    def test_removeduplicatetests_getargs_default(self):
        """
        Test that default input and output files are equal
        """
        sys.argv = [""]
        args = rdt.get_args()
        self.assertEqual(args.input_file, args.output_file)

    def test_removeduplicatetests_getargs_input(self):
        """
        Test that specified -i is read correctly
        """
        input_file = __file__
        sys.argv = [self.script, "-i", input_file]
        args = rdt.get_args()
        self.assertEqual(args.input_file, input_file)

    def test_removeduplicatetests_getargs_input_long(self):
        """
        Test that specified --input-file is read correctly
        """
        input_file = __file__
        sys.argv = [self.script, "--input-file", input_file]
        args = rdt.get_args()
        self.assertEqual(args.input_file, input_file)

    def test_removeduplicatetests_getargs_input_fails(self):
        """
        Test that missing input file triggers error
        """
        input_file = "wjefwrebrieugr.txt"
        sys.argv = [self.script, "-i", input_file]
        with self.assertRaisesRegex(
            FileNotFoundError, input_file
        ):
            rdt.get_args()

    def test_removeduplicatetests_getargs_output(self):
        """
        Test that specified -o is read correctly
        """
        output_file = "abc.txt"
        sys.argv = [self.script, "-o", output_file]
        args = rdt.get_args()
        self.assertEqual(args.output_file, output_file)

    def test_removeduplicatetests_getargs_output_long(self):
        """
        Test that specified --output-file is read correctly
        """
        output_file = "abc.txt"
        sys.argv = [self.script, "--output-file", output_file]
        args = rdt.get_args()
        self.assertEqual(args.output_file, output_file)

    def test_removeduplicatetests_frontmatter(self):
        """
        Test that front matter is read correctly
        """
        file_in = os.path.join(self.testdata_dir, "testlist_clm_target.xml")
        front_matter = rdt.get_front_matter(file_in)
        self.assertEqual(len(front_matter), 5)
        self.assertEqual(front_matter[-1], "-->\n")

    def test_removeduplicatetests_parsexml(self):
        """
        Test that parse_xml() doesn't error
        """
        file_in = os.path.join(self.testdata_dir, "testlist_clm_target.xml")
        rdt.parse_xml(file_in)

    def test_removeduplicatetests_reparse(self):
        """
        Test that reparse() doesn't error
        """
        file_in = os.path.join(self.testdata_dir, "testlist_clm_target.xml")
        tree = rdt.parse_xml(file_in)
        rdt.reparse(tree)

    def test_removeduplicatetests_xmltotext(self):
        """
        Test that xml_to_text() works correctly
        """
        file_in = os.path.join(self.testdata_dir, "testlist_clm_target.xml")
        tree = rdt.parse_xml(file_in)
        result = rdt.xml_to_text(tree)
        print(result)
        self.assertTrue(result.startswith("<testlist version="))
        self.assertTrue(result.endswith("</testlist>"))

if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
