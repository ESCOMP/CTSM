#!/usr/bin/env python3

"""Unit tests for add_super_testlists"""

import argparse
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from ctsm import unit_testing
from ctsm.add_super_testlists import (
    get_parser,
    process_and_check_args,
    make_new_machine_line_for_supertestlist,
    process_machines_block,
    add_super_testlist,
)

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestGetParser(unittest.TestCase):
    """Tests of add_super_testlists: get_parser"""

    def test_get_parser_defaults(self):
        """Tests that get_parser gives the expected defaults when no args are given"""
        parser = get_parser()
        args = parser.parse_args([])
        self.assertEqual(args.testlist_file, "testlist_clm.xml")
        self.assertIsNone(args.output_file)
        self.assertFalse(args.check_only)

    def test_get_parser_custom_args(self):
        """Tests that get_parser correctly parses explicitly-given arguments"""
        parser = get_parser()
        args = parser.parse_args(
            ["--testlist-file", "mylist.xml", "--output", "out.xml", "--check-only"]
        )
        self.assertEqual(args.testlist_file, "mylist.xml")
        self.assertEqual(args.output_file, "out.xml")
        self.assertTrue(args.check_only)


class TestProcessAndCheckArgs(unittest.TestCase):
    """Tests of add_super_testlists: process_and_check_args"""

    def setUp(self):
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_process_and_check_args_missing_file(self):
        """Tests that process_and_check_args aborts if testlist_file doesn't exist"""
        args = argparse.Namespace(
            testlist_file=os.path.join(self._testdir, "does_not_exist.xml"),
            output_file=None,
        )
        with self.assertRaisesRegex(SystemExit, "testlist file not found"):
            process_and_check_args(args)

    def test_process_and_check_args_default_output(self):
        """Tests that output_file defaults to '<testlist-file>.modified'"""
        testlist_path = os.path.join(self._testdir, "testlist_clm.xml")
        # pylint: disable=consider-using-with,unspecified-encoding
        open(testlist_path, "x").close()

        args = argparse.Namespace(testlist_file=testlist_path, output_file=None)
        args = process_and_check_args(args)

        self.assertEqual(args.testlist_file, Path(testlist_path))
        self.assertEqual(args.output_file, Path(testlist_path + ".modified"))

    def test_process_and_check_args_custom_output(self):
        """Tests that an explicitly-given output_file is respected"""
        testlist_path = os.path.join(self._testdir, "testlist_clm.xml")
        # pylint: disable=consider-using-with,unspecified-encoding
        open(testlist_path, "x").close()
        output_path = os.path.join(self._testdir, "custom_output.xml")

        args = argparse.Namespace(testlist_file=testlist_path, output_file=output_path)
        args = process_and_check_args(args)

        self.assertEqual(args.output_file, Path(output_path))


class TestMakeNewMachineLineForSupertestlist(unittest.TestCase):
    """Tests of add_super_testlists: make_new_machine_line_for_supertestlist"""

    def test_make_new_machine_line_basic(self):
        """Tests that a new machine line is built correctly, defaulting to ctsm_release"""
        attrs = ' name="derecho" compiler="intel" category="aux_clm"'
        indent = "      "
        result = make_new_machine_line_for_supertestlist(attrs, indent)
        expected = '      <machine name="derecho" compiler="intel" category="ctsm_release"/>'
        self.assertEqual(result, expected)

    def test_make_new_machine_line_custom_super_testlist(self):
        """Tests that a new machine line uses the given super_testlist as its category"""
        attrs = ' name="izumi" compiler="gnu" category="aux_clm"'
        indent = "  "
        result = make_new_machine_line_for_supertestlist(attrs, indent, super_testlist="ctsm_sci")
        self.assertIn('category="ctsm_sci"', result)
        self.assertIn('name="izumi"', result)
        self.assertIn('compiler="gnu"', result)

    def test_make_new_machine_line_missing_name(self):
        """Tests that missing a name attribute aborts"""
        attrs = ' compiler="intel" category="aux_clm"'
        with self.assertRaisesRegex(SystemExit, "doesn't have a name, compiler or category"):
            make_new_machine_line_for_supertestlist(attrs, "  ")

    def test_make_new_machine_line_missing_compiler(self):
        """Tests that missing a compiler attribute aborts"""
        attrs = ' name="derecho" category="aux_clm"'
        with self.assertRaisesRegex(SystemExit, "doesn't have a name, compiler or category"):
            make_new_machine_line_for_supertestlist(attrs, "  ")


class TestProcessMachinesBlock(unittest.TestCase):
    """Tests of add_super_testlists: process_machines_block"""

    def test_process_machines_block_already_present(self):
        """Tests that a block already containing super_testlist is returned unchanged"""
        block = (
            '\n      <machine name="derecho" compiler="intel" category="aux_clm"/>'
            '\n      <machine name="derecho" compiler="intel" category="ctsm_release"/>\n    '
        )
        result = process_machines_block(block, super_testlist="ctsm_release")
        self.assertEqual(result, block)

    def test_process_machines_block_no_machine_found(self):
        """Tests that a block with no <machine/> line aborts"""
        block = "\n      <notamachine/>\n    "
        with self.assertRaisesRegex(SystemExit, "block doesn't have a machine line"):
            process_machines_block(block, super_testlist="ctsm_release")

    def test_process_machines_block_duplicate_machines(self):
        """Tests that duplicated machine lines in a block abort"""
        line = '\n      <machine name="derecho" compiler="intel"/>'
        block = line + line + "\n    "
        with self.assertRaisesRegex(SystemExit, "block has duplicated tests"):
            process_machines_block(block, super_testlist="ctsm_release")

    def test_process_machines_block_testlist_not_present(self):
        """Tests that the block is left unchanged if the given testlist isn't in it"""
        block = '\n      <machine name="derecho" compiler="intel" category="aux_clm"/>\n    '
        result = process_machines_block(block, testlist="hillslope", super_testlist="ctsm_sci")
        self.assertEqual(result, block)

    def test_process_machines_block_not_in_testlist_present(self):
        """Tests that the block is left unchanged if not_in_testlist is present"""
        block = '\n      <machine name="derecho" compiler="intel" category="aux_clm"/>\n    '
        result = process_machines_block(
            block, not_in_testlist="aux_clm", super_testlist="ctsm_sci"
        )
        self.assertEqual(result, block)

    def test_process_machines_block_adds_missing_entry(self):
        """Tests that a block missing super_testlist gets a new <machine/> line added"""

        block = '\n      <machine name="derecho" compiler="intel" category="aux_clm"/>\n    '
        result = process_machines_block(block, super_testlist="ctsm_release")
        print( result )
        self.assertIn('category="ctsm_release"', result)
        self.assertIn('category="aux_clm"', result)


class TestAddSuperTestlist(unittest.TestCase):
    """Tests of add_super_testlists: add_super_testlist"""

    def test_add_super_testlist_no_change_needed(self):
        """Tests that text is returned unchanged when every block already has super_testlist"""
        xml = (
            "<testlist>\n"
            '  <test name="foo">\n'
            "    <machines>\n"
            '      <machine name="derecho" compiler="intel" category="ctsm_release"/>\n'
            "    </machines>\n"
            "  </test>\n"
            "</testlist>\n"
        )
        result = add_super_testlist(xml, "ctsm_release")
        self.assertEqual(result, xml)

    def test_add_super_testlist_adds_entry(self):
        """Tests that a missing super_testlist entry gets added to a <machines> block"""

        xml = (
            "<testlist>\n"
            '  <test name="foo">\n'
            "    <machines>\n"
            '      <machine name="derecho" compiler="intel" category="aux_clm"/>\n'
            "    </machines>\n"
            "  </test>\n"
            "</testlist>\n"
        )
        result = add_super_testlist(xml, "ctsm_release")
        self.assertIn('category="ctsm_release"', result)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
