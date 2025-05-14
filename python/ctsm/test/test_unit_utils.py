#!/usr/bin/env python3

"""Unit tests for utils and config_utils"""

import tempfile
import shutil
import unittest
import os

from ctsm import unit_testing
from ctsm.utils import fill_template_file, ensure_iterable
from ctsm.config_utils import _handle_config_value

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestUtilsFillTemplateFile(unittest.TestCase):
    """Tests of utils: fill_template_file"""

    def setUp(self):
        self._previous_dir = os.getcwd()
        self._testdir = tempfile.mkdtemp()

    def tearDown(self):
        os.chdir(self._previous_dir)
        shutil.rmtree(self._testdir, ignore_errors=True)

    def test_fillTemplateFile_basic(self):
        """Basic test of fill_template_file"""
        template_path = os.path.join(self._testdir, "template.txt")
        final_path = os.path.join(self._testdir, "final.txt")
        template_contents = """\
Hello
$foo
Goodbye
$bar
"""
        with open(template_path, "w") as f:
            f.write(template_contents)

        fillins = {"foo": "aardvark", "bar": "zyzzyva"}
        fill_template_file(template_path, final_path, fillins)

        expected_final_text = """\
Hello
aardvark
Goodbye
zyzzyva
"""
        with open(final_path) as f:
            final_contents = f.read()

        self.assertEqual(final_contents, expected_final_text)


class TestUtilsHandleConfigValue(unittest.TestCase):
    """Test of utils: _handle_config_value"""

    def test_handleConfigValue_UnsetCantBeUnset(self):
        """
        Tests the handling of UNSET variable read in from a .cfg file
        for which can_be_unset = False
        """
        val = "UNSET"
        item = "varname_in_cfg_file"
        default = None
        is_list = False
        convert_to_type = None
        can_be_unset = False
        allowed_values = None
        errmsg = "Must set a value for .cfg file variable: {}".format(item)

        with self.assertRaisesRegex(SystemExit, errmsg):
            val = _handle_config_value(
                var=val,
                default=default,
                item=item,
                is_list=is_list,
                convert_to_type=convert_to_type,
                can_be_unset=can_be_unset,
                allowed_values=allowed_values,
            )

    def test_handleConfigValue_UnsetCanBeUnset(self):
        """
        Tests the handling of UNSET variable read in from a .cfg file
        for which can_be_unset = True
        """
        val = "UNSET"
        item = "varname_in_cfg_file"
        default = [True, False, True]
        is_list = True
        convert_to_type = None
        can_be_unset = True
        allowed_values = None

        val = _handle_config_value(
            var=val,
            default=default,
            item=item,
            is_list=is_list,
            convert_to_type=convert_to_type,
            can_be_unset=can_be_unset,
            allowed_values=allowed_values,
        )

        self.assertEqual(val, default)

    def test_handleConfigValue_convertToBoolFail(self):
        """
        Tests the handling of misspelled boolean read in from a .cfg file
        Also test whether the code can read a list of booleans
        """
        val = "False Tree False"  # intentionally misspelled True
        item = "varname_in_cfg_file"
        default = None
        is_list = True
        convert_to_type = bool
        can_be_unset = False
        allowed_values = None
        errmsg = "Non-boolean value found for .cfg file variable: {}".format(item)

        with self.assertRaisesRegex(SystemExit, errmsg):
            val = _handle_config_value(
                var=val,
                default=default,
                item=item,
                is_list=is_list,
                convert_to_type=convert_to_type,
                can_be_unset=can_be_unset,
                allowed_values=allowed_values,
            )

    def test_handleConfigValue_convertToBoolPass(self):
        """
        Tests the handling of boolean read in from a .cfg file
        Also test whether the code can read a list of booleans
        """
        val = "yes no"
        item = "varname_in_cfg_file"
        default = None
        is_list = True
        convert_to_type = bool
        can_be_unset = False
        allowed_values = None

        val = _handle_config_value(
            var=val,
            default=default,
            item=item,
            is_list=is_list,
            convert_to_type=convert_to_type,
            can_be_unset=can_be_unset,
            allowed_values=allowed_values,
        )

        self.assertTrue(val[0])
        self.assertFalse(val[1])

    def test_handleConfigValue_convertToTypePass(self):
        """
        Tests the handling of non-boolean list from a .cfg file
        """
        val = "-9 0.001"
        item = "varname_in_cfg_file"
        default = None
        is_list = True
        convert_to_type = float
        can_be_unset = False
        allowed_values = None

        val = _handle_config_value(
            var=val,
            default=default,
            item=item,
            is_list=is_list,
            convert_to_type=convert_to_type,
            can_be_unset=can_be_unset,
            allowed_values=allowed_values,
        )

        self.assertEqual(val[0], -9)
        self.assertEqual(val[1], 0.001)

    def test_handleConfigValue_convertToTypeFail(self):
        """
        Tests the handling of an incorrectly entered list from a .cfg file
        """
        val = "1 2 3 x 5 6 7"
        item = "varname_in_cfg_file"
        default = None
        is_list = True
        convert_to_type = float
        can_be_unset = False
        allowed_values = None
        errmsg = "Wrong type for .cfg file variable: {}".format(item)

        with self.assertRaisesRegex(SystemExit, errmsg):
            val = _handle_config_value(
                var=val,
                default=default,
                item=item,
                is_list=is_list,
                convert_to_type=convert_to_type,
                can_be_unset=can_be_unset,
                allowed_values=allowed_values,
            )

    def test_handleConfigValue_allowedValsFail(self):
        """
        Tests that the code aborts if val does not include all allowed_values
        """
        val = "1 2 3 4.5 6 7"
        item = "varname_in_cfg_file"
        default = None
        is_list = True
        convert_to_type = float
        can_be_unset = False
        allowed_values = [1, 2, 3, 4, 5, 6, 7]
        v = 4.5  # v must equal the misstyped value in val
        errmsg = "{} is not an allowed value for {} in .cfg file. Check allowed_values".format(
            v, item
        )

        with self.assertRaisesRegex(SystemExit, errmsg):
            val = _handle_config_value(
                var=val,
                default=default,
                item=item,
                is_list=is_list,
                convert_to_type=convert_to_type,
                can_be_unset=can_be_unset,
                allowed_values=allowed_values,
            )

    def test_handleConfigValue_isListFail(self):
        """
        Tests that the code aborts if we forget to set is_list = True
        """
        val = "True False"
        item = "varname_in_cfg_file"
        default = None
        is_list = False
        convert_to_type = bool
        can_be_unset = False
        allowed_values = None
        errmsg = "More than 1 element found for .cfg file variable: {}".format(item)

        with self.assertRaisesRegex(SystemExit, errmsg):
            val = _handle_config_value(
                var=val,
                default=default,
                item=item,
                is_list=is_list,
                convert_to_type=convert_to_type,
                can_be_unset=can_be_unset,
                allowed_values=allowed_values,
            )

    def test_handleConfigValue_isListFalse(self):
        """
        Tests that the code works for a basic case of is_list = False
        """
        val_in = "0.5"
        item = "varname_in_cfg_file"
        default = None
        is_list = False
        convert_to_type = float
        can_be_unset = False
        allowed_values = None

        val_out = _handle_config_value(
            var=val_in,
            default=default,
            item=item,
            is_list=is_list,
            convert_to_type=convert_to_type,
            can_be_unset=can_be_unset,
            allowed_values=allowed_values,
        )

        self.assertEqual(val_out, float(val_in))


class TestUtilsEnsureIterable(unittest.TestCase):
    """Tests of utils: ensure_iterable"""

    def test_ensure_iterable_number(self):
        """
        Tests that ensure_iterable(NUMBER, 3) produces a list of 3 NUMBERs
        """
        val = 724.1987
        self.assertEqual(ensure_iterable(val, 3), [val, val, val])

    def test_ensure_iterable_none(self):
        """
        Tests that ensure_iterable(None, 2) produces a list of 2 Nones
        """
        val = None
        self.assertEqual(ensure_iterable(val, 2), [val, val])

    def test_ensure_iterable_already(self):
        """
        Tests that ensure_iterable() returns the input if it's already iterable
        """
        val = [11, 12]
        self.assertEqual(ensure_iterable(val, 2), val)

    def test_ensure_iterable_error_wrong_length(self):
        """
        Tests that ensure_iterable() errors if input is iterable but of the wrong length
        """
        with self.assertRaisesRegex(ValueError, "Input is iterable but wrong length"):
            ensure_iterable([11, 12], 3)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
