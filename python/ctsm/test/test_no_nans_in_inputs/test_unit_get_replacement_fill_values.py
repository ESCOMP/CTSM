#!/usr/bin/env python3
"""
Unit tests for get_replacement_fill_values.py script.
"""

import os
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import (
    ERR_STR_SKIP_FILE,
    ERR_STR_SKIP_VAR,
    USER_REQ_DELETE,
    USER_REQ_QUIT,
    USER_REQ_SKIP_FILE,
    USER_REQ_SKIP_VAR,
)
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    FillValueConfig,
    VarContext,
    find_user_nl_files,
    get_fill_value_from_user,
    get_var_info,
    replace_env_vars_in_netcdf_paths,
    var_data_has_nan,
    VARSTARTS_TO_DEFAULT_NEG999,
)
from ctsm.no_nans_in_inputs import get_replacement_fill_values


# Test constants used in multiple tests
TEST_VAR_NAME = "test_var"
DEFAULT_VAR_CONTEXT = VarContext(var_name=TEST_VAR_NAME, target_type=float)
INT_VAR_CONTEXT = VarContext(var_name=TEST_VAR_NAME, target_type=int)


class TestGetFillValueFromUser:
    """Test the get_fill_value_from_user function."""

    def test_numeric_input(self, monkeypatch):
        """Test that a numeric input is converted to the target type."""
        expected = 3.14
        monkeypatch.setattr("builtins.input", lambda _: str(expected))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == expected
        assert isinstance(result, float)

    def test_integer_input(self, monkeypatch):
        """Test that an integer input is converted correctly."""
        expected = 42
        monkeypatch.setattr("builtins.input", lambda _: str(expected))
        result = get_fill_value_from_user(INT_VAR_CONTEXT, FillValueConfig())
        assert result == expected
        assert isinstance(result, int)

    def test_delete_command(self, monkeypatch):
        """Test that 'delete' returns USER_REQ_DELETE."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_DELETE)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig(allow_delete=True))
        assert result == USER_REQ_DELETE

    def test_delete_case_insensitive(self, monkeypatch):
        """Test that 'DELETE' is recognized case-insensitively."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_DELETE.upper())
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig(allow_delete=True))
        assert result == USER_REQ_DELETE

    def test_delete_not_allowed(self, monkeypatch):
        """Test that 'delete' is rejected when allow_delete is False, then accepts a number."""
        fallback_value = -1e20
        inputs = iter([USER_REQ_DELETE, str(fallback_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig(allow_delete=False))
        assert result == fallback_value

    def test_skip_variable(self, monkeypatch):
        """Test that 'skip' raises ValueError with SKIP_VARIABLE."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_SKIP_VAR)
        with pytest.raises(ValueError, match=ERR_STR_SKIP_VAR):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_skip_file(self, monkeypatch):
        """Test that 'skipfile' raises ValueError with SKIP_FILE."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_SKIP_FILE)
        with pytest.raises(ValueError, match=ERR_STR_SKIP_FILE):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_quit(self, monkeypatch):
        """Test that 'quit' raises KeyboardInterrupt."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_QUIT)
        with pytest.raises(KeyboardInterrupt):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_default_value_on_empty_input(self, monkeypatch):
        """Test that empty input uses the default value."""
        default = 1.5e36
        monkeypatch.setattr("builtins.input", lambda _: "")
        result = get_fill_value_from_user(
            DEFAULT_VAR_CONTEXT, FillValueConfig(default_value=default)
        )
        assert result == default

    def test_default_delete_on_empty_input(self, monkeypatch):
        """Test that empty input uses delete as default."""
        monkeypatch.setattr("builtins.input", lambda _: "")
        config = FillValueConfig(default_value=USER_REQ_DELETE)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, config)
        assert result == USER_REQ_DELETE

    def test_auto_delete_when_default_is_delete(self):
        """Test that delete_if_none_filled auto-deletes without prompting."""
        config = FillValueConfig(default_value=USER_REQ_DELETE, delete_if_none_filled=True)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, config)
        assert result == USER_REQ_DELETE

    def test_auto_delete_only_when_default_is_delete(self, monkeypatch):
        """Test that delete_if_none_filled does NOT auto-delete when default is not delete."""
        entered_value = 0.0
        default = -888.0
        monkeypatch.setattr("builtins.input", lambda _: str(entered_value))
        config = FillValueConfig(default_value=default, delete_if_none_filled=True)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, config)
        assert result == entered_value

    def test_nan_input_rejected(self, monkeypatch):
        """Test that NaN input is rejected, then accepts a valid number."""
        valid_value = 42.5
        inputs = iter(["nan", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_invalid_input_then_valid(self, monkeypatch):
        """Test that invalid input is rejected, then accepts a valid number."""
        valid_value = -777.0
        inputs = iter(["not_a_number", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_empty_input_no_default_then_valid(self, monkeypatch):
        """Test that empty input with no default shows help, then accepts a value."""
        valid_value = 100.0
        inputs = iter(["", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_whitespace_input_treated_as_empty(self, monkeypatch):
        """Test that whitespace-only input is treated as empty."""
        valid_value = -555.0
        inputs = iter(["   ", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_ctrl_c_twice_raises(self, monkeypatch):
        """Test that pressing Ctrl-C twice raises KeyboardInterrupt."""

        def raise_keyboard_interrupt(_):
            raise KeyboardInterrupt

        monkeypatch.setattr("builtins.input", raise_keyboard_interrupt)
        with pytest.raises(KeyboardInterrupt):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_ctrl_c_then_valid_input(self, monkeypatch):
        """Test that Ctrl-C followed by valid input works."""
        valid_value = 123.0
        # call_count is defined in the outer (test method) scope
        call_count = 0

        def input_with_ctrl_c(_):
            # nonlocal lets us modify call_count from the enclosing scope
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                # First call: simulate Ctrl-C. get_fill_value_from_user will
                # catch this, show ncdump output, and loop back to prompt again.
                raise KeyboardInterrupt
            # Second call: return valid input, which should be accepted
            return str(valid_value)

        monkeypatch.setattr("builtins.input", input_with_ctrl_c)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_ctrl_c_then_quit(self, monkeypatch):
        """Test that Ctrl-C followed by 'quit' raises KeyboardInterrupt."""
        call_count = 0

        def input_with_ctrl_c_then_quit(_):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise KeyboardInterrupt
            return USER_REQ_QUIT

        monkeypatch.setattr("builtins.input", input_with_ctrl_c_then_quit)
        with pytest.raises(KeyboardInterrupt):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_dry_run(self):
        """With --dry-run, user input should be unneeded; result should always be to skip var"""
        var_context = DEFAULT_VAR_CONTEXT
        var_context.dry_run = True
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == USER_REQ_SKIP_VAR


class TestVarDataHasNan:
    """Test the var_data_has_nan function."""

    def test_true_when_nan_present(self):
        """Test detection of NaN values."""
        da = xr.DataArray([[1.0, 2.0], [np.nan, 4.0]])
        assert var_data_has_nan(da)

    def test_false_when_no_nan(self):
        """Test correct response for arrays without NaN values."""
        da = xr.DataArray([1.0, 2.0, 3.0, 4.0])
        assert not var_data_has_nan(da)

    def test_true_when_all_nan(self):
        """Test that arrays with all NaN values are detected."""
        da = xr.DataArray([np.nan, np.nan, np.nan])
        assert var_data_has_nan(da)

    def test_false_for_empty_array(self):
        """Test that empty arrays return False (no NaNs)."""
        da = xr.DataArray([])
        assert not var_data_has_nan(da)


# pylint: disable=too-few-public-methods
class TestGetVarInfo:
    """Test the get_var_info function."""

    @staticmethod
    def _create_test_dataset(data, var_name, long_name, units):
        """Helper to create a test dataset."""
        return xr.Dataset(
            {var_name: (("y", "x"), data)},
            attrs={"long_name": long_name, "units": units},
        )

    @pytest.mark.parametrize("dry_run", [False, True])
    def test_no_nan_in_data_suggests_delete(self, dry_run):
        """If data has no NaNs, default suggestion should be to delete the attribute."""
        data = np.array([[1.0, 2.0], [3.0, 4.0]])
        ds = self._create_test_dataset(data, TEST_VAR_NAME, "Test Variable", "test_units")
        _, config = get_var_info(
            TEST_VAR_NAME, ds, "/path/to/file", delete_if_none_filled=False, dry_run=dry_run
        )
        assert config.default_value == USER_REQ_DELETE
        assert config.allow_delete is True

    @pytest.mark.parametrize("dry_run", [False, True])
    def test_nan_in_data_prevents_delete(self, dry_run):
        """If data has NaNs, deleting the attribute should not be allowed."""
        data = np.array([[1.0, 2.0], [np.nan, 4.0]])
        ds = self._create_test_dataset(data, TEST_VAR_NAME, "Test Variable", "test_units")
        _, config = get_var_info(
            TEST_VAR_NAME, ds, "/path/to/file", delete_if_none_filled=False, dry_run=dry_run
        )
        assert config.default_value != USER_REQ_DELETE
        assert config.allow_delete is False

    @pytest.mark.parametrize("data_min", [10.0, 0.0, -1.0])
    def test_defaults_to_neg999_due_to_min(self, data_min):
        """For some minimum values, default should be -999."""
        data = np.array([[data_min, 20.0], [np.nan, 40.0]])
        ds = self._create_test_dataset(data, TEST_VAR_NAME, "Test Variable", "test_units")
        _, config = get_var_info(
            TEST_VAR_NAME, ds, "/path/to/file", delete_if_none_filled=False, dry_run=False
        )
        assert config.default_value == -999

    def test_special_varname_prefix_defaults_to_neg999(self):
        """Test if var name starts with special prefix like 'fertl_', default is -999."""
        var_name = VARSTARTS_TO_DEFAULT_NEG999[0] + "weuuewriebr"
        data = np.array([[-50.0, 20.0], [np.nan, 40.0]])
        ds = self._create_test_dataset(data, var_name, "Fertilizer", "g/m2")
        _, config = get_var_info(
            var_name, ds, "/path/to/file", delete_if_none_filled=False, dry_run=False
        )
        assert config.default_value == -999

    @pytest.mark.parametrize(
        "abs_path, expected",
        [
            ("/some/path/surfdata_map/file.nc", -999),
            ("/some/path/surfdata_map.nc", None),
        ],
    )
    def test_surfdata_map_defaults_to_neg999(self, abs_path, expected):
        """Test that a variable in a surfdata_map dir defaults to -999."""
        var_name = "grid1_to_grid2"
        data = np.array([[-50.0, 20.0], [np.nan, 40.0]])
        ds = self._create_test_dataset(data, var_name, "Mapping", "none")
        _, config = get_var_info(var_name, ds, abs_path, delete_if_none_filled=False, dry_run=False)
        assert config.default_value == expected

    @pytest.mark.parametrize("dry_run", [False, True])
    def test_no_special_case_no_default(self, dry_run):
        """If no special condition is met, there should be no default value."""
        data = np.array([[-50.0, -20.0], [np.nan, -10.0]])
        ds = self._create_test_dataset(data, TEST_VAR_NAME, "Test Variable", "test_units")
        _, config = get_var_info(
            TEST_VAR_NAME, ds, "/path/to/file", delete_if_none_filled=False, dry_run=dry_run
        )
        assert config.default_value is None

    def test_returns_correct_var_context(self):
        """Ensure VarContext is returned with correct info."""
        data = np.array([[1.0, 2.0]], dtype=np.float32)
        ds = self._create_test_dataset(data, TEST_VAR_NAME, "Test Variable", "test_units")
        file_path = "/path/to/file.nc"
        var_context, _ = get_var_info(
            TEST_VAR_NAME, ds, file_path, delete_if_none_filled=False, dry_run=False
        )
        assert var_context.var_name == TEST_VAR_NAME
        assert var_context.target_type == float
        assert var_context.file_path == file_path


class TestFindUserNlFiles:
    """Tests of find_user_nl_files()"""

    def test_find_user_nl_files(self, tmp_path):
        """Test find_user_nl_files()"""
        found_toplevel = f"{tmp_path}/user_nl_clm"
        Path(found_toplevel).touch()

        # Create files and dirs
        found_secondlevel = f"{tmp_path}/some_dir/user_nl_something"
        Path(found_secondlevel).parent.mkdir()
        Path(found_secondlevel).touch()
        notfound = f"{tmp_path}/some_dir/different_user_nl_confusing"
        Path(notfound).touch()

        # Get and check results
        results = find_user_nl_files(str(tmp_path))
        assert len(results) == 2
        assert found_toplevel in results
        assert found_secondlevel in results
        assert notfound not in results


class TestReplaceEnvVarsInNetcdfPath:
    """Tests of replace_env_vars_in_netcdf_paths()"""

    @pytest.mark.parametrize("nc_in", ["file_name.nc", "/abs/path/file.nc"])
    def test_replace_env_vars_in_netcdf_paths_nosubs(self, nc_in):
        """Test that replace_env_vars_in_netcdf_paths() works without any substitutions"""
        nc_out = replace_env_vars_in_netcdf_paths(nc_in)
        assert nc_out == nc_in

    @pytest.mark.parametrize("nc_in", ["$DIN_LOC_ROOT/file.nc", "${DIN_LOC_ROOT}/file.nc"])
    def test_replace_env_vars_in_netcdf_paths_dlr(self, nc_in):
        """Test that replace_env_vars_in_netcdf_paths() works when replacing DIN_LOC_ROOT"""
        nc_out = replace_env_vars_in_netcdf_paths(nc_in)

        expected = os.path.join(get_replacement_fill_values.INPUTDATA_PREFIX, "file.nc")

        assert nc_out == expected


class TestHowNetcdfIsReferencedInFile:
    """Tests of how_netcdf_is_referenced_in_file()."""

    def return_input(self, x):
        """Take one input argument and return it"""
        return x

    @pytest.fixture(autouse=True)
    def mock_convert_to_absolute_path(self, monkeypatch):
        """Mock convert_to_absolute_path() to just return what it was given"""
        mock = MagicMock(side_effect=lambda x, *args, **kwargs: x)
        monkeypatch.setattr(get_replacement_fill_values, "convert_to_absolute_path", mock)
        return mock

    @pytest.fixture(autouse=True)
    def mock_replace_env_vars_in_netcdf_paths(self, monkeypatch):
        """Mock replace_env_vars_in_netcdf_paths() to just return what it was given"""
        mock = MagicMock(side_effect=lambda x, *args, **kwargs: x)
        monkeypatch.setattr(get_replacement_fill_values, "replace_env_vars_in_netcdf_paths", mock)
        return mock

    def test_how_netcdf_is_referenced_in_file_1found_once(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """Test how_netcdf_is_referenced_in_file() for one netCDF file present once in text file"""
        nc_file = "file.nc"
        monkeypatch.setattr(
            get_replacement_fill_values,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: [nc_file],
        )
        set_of_how_this_netcdf_appears = (
            get_replacement_fill_values.how_netcdf_is_referenced_in_file("dummy", nc_file)
        )
        assert mock_convert_to_absolute_path.call_count == 2
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 1
        assert set_of_how_this_netcdf_appears == {nc_file}

    def test_how_netcdf_is_referenced_in_file_1found_twice_same(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """
        Test how_netcdf_is_referenced_in_file() for one netCDF file present twice in text file in
        the exact same way
        """
        nc_file = "file.nc"
        monkeypatch.setattr(
            get_replacement_fill_values,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: [nc_file, nc_file],
        )
        set_of_how_this_netcdf_appears = (
            get_replacement_fill_values.how_netcdf_is_referenced_in_file("dummy", nc_file)
        )
        assert mock_convert_to_absolute_path.call_count == 3
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 2
        assert set_of_how_this_netcdf_appears == {nc_file}

    def test_how_netcdf_is_referenced_in_file_1found_twice_diff(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """
        Test how_netcdf_is_referenced_in_file() for one netCDF file present twice in text file in
        different ways
        """
        nc_file = "file.nc"
        nc_file2 = "abc123" + nc_file
        monkeypatch.setattr(
            get_replacement_fill_values,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: [nc_file, nc_file2],
        )
        set_of_how_this_netcdf_appears = (
            get_replacement_fill_values.how_netcdf_is_referenced_in_file("dummy", nc_file)
        )
        assert mock_convert_to_absolute_path.call_count == 3
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 2
        assert set_of_how_this_netcdf_appears == {nc_file}

    def test_how_netcdf_is_referenced_in_file_2found(
        self, monkeypatch, mock_convert_to_absolute_path, mock_replace_env_vars_in_netcdf_paths
    ):
        """Test how_netcdf_is_referenced_in_file() for two netCDF files present in text file"""
        nc_file = "file.nc"
        nc_files = [nc_file, "file2.nc"]
        monkeypatch.setattr(
            get_replacement_fill_values,
            "extract_file_paths_from_file",
            lambda *args, **kwargs: nc_files,
        )
        set_of_how_this_netcdf_appears = (
            get_replacement_fill_values.how_netcdf_is_referenced_in_file("dummy", nc_file)
        )
        assert mock_convert_to_absolute_path.call_count == 3
        assert mock_replace_env_vars_in_netcdf_paths.call_count == 2
        assert set_of_how_this_netcdf_appears == {nc_file}



