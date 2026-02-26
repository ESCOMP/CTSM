"""System tests of netcdf_utils (anything touching filesystem)"""

# pylint: disable=too-few-public-methods

import pytest
import numpy as np
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR, USER_REQ_DELETE, VARSTARTS_TO_DEFAULT_NEG999
from ctsm.no_nans_in_inputs.netcdf_utils import (
    get_var_info,
    get_vars_with_nan_fills,
    show_ncdump_for_variable,
    var_data_has_nan,
)


TEST_VAR_NAME = "test_var"


class TestGetVarsWithNanFills:
    """Test the get_vars_with_nan_fills function."""

    @pytest.mark.parametrize(
        "fill_value, expected",
        [
            (np.float32(np.nan), True),
            (np.float64(np.nan), True),
            (np.float32(-999.0), False),
            (np.int32(-999), False),
            (None, False),
        ],
    )
    def test_fill_value_detection(self, tmp_path, fill_value, expected):
        """Test that get_vars_with_nan_fills correctly detects NaN/non-NaN/absent fill values."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        dtype = np.float32 if fill_value is None else type(fill_value)
        ds = xr.Dataset(
            {
                var_name: xr.DataArray(
                    np.array([1, 2], dtype=dtype),
                    dims=["time"],
                ),
            }
        )
        # Use encoding to set (or suppress) the _FillValue
        encoding = {var_name: {ATTR: fill_value}}
        ds.to_netcdf(str(test_file), encoding=encoding)
        vars_with_nan_fills = get_vars_with_nan_fills(str(test_file))
        assert (var_name in vars_with_nan_fills) == expected


class TestShowNcdumpForVariable:
    """Test the show_ncdump_for_variable function."""

    TEST_VAR = "temp"

    @pytest.fixture
    def test_netcdf_file(self, tmp_path):
        """Create a temporary NetCDF file for testing."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                self.TEST_VAR: xr.DataArray(
                    np.array([1.0, 2.0, 3.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan), "long_name": "temperature"},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds.close()
        return str(test_file)

    def test_prints_matching_lines(self, test_netcdf_file, capsys):
        """Test that matching lines from ncdump are printed."""
        show_ncdump_for_variable(test_netcdf_file, self.TEST_VAR)
        captured = capsys.readouterr()
        assert self.TEST_VAR in captured.out
        assert "Lines matching" in captured.out

    def test_prints_no_match_message(self, test_netcdf_file, capsys):
        """Test that a message is printed when no lines match."""
        show_ncdump_for_variable(test_netcdf_file, "nonexistent_var")
        captured = capsys.readouterr()
        assert "No lines found matching" in captured.out

    def test_none_file_path(self, capsys):
        """Test that None file path prints a message and returns."""
        show_ncdump_for_variable(None, self.TEST_VAR)
        captured = capsys.readouterr()
        assert "No file path available" in captured.out

    def test_nonexistent_file(self, capsys):
        """Test that a nonexistent file prints an error."""
        show_ncdump_for_variable("/nonexistent/file.nc", self.TEST_VAR)
        captured = capsys.readouterr()
        assert "Error running ncdump" in captured.out


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
