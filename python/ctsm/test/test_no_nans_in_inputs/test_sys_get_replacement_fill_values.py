#!/usr/bin/env python3
"""
System tests for get_replacement_fill_values.py script.

Tests functions that require file I/O.
"""
# pylint: disable=too-few-public-methods

import os
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import (
    ATTR,
    OPEN_DS_KWARGS,
    USER_REQ_DELETE,
    USER_REQ_SKIP_FILE,
    USER_REQ_QUIT,
)
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    collect_new_fill_values,
    extract_file_paths_from_xml,
    file_has_nan_fill,
    get_vars_with_nan_fills,
    load_progress,
    main as main_func,
    create_empty_progress_dict_onefile,
    save_progress,
    show_ncdump_for_variable,
)
from ctsm.no_nans_in_inputs import get_replacement_fill_values


# Test constants
TEST_PATH_PARAM = "lnd/clm2/paramdata/test_params.nc"
TEST_PATH_SURF = "lnd/clm2/surfdata/test_surf.nc"
TEST_PATH_INIT = "lnd/clm2/initdata/test_init.nc"
TEST_PATH_OTHER = "share/meshes/test_mesh.nc"


class TestExtractFilePathsFromXml:
    """Test the extract_file_paths_from_xml function."""

    def test_extracts_multiple_paths(self, create_mock_xml_file):
        """Test extracting multiple file paths from XML."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <surfdata>{TEST_PATH_SURF}</surfdata>
    <initdata>{TEST_PATH_INIT}</initdata>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM, TEST_PATH_SURF, TEST_PATH_INIT}

    def test_ignores_non_lnd_clm2_paths(self, create_mock_xml_file):
        """Test that paths not containing lnd/clm2/ are ignored."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <meshfile>{TEST_PATH_OTHER}</meshfile>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM}
        assert TEST_PATH_OTHER not in result

    def test_ignores_non_path_text(self, create_mock_xml_file):
        """Test that non-path text content is ignored."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <some_setting>42</some_setting>
    <another_setting>.true.</another_setting>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM}

    def test_deduplicates_paths(self, create_mock_xml_file):
        """Test that duplicate paths are returned only once."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="clm6_0">{TEST_PATH_PARAM}</paramfile>
    <paramfile phys="clm5_0">{TEST_PATH_PARAM}</paramfile>
    <something_else>{TEST_PATH_PARAM}</something_else>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM}

    def test_handles_whitespace_around_path(self, create_mock_xml_file):
        """Test that leading/trailing whitespace around paths is handled."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>  {TEST_PATH_PARAM}  </paramfile>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM}

    def test_handles_multiline_path(self, create_mock_xml_file):
        """Test extracting a path that spans multiple lines (with leading newline)."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>
{TEST_PATH_PARAM}
</paramfile>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM}

    def test_empty_xml(self, create_mock_xml_file):
        """Test with XML that has no file paths."""
        xml_path = create_mock_xml_file(
            """<?xml version="1.0"?>
<namelist_defaults>
    <some_setting>42</some_setting>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == set()

    def test_file_not_found_exits(self, capsys):
        """Test that a missing XML file causes SystemExit."""
        nonexistent_path = "/nonexistent/file.xml"
        with pytest.raises(SystemExit):
            extract_file_paths_from_xml(nonexistent_path)
        captured = capsys.readouterr()
        assert nonexistent_path in captured.err
        assert "not found" in captured.err

    def test_handles_elements_with_attributes(self, create_mock_xml_file):
        """Test extracting paths from elements that have attributes."""
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="clm6_0">{TEST_PATH_PARAM}</paramfile>
    <paramfile phys="clm5_0">{TEST_PATH_SURF}</paramfile>
</namelist_defaults>
"""
        )
        result = extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM, TEST_PATH_SURF}


class TestFileHasNanFill:
    """Test the file_has_nan_fill function."""

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
        """Test that file_has_nan_fill correctly detects NaN vs non-NaN vs absent fill values."""
        test_file = tmp_path / "test.nc"
        dtype = np.float32 if fill_value is None else type(fill_value)
        ds = xr.Dataset(
            {
                "temp": xr.DataArray(
                    np.array([1, 2], dtype=dtype),
                    dims=["time"],
                ),
            }
        )
        # Use encoding to set (or suppress) the _FillValue
        encoding = {"temp": {ATTR: fill_value}}
        ds.to_netcdf(str(test_file), encoding=encoding)
        any_nan_fill, vars_with_nan_fills = file_has_nan_fill(str(test_file))
        assert any_nan_fill == expected
        assert len(vars_with_nan_fills) == int(expected)


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


class TestMain:
    """Test the main function of get_replacement_fill_values.py."""

    @patch("sys.argv", ["get_replacement_fill_values.py"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.file_has_nan_fill",
        return_value=(True, ["temp"]),
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
    def test_main_happy_path(
        self,
        mock_collect,
        mock_file_has_nan,
        mock_exists,
        mock_extract,
        mock_check_write,
        tmp_path,
        mock_xml_file_path,
    ):  # pylint: disable=unused-argument
        """Test main function with a single matching file."""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()

        assert result == 0
        mock_extract.assert_called_once()
        path_to_test_file_rel = "lnd/clm2/test.nc"
        path_to_test_file_abs = os.path.join(str(tmp_path), path_to_test_file_rel)
        mock_collect.assert_called_once()
        mock_check_write.assert_called_once()

        # Check the arguments passed to collect_new_fill_values
        args, kwargs = mock_collect.call_args
        progress = args[0]
        assert len(progress) == 1
        expected_dict = create_empty_progress_dict_onefile()
        expected_dict["found_in_files"] = {mock_xml_file_path: path_to_test_file_rel}
        expected_dict["vars_with_nan_fills"] = ["temp"]
        assert progress[path_to_test_file_abs] == expected_dict
        assert "delete_if_none_filled" in kwargs
        assert not kwargs["delete_if_none_filled"]
        assert "dry_run" in kwargs
        assert not kwargs["dry_run"]

    @patch("sys.argv", ["get_replacement_fill_values.py", "--dry-run"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.get_vars_with_nan_fills",
        return_value=["temp"],
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
    def test_main_dry_run(
        self,
        mock_collect,
        mock_get_vars_with_nan,
        mock_exists,
        mock_extract,
        mock_check_write,
        capsys,
        tmp_path,
        mock_xml_file_path,
    ):  # pylint: disable=unused-argument
        """Test main function with a single matching file under --dry-run"""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()

        assert result == 0
        mock_extract.assert_called_once()
        path_to_test_file_rel = "lnd/clm2/test.nc"
        path_to_test_file_abs = os.path.join(str(tmp_path), path_to_test_file_rel)
        mock_get_vars_with_nan.assert_called_with(path_to_test_file_abs)
        mock_collect.assert_called_once()

        # Shouldn't check write access in dry run
        mock_check_write.assert_not_called()

        # Check the arguments passed to collect_new_fill_values
        args, kwargs = mock_collect.call_args
        progress = args[0]
        expected_dict = create_empty_progress_dict_onefile()
        expected_dict["found_in_files"] = {mock_xml_file_path: path_to_test_file_rel}
        expected_dict["vars_with_nan_fills"] = ["temp"]
        assert len(progress) == 1
        assert progress[path_to_test_file_abs] == expected_dict
        assert "delete_if_none_filled" in kwargs
        assert not kwargs["delete_if_none_filled"]
        assert "dry_run" in kwargs
        assert kwargs["dry_run"]

        # Check stdout
        stdout = capsys.readouterr().out
        assert "Checking write access" not in stdout

    @patch("sys.argv", ["get_replacement_fill_values.py"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.get_vars_with_nan_fills",
        return_value=["temp"],
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
    def test_main_missing_file(
        self,
        mock_collect,
        mock_get_vars_with_nan,
        mock_extract,
        mock_check_write,
        capsys,
        tmp_path,
        mock_xml_file_path,
    ):  # pylint: disable=unused-argument
        """Test main function with a file that doesn't exist"""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()

        assert result == 0
        mock_extract.assert_called_once()
        path_to_test_file_rel = "lnd/clm2/test.nc"
        path_to_test_file_abs = os.path.join(str(tmp_path), path_to_test_file_rel)
        mock_get_vars_with_nan.assert_not_called()

        # Check the arguments passed to collect_new_fill_values
        mock_collect.assert_called_once()
        args, _ = mock_collect.call_args
        progress = args[0]
        assert progress == {}

        # Check stdout
        stdout = capsys.readouterr().out
        assert "1\tFiles not found" in stdout
        assert f"Not found: '{path_to_test_file_abs}'" in stdout

    @patch("sys.argv", ["get_replacement_fill_values.py", "--delete-if-none-filled"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.file_has_nan_fill",
        return_value=(True, ["temp"]),
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
    def test_main_with_delete_flag(
        self,
        mock_collect,
        mock_file_has_nan,
        mock_exists,
        mock_extract,
        mock_check_write,
    ):  # pylint: disable=unused-argument
        """Test main function with --delete-if-none-filled flag."""
        mock_ds = MagicMock(spec=xr.Dataset)
        mock_ds.__iter__.return_value = iter(["temp"])

        result = main_func()
        assert result == 0

        mock_collect.assert_called_once()
        _, kwargs = mock_collect.call_args
        assert "delete_if_none_filled" in kwargs
        assert kwargs["delete_if_none_filled"]

    @patch("sys.argv", ["get_replacement_fill_values.py"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=False,
    )
    def test_main_no_write_access_exits(
        self, mock_check_write, capsys
    ):  # pylint: disable=unused-argument
        """Test that main exits if there's no write access."""
        with pytest.raises(SystemExit) as e:
            main_func()
        assert e.value.code == 1
        captured = capsys.readouterr()
        assert "No write access" in captured.err


class TestProgressFunctions:
    """Test the save_progress and load_progress functions."""

    def test_save_and_load_progress(self, tmp_path):
        """Test that data is saved and loaded correctly."""
        progress_file = tmp_path / "progress.json"
        data_to_save = {
            "/path/to/file1.nc": {"var1": -999, "var2": USER_REQ_DELETE},
            "/path/to/file2.nc": {"var3": -999.0},
        }

        # Save data
        save_progress(data_to_save, str(progress_file))

        # Load data
        loaded_data = load_progress(str(progress_file))

        # Check that loaded data is identical to saved data
        assert loaded_data == data_to_save

    def test_load_nonexistent_file(self):
        """Test loading a nonexistent progress file returns an empty dict."""
        loaded_data = load_progress("/nonexistent/path/progress.json")
        assert loaded_data == {}

    def test_load_corrupted_json(self, tmp_path, capsys):
        """Test loading a corrupted JSON file returns an empty dict and warns."""
        progress_file = tmp_path / "progress.json"
        with open(progress_file, "w", encoding="utf-8") as f:
            f.write("this is not valid json")

        loaded_data = load_progress(str(progress_file))

        assert loaded_data == {}
        captured = capsys.readouterr()
        assert "Warning: Could not load progress file" in captured.err

    def test_save_io_error(self, capsys):
        """Test that an IOError during save is caught and a warning is printed."""
        with patch("builtins.open", side_effect=IOError("Permission denied")):
            save_progress({"some": "data"}, "/readonly/path/progress.json")
            captured = capsys.readouterr()
            assert "Warning: Could not save progress" in captured.err


class TestCollectNewFillValues:
    """Test the collect_new_fill_values function using real files."""

    def _create_test_netcdf(self, path, var_dict):
        """Helper function to create a NetCDF file for testing."""
        ds = xr.Dataset()
        for var_name, properties in var_dict.items():
            ds[var_name] = xr.DataArray(
                properties["data"],
                dims=["time"],
                attrs=properties["attrs"],
            )
        ds.to_netcdf(path)
        ds.close()

    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.save_progress")
    @patch("builtins.input", side_effect=["-999"])
    def test_user_provides_number(self, mock_input, mock_save, tmp_path, capsys):
        """Test happy path where user provides a numeric fill value."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {
                var_name: {
                    "data": np.array([1.0, 2.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan), "long_name": "temperature"},
                }
            },
        )
        expected_fill_value = -999.0
        expected_result = {str(test_file): create_empty_progress_dict_onefile()}
        assert not expected_result[str(test_file)]["new_fill_values"]

        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        result = collect_new_fill_values(progress)

        # Check final result
        expected_result[str(test_file)]["new_fill_values"] = {var_name: expected_fill_value}
        expected_result[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert result == expected_result

        # Check what was saved
        mock_save.assert_called_once()
        saved_data = mock_save.call_args[0][0]
        saved_fill_value = saved_data[str(test_file)]["new_fill_values"][var_name]
        assert saved_fill_value == expected_fill_value
        assert isinstance(saved_fill_value, float)

    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.save_progress")
    @patch("builtins.input", side_effect=[USER_REQ_SKIP_FILE])
    def test_user_skips_file(self, mock_input, mock_save, tmp_path):
        """Test that requesting skipfile skips the current file."""
        test_file = tmp_path / "test.nc"
        self._create_test_netcdf(
            test_file,
            {
                "temp": {
                    "data": np.array([1.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan)},
                },
                "precip": {
                    "data": np.array([2.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan)},
                },
            },
        )
        progress = {str(test_file): create_empty_progress_dict_onefile()}
        result = collect_new_fill_values(progress)

        # 'temp' is processed, user skips, neither it nor 'precip' are in results
        expected_result = {str(test_file): create_empty_progress_dict_onefile()}
        assert not expected_result[str(test_file)]["new_fill_values"]
        assert result == expected_result
        # Progress was never saved
        mock_save.assert_not_called()

    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.save_progress")
    @patch("builtins.input", side_effect=["-100", USER_REQ_SKIP_FILE])
    def test_user_enters_then_skips_file(self, mock_input, mock_save, tmp_path):
        """Test that requesting skipfile skips the current file but saves entered results."""
        test_file = tmp_path / "test.nc"
        data = {
            "temp": {
                "data": np.array([1.0], dtype=np.float32),
                "attrs": {ATTR: np.float32(np.nan)},
            },
            "precip": {
                "data": np.array([2.0], dtype=np.float32),
                "attrs": {ATTR: np.float32(np.nan)},
            },
        }
        self._create_test_netcdf(
            test_file,
            data,
        )
        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = data.keys()
        result = collect_new_fill_values(progress)

        # 'temp' is processed, user enters so it saves, then precip is skipped
        expected_result = {str(test_file): create_empty_progress_dict_onefile()}
        expected_result[str(test_file)]["new_fill_values"] = {"temp": -100.0}
        expected_result[str(test_file)]["vars_with_nan_fills"] = data.keys()
        assert result == expected_result
        # Progress was saved only after temp
        mock_save.assert_called_once()

    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.save_progress")
    @patch("builtins.input", side_effect=[USER_REQ_QUIT])
    def test_user_quits(self, mock_input, mock_save, tmp_path, mock_progress_file):
        """Test that requesting quit exits the collection loop."""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {var_name: {"data": [1.0], "attrs": {ATTR: np.nan}}},
        )
        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert not progress[str(test_file)]["new_fill_values"]

        assert not os.path.exists(str(mock_progress_file))

        with pytest.raises(SystemExit) as e:
            collect_new_fill_values(progress)
        assert e.value.code == 0

        # No value was collected, so save_progress is not called
        mock_save.assert_not_called()

    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.save_progress")
    def test_dry_run(self, mock_save, tmp_path, capsys):
        """Test --dry-run"""
        test_file = tmp_path / "test.nc"
        var_name = "temp"
        self._create_test_netcdf(
            test_file,
            {
                var_name: {
                    "data": np.array([1.0, 2.0], dtype=np.float32),
                    "attrs": {ATTR: np.float32(np.nan), "long_name": "temperature"},
                }
            },
        )

        progress = {str(test_file): create_empty_progress_dict_onefile()}
        progress[str(test_file)]["vars_with_nan_fills"] = [var_name]
        result = collect_new_fill_values(progress, dry_run=True)

        # Check that nothing was saved
        mock_save.assert_not_called()

        # Check that result is what we expect
        expected = {str(test_file): create_empty_progress_dict_onefile()}
        expected[str(test_file)]["vars_with_nan_fills"] = [var_name]
        assert result == expected

        # Check that result does print some info but not everything
        stdout = capsys.readouterr().out
        assert "Variable: temp" in stdout
        assert "new fill value(s)" not in stdout
