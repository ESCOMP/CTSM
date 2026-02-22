#!/usr/bin/env python3
"""
System tests for get_replacement_fill_values.py script.

Tests functions that require file I/O.
"""

import sys
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    extract_file_paths_from_xml,
    main as main_func,
    show_ncdump_for_variable,
    var_has_nan_fill,
)


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


class TestVarHasNanFill:
    """Test the var_has_nan_fill function."""

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
        """Test that var_has_nan_fill correctly detects NaN vs non-NaN vs absent fill values."""
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
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert var_has_nan_fill(ds_read, "temp") == expected
        ds_read.close()

    def test_checks_correct_variable(self, tmp_path):
        """Test that only the specified variable is checked."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                "has_nan": xr.DataArray(
                    np.array([1.0, 2.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
                "no_nan": xr.DataArray(
                    np.array([3.0, 4.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(-999.0)},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert var_has_nan_fill(ds_read, "has_nan")
        assert not var_has_nan_fill(ds_read, "no_nan")
        ds_read.close()


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
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.load_bad_files",
        return_value={"/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.xr.open_dataset")
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.var_has_nan_fill",
        return_value=True,
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
    def test_main_happy_path(
        self,
        mock_collect,
        mock_var_has_nan,
        mock_open_dataset,
        mock_exists,
        mock_load_bad,
        mock_extract,
        mock_check_write,
    ):  # pylint: disable=unused-argument
        """Test main function with a single matching file."""
        # Setup mock dataset
        mock_ds = MagicMock(spec=xr.Dataset)
        # Make the mock dataset iterable
        mock_ds.__iter__.return_value = iter(["temp"])

        # Setup context manager mock for open_dataset
        mock_open_dataset.return_value = mock_ds

        result = main_func()

        assert result == 0
        mock_extract.assert_called_once()
        mock_load_bad.assert_called_once()
        mock_open_dataset.assert_called_once_with(
            "/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/test.nc",
            decode_cf=False,
            decode_timedelta=False,
            decode_times=False,
        )
        mock_var_has_nan.assert_called_with(mock_ds, "temp")
        mock_collect.assert_called_once()

        # Check the arguments passed to collect_new_fill_values
        args, kwargs = mock_collect.call_args
        matches = args[0]
        assert len(matches) == 1
        assert matches[0] == (
            "lnd/clm2/test.nc",
            "/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/test.nc",
        )
        assert "delete_if_none_filled" in kwargs
        assert not kwargs["delete_if_none_filled"]

    @patch("sys.argv", ["get_replacement_fill_values.py", "--delete-if-none-filled"])
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.check_write_access",
        return_value=True,
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.extract_file_paths_from_xml",
        return_value={"lnd/clm2/test.nc"},
    )
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.load_bad_files",
        return_value={"/glade/campaign/cesm/cesmdata/cseg/inputdata/lnd/clm2/test.nc"},
    )
    @patch("os.path.exists", return_value=True)
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.xr.open_dataset")
    @patch(
        "ctsm.no_nans_in_inputs.get_replacement_fill_values.var_has_nan_fill",
        return_value=True,
    )
    @patch("ctsm.no_nans_in_inputs.get_replacement_fill_values.collect_new_fill_values")
    def test_main_with_delete_flag(
        self,
        mock_collect,
        mock_var_has_nan,
        mock_open_dataset,
        mock_exists,
        mock_load_bad,
        mock_extract,
        mock_check_write,
    ):  # pylint: disable=unused-argument
        """Test main function with --delete-if-none-filled flag."""
        mock_ds = MagicMock(spec=xr.Dataset)
        mock_ds.__iter__.return_value = iter(["temp"])
        mock_open_dataset.return_value = mock_ds

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
