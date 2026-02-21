#!/usr/bin/env python3
"""
System tests for get_replacement_fill_values.py script.

Tests functions that require file I/O.
"""

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import ATTR
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    extract_file_paths_from_xml,
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

    def test_true_for_nan_fill_float32(self, tmp_path):
        """Test that a float32 variable with NaN fill value returns True."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                "temp": xr.DataArray(
                    np.array([1.0, 2.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(np.nan)},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert var_has_nan_fill(ds_read, "temp")
        ds_read.close()

    def test_true_for_nan_fill_float64(self, tmp_path):
        """Test that a float64 variable with NaN fill value returns True."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                "pressure": xr.DataArray(
                    np.array([1000.0, 1010.0], dtype=np.float64),
                    dims=["time"],
                    attrs={ATTR: np.float64(np.nan)},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert var_has_nan_fill(ds_read, "pressure")
        ds_read.close()

    def test_false_for_numeric_fill(self, tmp_path):
        """Test that a variable with a numeric (non-NaN) fill value returns False."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                "temp": xr.DataArray(
                    np.array([1.0, 2.0], dtype=np.float32),
                    dims=["time"],
                    attrs={ATTR: np.float32(-999.0)},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert not var_has_nan_fill(ds_read, "temp")
        ds_read.close()

    def test_false_for_no_fill_attr(self, tmp_path):
        """Test that a variable without a _FillValue attribute returns False."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                "temp": xr.DataArray(
                    np.array([1.0, 2.0], dtype=np.float32),
                    dims=["time"],
                ),
            }
        )
        # Use encoding to prevent xarray from adding a default _FillValue
        ds.to_netcdf(str(test_file), encoding={"temp": {ATTR: None}})
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert not var_has_nan_fill(ds_read, "temp")
        ds_read.close()

    def test_false_for_integer_fill(self, tmp_path):
        """Test that an integer variable with a fill value returns False (can't be NaN)."""
        test_file = tmp_path / "test.nc"
        ds = xr.Dataset(
            {
                "count": xr.DataArray(
                    np.array([1, 2, 3], dtype=np.int32),
                    dims=["time"],
                    attrs={ATTR: np.int32(-999)},
                ),
            }
        )
        ds.to_netcdf(str(test_file))
        ds_read = xr.open_dataset(
            str(test_file), decode_cf=False, decode_timedelta=False, decode_times=False
        )
        assert not var_has_nan_fill(ds_read, "count")
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
