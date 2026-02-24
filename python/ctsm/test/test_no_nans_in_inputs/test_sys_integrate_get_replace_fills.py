"""
Tests of the integrated get_replacement_fill_values.py -> replace_fill_values.py pipeline
"""

import os
from unittest.mock import patch
import xml.etree.ElementTree as ET

import pytest
import numpy as np
import xarray as xr

from ctsm.no_nans_in_inputs.constants import (
    ATTR,
    NEW_FILLVALUES_FILE,
    OPEN_DS_KWARGS,
    USER_REQ_DELETE,
)
from ctsm.no_nans_in_inputs import get_replacement_fill_values
from ctsm.no_nans_in_inputs.replace_fill_values import main as replace_fill_values
from ctsm.no_nans_in_inputs.replace_fill_values import get_output_filename

# Test constants
TEST_VAR_TEMP = "temp"
TEST_VAR_PRESSURE = "pressure"
TEST_OUTPUT_FILE = "output.nc"
TEST_FILL_VALUE = -123.4


@pytest.fixture(name="test_netcdf_file")
def fixture_test_netcdf_file(tmp_path):
    """Create a temporary NetCDF file for testing."""
    test_file = tmp_path / "lnd" / "clm2" / "test.nc"
    os.makedirs(os.path.dirname(str(test_file)))

    # Create a simple NetCDF file with float variables that have NaN fill values
    # (NetCDF doesn't allow NaN for integer types, and our scripts only work on
    # variables that already have NaN fill values)
    ds = xr.Dataset(
        {
            TEST_VAR_TEMP: xr.DataArray(
                np.array([np.nan, 2.0, 3.0], dtype=np.float32),
                dims=["time"],
                attrs={ATTR: np.float32(np.nan)},
            ),
            TEST_VAR_PRESSURE: xr.DataArray(
                np.array([1000.0, 1010.0, 1020.0], dtype=np.float64),
                dims=["time"],
                attrs={ATTR: np.float64(np.nan)},
            ),
        }
    )
    ds.to_netcdf(str(test_file))
    ds.close()

    yield str(test_file)


@pytest.mark.parametrize("abs_or_rel", ["abs", "rel"])
def test_integrate_get_replace(tmp_path, test_netcdf_file, create_mock_xml_file, abs_or_rel):
    """Test the integrated get -> replace pipeline"""

    # Get the path to put in the XML
    if abs_or_rel == "abs":
        netcdf_path_for_xml = test_netcdf_file
    elif abs_or_rel == "rel":
        netcdf_path_for_xml = os.path.relpath(test_netcdf_file, start=tmp_path)
    else:
        raise RuntimeError(f"Unrecognized {abs_or_rel=}")

    # Write the XML file
    xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{netcdf_path_for_xml}</paramfile>
</namelist_defaults>
"""
    xml_file = create_mock_xml_file(xml_content)

    # Call get_replacement_fill_values.py
    with patch("sys.argv", ["get_replacement_fill_values.py"]):
        with patch(
            "builtins.input",
            side_effect=[USER_REQ_DELETE, str(TEST_FILL_VALUE)],  # alphabetically 1st, then 2nd var
        ):
            get_replacement_fill_values.main()
            print(f"{NEW_FILLVALUES_FILE=}")

    # Call replace_fill_values.py
    with patch("sys.argv", ["replace_fill_values.py"]):
        replace_fill_values()

    # Check the output file
    output_file = get_output_filename(str(test_netcdf_file))
    assert os.path.exists(output_file)
    ds = xr.open_dataset(output_file, **OPEN_DS_KWARGS)
    assert ATTR in ds["temp"].encoding
    assert ds["temp"].encoding[ATTR] == TEST_FILL_VALUE
    assert np.isnan(ds["temp"].values[0])
    assert ATTR not in ds["pressure"].encoding

    # Check that the XML points to the output file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    paramfile = root.find("paramfile")
    assert paramfile is not None
    assert paramfile.text == get_output_filename(netcdf_path_for_xml)
