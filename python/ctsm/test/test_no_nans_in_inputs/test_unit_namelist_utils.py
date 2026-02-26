"""Unit tests of functions in namelist_utils (anything not touching filesystem)"""

import os

import pytest

from ctsm.no_nans_in_inputs.namelist_utils import _replace_env_vars_in_netcdf_paths


class TestReplaceEnvVarsInNetcdfPath:
    """Tests of _replace_env_vars_in_netcdf_paths()"""

    @pytest.mark.parametrize("nc_in", ["file_name.nc", "/abs/path/file.nc"])
    def test_replace_env_vars_in_netcdf_paths_nosubs(self, nc_in):
        """Test that _replace_env_vars_in_netcdf_paths() works without any substitutions"""
        nc_out = _replace_env_vars_in_netcdf_paths(nc_in)
        assert nc_out == nc_in

    @pytest.mark.parametrize("nc_in", ["$DIN_LOC_ROOT/file.nc", "${DIN_LOC_ROOT}/file.nc"])
    def test_replace_env_vars_in_netcdf_paths_dlr(self, tmp_path, nc_in):
        """Test that _replace_env_vars_in_netcdf_paths() works when replacing DIN_LOC_ROOT"""
        nc_out = _replace_env_vars_in_netcdf_paths(nc_in)

        # Because we've patched INPUTDATA_PREFIX to be tmp_path
        expected = os.path.join(tmp_path, "file.nc")

        assert nc_out == expected
