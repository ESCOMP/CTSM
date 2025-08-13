#!/usr/bin/env python3

"""System tests for set_paramfile"""

import unittest
import os
import sys
import shutil
import tempfile
import numpy as np
import xarray as xr
from netCDF4 import Dataset  # pylint: disable=no-name-in-module

from ctsm import unit_testing

from ctsm.netcdf_utils import get_netcdf_format
from ctsm.param_utils import set_paramfile as sp
from ctsm.param_utils.paramfile_shared import open_paramfile

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


class TestSysSetParamfile(unittest.TestCase):
    """System tests of set_paramfile"""

    def setUp(self):
        self.orig_argv = sys.argv
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        sys.argv = self.orig_argv
        shutil.rmtree(self.tempdir, ignore_errors=True)

    def test_set_paramfile_copyfile(self):
        """Test that set_paramfile can straight-up copy to a new file"""
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = ["set_paramfile", "-i", PARAMFILE, "-o", output_path]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that contents are functionally identical
        self.assertEqual(ds_in, ds_out)

        # Check that both are the same kind of netCDF
        self.assertEqual(get_netcdf_format(PARAMFILE), get_netcdf_format(output_path))

    def test_set_paramfile_extractvars(self):
        """Test that set_paramfile can copy to a new file with only some requested variables"""
        output_path = os.path.join(self.tempdir, "output.nc")
        vars_to_include = ["a_coef", "bgc_cn_s2"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-v",
            ",".join(vars_to_include),
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that output file variables/coords are what we expect
        # (use set() so order doesn't matter)
        expected_var_list = set(vars_to_include + list(ds_in.coords))
        self.assertEqual(set(ds_out.variables), expected_var_list)

        # Check that included variables/coords match
        for var in vars_to_include:
            self.assertEqual(ds_in[var], ds_out[var])

        # Check that both are the same kind of netCDF
        self.assertEqual(get_netcdf_format(PARAMFILE), get_netcdf_format(output_path))

    def test_set_paramfile_extractpfts(self):
        """Test that set_paramfile can copy to a new file with only some requested PFTs"""
        output_path = os.path.join(self.tempdir, "output.nc")
        pfts_to_include = ["not_vegetated", "needleleaf_evergreen_temperate_tree"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            ",".join(pfts_to_include),
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that included variables/coords match
        for var in ds_in.variables:
            if sp.PFTNAME_VAR in ds_in[var].coords:
                self.assertTrue(ds_in[var].isel(pft=[0, 1]).equals(ds_out[var]))
            else:
                self.assertTrue(ds_in[var].equals(ds_out[var]))

    def test_set_paramfile_changeparams_scalar_errors_given_list(self):
        """Test that set_paramfile errors if given a list for a scalar parameter"""
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "a_coef=0.87,0.91",
        ]
        with self.assertRaisesRegex(RuntimeError, "Incorrect N dims"):
            sp.main()

    def test_set_paramfile_changeparam_1d_errors_given_scalar(self):
        """Test that set_paramfile errors if given a scalar for a 1-d parameter"""
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "xl=0.724",
        ]
        with self.assertRaisesRegex(RuntimeError, "Incorrect N dims"):
            sp.main()

    def test_set_paramfile_changeparams_scalar_double(self):
        """Test that set_paramfile can copy to a new file with some scalar double params changed"""
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "a_coef=0.87",
            "bgc_cn_s2=87",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        for var in ds_in.variables:
            # Check that all variables/coords are equal except the ones we changed, which should be
            # set to what we asked
            if var == "a_coef":
                self.assertTrue(ds_in[var].values == 0.13)
                self.assertTrue(ds_out[var].values == 0.87)
            elif var == "bgc_cn_s2":
                self.assertTrue(ds_in[var].values == 11)
                self.assertTrue(ds_out[var].values == 87)
            else:
                self.assertTrue(ds_in[var].equals(ds_out[var]))

            # Check that data type hasn't changed
            self.assertTrue(ds_in[var].dtype == ds_out[var].dtype)

            # Check that fill value hasn't changed
            if "_FillValue" in ds_in[var].encoding:
                fv_in = ds_in[var].encoding["_FillValue"]
                fv_out = ds_out[var].encoding["_FillValue"]
                if isinstance(fv_in, bytes):
                    self.assertTrue(isinstance(fv_out, bytes))
                    self.assertEqual(fv_in, fv_out)
                else:
                    self.assertEqual(np.isnan(fv_in), np.isnan(fv_out))
                    if not np.isnan(fv_in):
                        self.assertEqual(fv_in, fv_out)

    def test_set_paramfile_changeparams_scalar_int(self):
        """Test that set_paramfile can copy to a new file with a scalar integer param changed"""
        output_path = os.path.join(self.tempdir, "output.nc")
        this_var = "upplim_destruct_metamorph"
        new_value = 1987
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            f"upplim_destruct_metamorph={new_value}",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that the variable in question is actually an integer to begin with
        self.assertTrue(sp.is_integer(ds_in[this_var].values))
        # Also check that it actually differs from our new value
        self.assertTrue(ds_in[this_var].values != new_value)

        for var in ds_in.variables:
            # Check that all variables/coords are equal except the one we changed, which should be
            # set to what we asked
            if var == this_var:
                self.assertTrue(ds_out[var].values == new_value)
            else:
                self.assertTrue(ds_in[var].equals(ds_out[var]))

            # Check that data type hasn't changed
            self.assertTrue(ds_in[var].dtype == ds_out[var].dtype)

    def test_set_paramfile_extractpfts_changeparam_dbl(self):
        """
        Test that set_paramfile can (1) copy to a new file with only some requested PFTs and (2)
        change the values of double parameters of those PFTs
        """
        output_path = os.path.join(self.tempdir, "output.nc")
        pfts_to_include = ["not_vegetated", "needleleaf_evergreen_temperate_tree"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            ",".join(pfts_to_include),
            "xl=0.724,0.87",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that included variables/coords match as expected
        for var in ds_in.variables:
            if var == "xl":
                self.assertTrue(np.array_equal(np.array([0.724, 0.87]), ds_out[var].values))
            elif sp.PFTNAME_VAR in ds_in[var].coords:
                self.assertTrue(ds_in[var].isel(pft=[0, 1]).equals(ds_out[var]))
            else:
                self.assertTrue(ds_in[var].equals(ds_out[var]))

    def test_set_paramfile_extractpfts_changeparam_int(self):
        """
        Test that set_paramfile can (1) copy to a new file with only some requested PFTs and (2)
        change the values of integer parameters of those PFTs
        """
        output_path = os.path.join(self.tempdir, "output.nc")
        pfts_to_include = ["not_vegetated", "needleleaf_evergreen_temperate_tree"]
        this_var = "max_NH_planting_date"
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            ",".join(pfts_to_include),
            f"{this_var}=1986,1987",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that the variable in question is actually an integer to begin with
        self.assertTrue(sp.is_integer(ds_in[this_var].values))
        # Also check that it actually differs from our new values
        self.assertTrue(ds_in[this_var].values[0] != 1986)
        self.assertTrue(ds_in[this_var].values[1] != 1987)

        # Check that included variables/coords match as expected
        for var in ds_in.variables:
            if var == this_var:
                self.assertTrue(np.array_equal(np.array([1986, 1987]), ds_out[var].values))
            elif sp.PFTNAME_VAR in ds_in[var].coords:
                self.assertTrue(ds_in[var].isel(pft=[0, 1]).equals(ds_out[var]))
            else:
                self.assertTrue(ds_in[var].equals(ds_out[var]))

    def test_set_paramfile_setparams_scalar_double_tonan_with_nan(self):
        """Test setting scalar double to NaN using 'nan'"""
        output_path = os.path.join(self.tempdir, "output.nc")
        this_var = "a_coef"
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            f"{this_var}=nan",
        ]

        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_out = open_paramfile(output_path, mask_and_scale=True)
        self.assertTrue(np.isnan(ds_out[this_var]))

    def test_set_paramfile_setparams_scalar_double_tonan_with_nancaps(self):
        """Test setting scalar double to NaN using 'NaN'"""
        output_path = os.path.join(self.tempdir, "output.nc")
        this_var = "a_coef"
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            f"{this_var}=NaN",
        ]

        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_out = open_paramfile(output_path, mask_and_scale=True)
        self.assertTrue(np.isnan(ds_out[this_var]))

    def test_set_paramfile_setparams_pft_double_tonan_with_nan(self):
        """Test setting PFT-dimensioned double to NaN using 'nan'"""
        output_path = os.path.join(self.tempdir, "output.nc")
        this_var = "xl"
        pfts_to_include = ["not_vegetated", "needleleaf_evergreen_temperate_tree"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            ",".join(pfts_to_include),
            f"{this_var}=nan,nan",
        ]

        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_out = open_paramfile(output_path, mask_and_scale=True)
        self.assertTrue(all(np.isnan(ds_out[this_var])))

    def test_set_paramfile_setparams_nan_but_no_fillvalue(self):
        """Test that NotImplementedError is given if trying to set NaN but param has no FillValue"""

        # Create paramfile with a double variable without fill value
        input_path = os.path.join(self.tempdir, "input.nc")
        ds = open_paramfile(PARAMFILE, mask_and_scale=True)
        new_param_name = "new_param_abc123"
        ds[new_param_name] = xr.DataArray(data=np.array(3.14))
        ds[new_param_name].encoding["_FillValue"] = None
        self.assertTrue(new_param_name in ds)
        ds.to_netcdf(input_path)

        # Check that it doesn't have fill value
        ds_in = open_paramfile(input_path, mask_and_scale=True)
        self.assertFalse("_FillValue" in ds_in[new_param_name].encoding)

        # Call set_paramfile, trying to set it to NaN
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = [
            "set_paramfile",
            "-i",
            input_path,
            "-o",
            output_path,
            f"{new_param_name}=nan",
        ]
        with self.assertRaisesRegex(
            NotImplementedError, "Not able to set NaN if parameter doesn't already have fill value:"
        ):
            sp.main()

    def test_set_paramfile_setparams_scalar_int_tonan_with_nan(self):
        """Test that NotImplementedError is given if trying to set NaN for an integer"""
        output_path = os.path.join(self.tempdir, "output.nc")
        this_var = "upplim_destruct_metamorph"
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            f"{this_var}=nan",
        ]

        with self.assertRaisesRegex(
            NotImplementedError, "Not able to set NaN for integer parameters:"
        ):
            sp.main()

    # TODO: Add test changing vector int to NaN using nan
    # TODO: Add test changing scalar double to NaN using _
    # TODO: Add test changing scalar int to NaN using _
    # TODO: Add test changing vector double to NaN using _
    # TODO: Add test changing vector int to NaN using _

    # TODO: Add test overwriting existing scalar double NaN using nan (lowercase)

    # TODO: Test changing param when extracting just one PFT

    # TODO: Test changing PFT name

    def test_set_paramfile_changeparam_multidim_errors(self):
        """
        Test that set_paramfile errors if requesting change of a multi-dimensional parameter. This
        test will obviously need to be replaced once that functionality is added.
        """
        output_path = os.path.join(self.tempdir, "output.nc")
        pfts_to_include = ["not_vegetated", "needleleaf_evergreen_temperate_tree"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "mimics_till_decompk_multipliers=dummy",
        ]

        with self.assertRaises(NotImplementedError):
            sp.main()


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
