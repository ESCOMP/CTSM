#!/usr/bin/env python3

"""System tests for set_paramfile"""

import unittest
import os
import sys
import shutil
import tempfile
import numpy as np
import xarray as xr

from ctsm import unit_testing

from ctsm.netcdf_utils import get_netcdf_format
from ctsm.param_utils import set_paramfile as sp
from ctsm.param_utils.paramfile_shared import open_paramfile, are_paramfile_dataarrays_identical
from ctsm.param_utils.paramfile_shared import check_pfts_in_paramfile, get_selected_pft_indices

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


PARAMFILE = os.path.join(
    os.path.dirname(__file__), "testinputs", "ctsm5.3.041.Nfix_params.v13.c250221_upplim250.nc"
)


class TestSysSetParamfile(unittest.TestCase):
    """System tests of set_paramfile"""

    # pylint: disable=too-many-public-methods

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

    def test_set_paramfile_copy_without_adding_fillvalue(self):
        """Test that set_paramfile can copy to a new file without adding _FillValue"""
        input_path = os.path.join(self.tempdir, "input.nc")
        output_path = os.path.join(self.tempdir, "output.nc")
        param_name = "param0"

        # Save a test paramfile without _FillValue
        da = xr.DataArray(data=np.float64([1, 2, 3]))
        da.encoding["_FillValue"] = None
        ds = xr.Dataset(data_vars={param_name: da})
        ds.to_netcdf(input_path, encoding={param_name: {"_FillValue": None}})
        ds_in = open_paramfile(input_path)
        self.assertFalse("_FillValue" in ds_in[param_name].encoding)
        self.assertFalse("_FillValue" in ds_in[param_name].attrs)

        # Use set_paramfile to copy to new file
        sys.argv = ["set_paramfile", "-i", input_path, "-o", output_path]
        sp.main()

        # Check that _FillValue wasn't added
        ds_out = open_paramfile(output_path)
        self.assertFalse("_FillValue" in ds_out[param_name].encoding)
        self.assertFalse("_FillValue" in ds_out[param_name].attrs)

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
            "--drop-other-pfts",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        # Check that included variables/coords match
        for var in ds_in.variables:
            actual = ds_out[var]
            if sp.PFTNAME_VAR in ds_in[var].coords:
                expected = ds_in[var].isel(pft=[0, 1])
            else:
                expected = ds_in[var]
            self.assertTrue(are_paramfile_dataarrays_identical(expected, actual))

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
                self.assertTrue(are_paramfile_dataarrays_identical(ds_in[var], ds_out[var]))

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

    def test_set_paramfile_changeparams_1d_double(self):
        """
        Test that set_paramfile can copy to a new file with a 1-d double param changed (not PFT-
        dimensioned)
        """
        output_path = os.path.join(self.tempdir, "output.nc")
        this_var = "mimics_fmet"
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            f"{this_var}=0.1,0.2,0.3,0.4",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_in = open_paramfile(PARAMFILE)
        ds_out = open_paramfile(output_path)

        for var in ds_in.variables:
            # Check that all variables/coords are equal except the ones we changed, which should be
            # set to what we asked
            if var == this_var:
                self.assertTrue(
                    np.array_equal(ds_in[var].values, np.array([0.75, 0.85, 0.013, 40]))
                )
                self.assertTrue(np.array_equal(ds_out[var].values, np.array([0.1, 0.2, 0.3, 0.4])))
            else:
                self.assertTrue(are_paramfile_dataarrays_identical(ds_in[var], ds_out[var]))

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
            "--drop-other-pfts",
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

    def test_set_paramfile_changeparam_dbl_onlysomepfts(self):
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

                # Changed values (first 2)
                this_slice = slice(0, 2)
                expected = np.array([0.724, 0.87])
                result = ds_out[var].isel(pft=this_slice).values
                self.assertTrue(np.array_equal(expected, result))

                # Preserved values (everything but the first 2)
                this_slice = slice(2, None)
                expected = ds_in["xl"].isel(pft=this_slice)
                result = ds_out["xl"].isel(pft=this_slice)
                self.assertTrue(are_paramfile_dataarrays_identical(expected, result))
            else:
                self.assertTrue(are_paramfile_dataarrays_identical(ds_in[var], ds_out[var]))

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
            "--drop-other-pfts",
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
                self.assertTrue(
                    are_paramfile_dataarrays_identical(ds_in[var].isel(pft=[0, 1]), ds_out[var])
                )
            else:
                self.assertTrue(are_paramfile_dataarrays_identical(ds_in[var], ds_out[var]))

    def test_set_paramfile_fill_value_scalar_double_nan(self):
        """
        Test that setting scalar double to fill value writes a literal NaN if that's the _FillValue
        """
        # Create paramfile with a double variable with fill value NaN
        input_path = os.path.join(self.tempdir, "input.nc")
        ds = open_paramfile(PARAMFILE, mask_and_scale=True)
        new_param_name = "new_param_abc123"
        ds[new_param_name] = xr.DataArray(data=np.array(3.14))
        ds[new_param_name].encoding["_FillValue"] = np.nan
        self.assertTrue(new_param_name in ds)
        ds.to_netcdf(input_path)

        # Check that its fill value is NaN
        ds_in = open_paramfile(input_path, mask_and_scale=True)
        self.assertTrue("_FillValue" in ds_in[new_param_name].encoding)
        self.assertTrue(np.isnan(ds_in[new_param_name].encoding["_FillValue"]))

        # Ask to set it to the FillValue
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = [
            "set_paramfile",
            "-i",
            input_path,
            "-o",
            output_path,
            f"{new_param_name}=nan",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))

        # Ensure it wrote a literal NaN
        ds_out = open_paramfile(output_path, mask_and_scale=False)
        self.assertTrue(np.isnan(ds_out[new_param_name]))

        # Ensure it preserved NaN FillValue
        ds_out = open_paramfile(output_path, mask_and_scale=True)
        self.assertTrue(np.isnan(ds_out[new_param_name].encoding["_FillValue"]))

    def test_set_paramfile_fill_value_scalar_double_real(self):
        """
        Test that setting scalar double to fill value does NOT write NaN if that's not the
        _FillValue
        """
        # Create paramfile with a double variable with fill value -999.9
        input_path = os.path.join(self.tempdir, "input.nc")
        ds = open_paramfile(PARAMFILE, mask_and_scale=True)
        new_param_name = "new_param_abc123"
        ds[new_param_name] = xr.DataArray(data=np.array(3.14))
        fill_value = -999.9
        ds[new_param_name].encoding["_FillValue"] = fill_value
        self.assertTrue(new_param_name in ds)
        ds.to_netcdf(input_path)

        # Check that its fill value is what we asked for
        ds_in = open_paramfile(input_path, mask_and_scale=True)
        self.assertTrue("_FillValue" in ds_in[new_param_name].encoding)
        self.assertEqual(fill_value, ds_in[new_param_name].encoding["_FillValue"])

        # Ask to set it to the FillValue
        output_path = os.path.join(self.tempdir, "output.nc")
        sys.argv = [
            "set_paramfile",
            "-i",
            input_path,
            "-o",
            output_path,
            f"{new_param_name}=nan",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))

        # Ensure it preserved FillValue
        ds_out = open_paramfile(output_path, mask_and_scale=True)
        self.assertEqual(fill_value, ds_out[new_param_name].encoding["_FillValue"])
        self.assertTrue(np.isnan(ds_out[new_param_name].values))

        # Ensure it wrote the FillValue and not a literal NaN
        ds_out = open_paramfile(output_path, mask_and_scale=False)
        self.assertEqual(fill_value, ds_out[new_param_name].values)

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
        pfts_to_include = ["not_vegetated", "needleleaf_evergreen_temperate_tree"]
        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            ",".join(pfts_to_include),
            "--drop-other-pfts",
            "xl=nan,nan",
            "planting_temp=nan,nan",
        ]

        # Check that planting_temp is already nan and xl isn't
        ds_in = open_paramfile(PARAMFILE, mask_and_scale=True)
        self.assertTrue(all(np.isnan(ds_in["planting_temp"].isel(pft=[0, 1]))))
        self.assertTrue(not any(np.isnan(ds_in["xl"].isel(pft=[0, 1]))))

        sp.main()
        self.assertTrue(os.path.exists(output_path))
        ds_out = open_paramfile(output_path, mask_and_scale=True)
        self.assertTrue(all(np.isnan(ds_out["xl"])))
        self.assertTrue(all(np.isnan(ds_out["planting_temp"])))

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
            NotImplementedError, "Can't set parameter to fill value if it doesn't already have one:"
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
            NotImplementedError, "Can't set integer parameter to fill value:"
        ):
            sp.main()

    # TODO: Test changing param when extracting just one PFT

    # TODO: Test changing PFT name

    def test_set_paramfile_changeparam_multidim_errors(self):
        """
        Test that set_paramfile errors if requesting change of a multi-dimensional parameter. This
        test will obviously need to be replaced once that functionality is added.
        """
        output_path = os.path.join(self.tempdir, "output.nc")
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

    def test_set_paramfile_setparams_just_one_pft(self):
        """Test changing just one PFT's value of something without dropping others"""
        output_path = os.path.join(self.tempdir, "output.nc")
        pft_to_include = "needleleaf_deciduous_boreal_tree"
        this_var = "rswf_max"
        new_value = 0.7

        # Ensure it wasn't new_value before
        ds_in = open_paramfile(PARAMFILE)
        pft_names = check_pfts_in_paramfile([pft_to_include], ds_in)
        pft_index = get_selected_pft_indices([pft_to_include], pft_names)[0]
        self.assertFalse(ds_in[this_var].values[pft_index] == new_value)

        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            pft_to_include,
            f"{this_var}={new_value}",
        ]
        sp.main()

        ds_out = open_paramfile(output_path)
        for i, value in enumerate(ds_out[this_var]):
            if i == pft_index:
                self.assertTrue(value == new_value)
            else:
                self.assertTrue(value == ds_in[this_var].values[i])

    def test_set_paramfile_setparams_just_one_pft_dropothers_noset(self):
        """Test dropping all but one PFT without changing any parameters"""
        output_path = os.path.join(self.tempdir, "output.nc")
        pft_to_include = "needleleaf_deciduous_boreal_tree"

        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            pft_to_include,
            "--drop-other-pfts",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))

        # Check that the file is just what you get if you drop all but the one PFT
        ds_in = open_paramfile(PARAMFILE)
        ds_in_1pft = sp.drop_other_pfts([pft_to_include], ds_in)
        ds_out = open_paramfile(output_path)
        self.assertTrue(set(ds_in_1pft.variables) == set(ds_out.variables))
        for var in ds_in_1pft:
            self.assertTrue(are_paramfile_dataarrays_identical(ds_in_1pft[var], ds_out[var]))
        self.assertTrue(ds_in_1pft.equals(ds_out))
        self.assertEqual(ds_in_1pft.sizes["pft"], 1)
        self.assertEqual(ds_out.sizes["pft"], 1)

    def test_set_paramfile_setparams_just_one_pft_dropothers_doset(self):
        """Test dropping all but one PFT, changing one parameter"""
        output_path = os.path.join(self.tempdir, "output.nc")
        pft_to_include = "needleleaf_deciduous_boreal_tree"
        this_var = "rswf_max"
        new_value = 0.7

        # Ensure it wasn't new_value before
        ds_in = open_paramfile(PARAMFILE)
        pft_names = check_pfts_in_paramfile([pft_to_include], ds_in)
        pft_index = get_selected_pft_indices([pft_to_include], pft_names)[0]
        self.assertFalse(ds_in[this_var].values[pft_index] == new_value)

        sys.argv = [
            "set_paramfile",
            "-i",
            PARAMFILE,
            "-o",
            output_path,
            "-p",
            pft_to_include,
            "--drop-other-pfts",
            f"{this_var}={new_value}",
        ]
        sp.main()
        self.assertTrue(os.path.exists(output_path))

        # Check that all variables match except for the one we changed
        ds_out = open_paramfile(output_path)
        ds_in_1pft = sp.drop_other_pfts([pft_to_include], ds_in)
        for var in ds_in.variables:
            da_in = ds_in_1pft[var]
            da_out = ds_out[var]
            if var == this_var:
                self.assertFalse(are_paramfile_dataarrays_identical(da_in, da_out))
            else:
                self.assertTrue(are_paramfile_dataarrays_identical(da_in, da_out))

    def test_set_paramfile_int_errors_given_float_point0(self):
        """
        Test that set_paramfile errors if given float value ending in .0 for an int field
        """
        input_path = os.path.join(self.tempdir, "input.nc")
        output_path = os.path.join(self.tempdir, "output.nc")
        param_name = "param0"

        # Save a test paramfile with an int field
        da = xr.DataArray(data=np.int32(3))
        ds = xr.Dataset(data_vars={param_name: da})
        ds.to_netcdf(input_path, encoding={param_name: {"_FillValue": None}})
        ds_in = open_paramfile(input_path)
        self.assertTrue(sp.is_integer(ds_in[param_name].values))

        # Try giving it a value ending in .0
        sys.argv = ["set_paramfile", "-i", input_path, "-o", output_path, f"{param_name}=4.0"]
        with self.assertRaisesRegex(ValueError, "Invalid assignment to an integer parameter:"):
            sp.main()

    def test_set_paramfile_int_errors_given_float_point1(self):
        """
        Test that set_paramfile errors if given float value NOT ending in .0 for an int field
        """
        input_path = os.path.join(self.tempdir, "input.nc")
        output_path = os.path.join(self.tempdir, "output.nc")
        param_name = "param0"

        # Save a test paramfile with an int field
        da = xr.DataArray(data=np.int32(3))
        ds = xr.Dataset(data_vars={param_name: da})
        ds.to_netcdf(input_path, encoding={param_name: {"_FillValue": None}})
        ds_in = open_paramfile(input_path)
        self.assertTrue(sp.is_integer(ds_in[param_name].values))

        # Try giving it a value ending in .1
        sys.argv = ["set_paramfile", "-i", input_path, "-o", output_path, f"{param_name}=4.1"]
        with self.assertRaisesRegex(ValueError, "Invalid assignment to an integer parameter:"):
            sp.main()

    def test_set_paramfile_double_ok_given_int(self):
        """
        Test that set_paramfile works if given int value for a double field
        """
        input_path = os.path.join(self.tempdir, "input.nc")
        output_path = os.path.join(self.tempdir, "output.nc")
        param_name = "param0"

        # Save a test paramfile with a double field
        da = xr.DataArray(data=np.float32(3.14))
        ds = xr.Dataset(data_vars={param_name: da})
        ds.to_netcdf(input_path, encoding={param_name: {"_FillValue": None}})
        ds_in = open_paramfile(input_path)
        self.assertFalse(sp.is_integer(ds_in[param_name].values))

        # Give it an integer
        sys.argv = ["set_paramfile", "-i", input_path, "-o", output_path, f"{param_name}=4"]
        sp.main()

        # Check that it's still a double after saving
        ds_out = open_paramfile(output_path)
        self.assertFalse(sp.is_integer(ds_out[param_name].values))


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
