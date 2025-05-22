#!/usr/bin/env python3
"""
Unit tests for creating and modifying surface datasets in SinglePointCase

for the rest of SinglePointCase tests please see : test_unit_singlept_data

You can run this by:
    python -m unittest test_unit_singlept_data_surfdata.py
"""

import unittest
import os
import sys


import numpy as np
import xarray as xr

# -- add python/ctsm  to path (needed if we want to run the test stand-alone)
_CTSM_PYTHON = os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, os.pardir)
sys.path.insert(1, _CTSM_PYTHON)

# pylint: disable=wrong-import-position
from ctsm import unit_testing
from ctsm.site_and_regional.single_point_case import SinglePointCase

# pylint: disable=invalid-name
# pylint: disable=too-many-lines


class TestSinglePointCaseSurfaceNoCrop(unittest.TestCase):
    """
    Basic class for testing creating and modifying surface dataset for
    non-crop cases (aka using 16 pft dataset) in SinglePointCase class in single_point_case.py.

    """

    plat = 20.1
    plon = 50.5
    site_name = None
    create_domain = True
    create_surfdata = True
    create_landuse = True
    create_datm = True
    create_user_mods = True
    dom_pft = [8]
    evenly_split_cropland = False
    pct_pft = None
    num_pft = 16
    cth = 0.9
    cbh = 0.1
    include_nonveg = False
    uni_snow = True
    cap_saturation = True
    out_dir = os.getcwd()
    overwrite = False

    # -- dimensions of xarray dataset
    lsmlat = [plat]
    lsmlon = [plon]
    months = np.arange(1, 13, 1, dtype=int)
    lsmpft = np.arange(0, 79, 1, dtype=int)
    natpft = np.arange(0, 15, 1, dtype=int)
    cft = np.arange(15, 17, 1, dtype=int)
    numurbl = np.arange(0, 3, 1, dtype=int)

    ds_test = xr.Dataset(
        {
            "PCT_NATVEG": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={
                    "long_name": "total percent natural vegetation landunit",
                    "units": "unitless",
                },
            ),
            "PCT_CROP": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "total percent crop landunit", "units": "unitless"},
            ),
            "PCT_NAT_PFT": xr.DataArray(
                data=np.random.rand(1, 1, 15),
                dims=["lsmlat", "lsmlon", "natpft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "natpft": natpft},
                attrs={
                    "long_name": "percent plant functional type on the natural veg landunit",
                    "units": "unitless",
                },
            ),
            "PCT_CFT": xr.DataArray(
                data=np.random.rand(1, 1, 2),
                dims=["lsmlat", "lsmlon", "cft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "cft": cft},
                attrs={
                    "long_name": "percent crop functional type on the crop landunit",
                    "units": "unitless",
                },
            ),
            "PCT_LAKE": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent lake", "units": "unitless"},
            ),
            "PCT_WETLAND": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent wetland", "units": "unitless"},
            ),
            "PCT_OCEAN": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent ocean", "units": "unitless"},
            ),
            "PCT_URBAN": xr.DataArray(
                data=np.random.rand(1, 1, 3),
                dims=["lsmlat", "lsmlon", "numurbl"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "numurbl": numurbl},
                attrs={
                    "long_name": "percent urban for each density type",
                    "units": "unitless",
                },
            ),
            "PCT_GLACIER": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent glacier", "units": "unitless"},
            ),
            "STD_ELEV": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "standard deviation of elevation", "units": "m"},
            ),
            "FMAX": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={
                    "long_name": "maximum fractional saturated area",
                    "units": "unitless",
                },
            ),
            "MONTHLY_HEIGHT_TOP": xr.DataArray(
                data=np.random.rand(1, 1, months[-1], lsmpft[-1] + 1),
                dims=["lsmlat", "lsmlon", "time", "lsmpft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "time": months, "lsmpft": lsmpft},
                attrs={
                    "long_name": "monthly height top by pft and month",
                    "units": "m",
                },
            ),
            "MONTHLY_HEIGHT_BOT": xr.DataArray(
                data=np.random.rand(1, 1, months[-1], lsmpft[-1] + 1),
                dims=["lsmlat", "lsmlon", "time", "lsmpft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "time": months, "lsmpft": lsmpft},
                attrs={
                    "long_name": "monthly height bottom by pft and month",
                    "units": "m",
                },
            ),
        },
        attrs={"Conventions": "test data only"},
    )

    def test_modify_surfdata_atpoint_nocrop_1pft_pctnatpft(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NAT_PFT for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [5]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.zeros((1, 1, 15))
        expected_out[:, :, 5] = 100

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_NAT_PFT"].data, expected_out)

    def test_modify_surfdata_atpoint_nocrop_1pft_pctnatveg(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NATVEG for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [5]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_NATVEG"].data[:, :], 100)

    def test_modify_surfdata_atpoint_nocrop_1pft_pctcrop(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_CROP for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [5]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_CROP"].data[:, :], 0)

    def test_modify_surfdata_atpoint_nocrop_1pft_glacier(self):
        """
        Test modify_surfdata_atpoint
        Checks GLACIER for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = False
        single_point.dom_pft = [5]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_GLACIER"].data[:, :], 0)

    def test_modify_surfdata_atpoint_nocrop_1pft_wetland(self):
        """
        Test modify_surfdata_atpoint
        Checks WETLAND for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [5]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_WETLAND"].data[:, :], 0)

    def test_modify_surfdata_atpoint_nocrop_1pft_lake(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_LAKE for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [5]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_LAKE"].data[:, :], 0)

    def test_modify_surfdata_atpoint_nocrop_1pft_unisnow(self):
        """
        Test modify_surfdata_atpoint
        Checks STD_ELV for one pft and unisnow
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = False
        single_point.dom_pft = [5]
        single_point.uni_snow = True
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["STD_ELEV"].data[:, :], 20)

    def test_modify_surfdata_atpoint_nocrop_1pft_capsat(self):
        """
        Test modify_surfdata_atpoint
        Checks FMAX for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = False
        single_point.dom_pft = [5]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)
        single_point.cap_saturation = True

        self.assertEqual(ds_out["FMAX"].data[:, :], 0)

    def test_modify_surfdata_atpoint_nocrop_multipft(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NAT_PFT for multi pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = False
        single_point.dom_pft = [1, 3, 5]
        single_point.pct_pft = [0.5, 0.4, 0.1]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.zeros((1, 1, 15))
        expected_out[:, :, 1] = 0.5
        expected_out[:, :, 3] = 0.4
        expected_out[:, :, 5] = 0.1

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_NAT_PFT"].data, expected_out)

    def test_modify_surfdata_atpoint_nocrop_urban_nononveg(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks URBAN for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = False
        single_point.dom_pft = [7]
        single_point.plat = [34.05]
        single_point.plon = [118.25]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.zeros((1, 1, 3))

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_URBAN"].data, expected_out)

    def test_modify_surfdata_atpoint_nocrop_urban_include_nonveg(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks URBAN for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = True
        single_point.dom_pft = [7]
        single_point.plat = [34.05]
        single_point.plon = [118.25]

        # -- change it to something known
        self.ds_test["PCT_URBAN"][:, :, :] = 1
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.ones((1, 1, 3))

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_URBAN"].data, expected_out)

    def test_modify_surfdata_atpoint_nocrop_wetland_include_nonveg(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks PCT_WETLAND for one pft to make sure it is not zerod-out
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = True
        single_point.dom_pft = [7]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertNotEqual(ds_out["PCT_WETLAND"].data[:, :], 0)

    def test_modify_surfdata_atpoint_nocrop_nopft_zero_nonveg(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_CROP for no pft and zero nonveg
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = None
        single_point.include_nonveg = False
        self.ds_test["PCT_CROP"].values = [[40]]
        self.ds_test["PCT_LAKE"].values = [[10]]
        self.ds_test["PCT_WETLAND"].values = [[10]]
        self.ds_test["PCT_NATVEG"].values = [[40]]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_CROP"].data[:, :], 50)

    def test_modify_surfdata_atpoint_nocrop_nopft_include_nonveg(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_CROP for no pft and include nonveg
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = None
        single_point.include_nonveg = True
        self.ds_test["PCT_CROP"].values = [[40]]
        self.ds_test["PCT_LAKE"].values = [[10]]
        self.ds_test["PCT_WETLAND"].values = [[10]]
        self.ds_test["PCT_NATVEG"].values = [[40]]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_CROP"].data[:, :], 40)


class TestSinglePointCaseSurfaceCrop(unittest.TestCase):
    """
    Basic class for testing creating and modifying surface dataset for
    crop cases (aka using 78 pft dataset) in SinglePointCase class in single_point_case.py.
    """

    plat = 20.1
    plon = 50.5
    site_name = None
    create_domain = True
    create_surfdata = True
    create_landuse = True
    create_datm = True
    create_user_mods = True
    dom_pft = [17]
    evenly_split_cropland = False
    pct_pft = None
    num_pft = 78
    cth = 0.9
    cbh = 0.1
    include_nonveg = False
    uni_snow = False
    cap_saturation = False
    out_dir = os.getcwd()
    overwrite = False

    # -- dimensions of xarray dataset
    lsmlat = [plat]
    lsmlon = [plon]
    months = np.arange(1, 12, 1, dtype=int)
    lsmpft = np.arange(0, 79, 1, dtype=int)
    natpft = np.arange(0, 15, 1, dtype=int)
    cft = np.arange(15, 79, 1, dtype=int)
    numurbl = np.arange(0, 3, 1, dtype=int)

    ds_test = xr.Dataset(
        {
            "PCT_NATVEG": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={
                    "long_name": "total percent natural vegetation landunit",
                    "units": "unitless",
                },
            ),
            "PCT_CROP": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "total percent crop landunit", "units": "unitless"},
            ),
            "PCT_NAT_PFT": xr.DataArray(
                data=np.random.rand(1, 1, 15),
                dims=["lsmlat", "lsmlon", "natpft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "natpft": natpft},
                attrs={
                    "long_name": "percent plant functional type on the natural veg landunit",
                    "units": "unitless",
                },
            ),
            "PCT_CFT": xr.DataArray(
                data=np.random.rand(1, 1, 64),
                dims=["lsmlat", "lsmlon", "cft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "cft": cft},
                attrs={
                    "long_name": "percent crop functional type on the crop landunit",
                    "units": "unitless",
                },
            ),
            "PCT_LAKE": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent lake", "units": "unitless"},
            ),
            "PCT_WETLAND": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent wetland", "units": "unitless"},
            ),
            "PCT_OCEAN": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent ocean", "units": "unitless"},
            ),
            "PCT_URBAN": xr.DataArray(
                data=np.random.rand(1, 1, 3),
                dims=["lsmlat", "lsmlon", "numurbl"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "numurbl": numurbl},
                attrs={
                    "long_name": "percent urban for each density type",
                    "units": "unitless",
                },
            ),
            "PCT_GLACIER": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "percent glacier", "units": "unitless"},
            ),
            "STD_ELEV": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={"long_name": "standard deviation of elevation", "units": "m"},
            ),
            "FMAX": xr.DataArray(
                data=np.random.rand(1, 1),
                dims=["lsmlat", "lsmlon"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon},
                attrs={
                    "long_name": "maximum fractional saturated area",
                    "units": "unitless",
                },
            ),
            "MONTHLY_HEIGHT_TOP": xr.DataArray(
                data=np.random.rand(1, 1, months[-1], lsmpft[-1] + 1),
                dims=["lsmlat", "lsmlon", "time", "lsmpft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "time": months, "lsmpft": lsmpft},
                attrs={
                    "long_name": "monthly height top by pft and month",
                    "units": "m",
                },
            ),
            "MONTHLY_HEIGHT_BOT": xr.DataArray(
                data=np.random.rand(1, 1, months[-1], lsmpft[-1] + 1),
                dims=["lsmlat", "lsmlon", "time", "lsmpft"],
                coords={"lsmlat": lsmlat, "lsmlon": lsmlon, "time": months, "lsmpft": lsmpft},
                attrs={
                    "long_name": "monthly height bottom by pft and month",
                    "units": "m",
                },
            ),
        },
        attrs={"Conventions": "test data only"},
    )

    def test_modify_surfdata_atpoint_crop_1pft_pctnatpft(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NAT_PFT for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [19]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.zeros((1, 1, 64))
        expected_out[:, :, 4] = 100

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_CFT"].data, expected_out)

    def test_modify_surfdata_atpoint_crop_1pft_pctnatveg(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NATVEG for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_NATVEG"].data[:, :], 0)

    def test_modify_surfdata_atpoint_crop_1pft_pctcrop(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_CROP for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_CROP"].data[:, :], 100)

    def test_modify_surfdata_atpoint_crop_1pft_glacier(self):
        """
        Test modify_surfdata_atpoint
        Checks GLACIER for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_GLACIER"].data[:, :], 0)

    def test_modify_surfdata_atpoint_crop_1pft_wetland(self):
        """
        Test modify_surfdata_atpoint
        Checks WETLAND for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_WETLAND"].data[:, :], 0)

    def test_modify_surfdata_atpoint_crop_1pft_lake(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_LAKE for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_LAKE"].data[:, :], 0)

    def test_modify_surfdata_atpoint_crop_1pft_unisnow(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks STD_ELV for one pft and unisnow
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17]
        single_point.uni_snow = True
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["STD_ELEV"].data[:, :], 20)

    def test_modify_surfdata_atpoint_crop_1pft_capsat(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks FMAX for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.cap_saturation = True
        single_point.dom_pft = [22]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)
        single_point.cap_saturation = True

        self.assertEqual(ds_out["FMAX"].data[:, :], 0)

    def test_modify_surfdata_atpoint_crop_multipft(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks PCT_NAT_PFT for multi pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = [17, 22]
        single_point.pct_pft = [0.6, 0.4]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.zeros((1, 1, 64))
        expected_out[:, :, 2] = 0.6
        expected_out[:, :, 7] = 0.4

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_CFT"].data, expected_out)

    def test_modify_surfdata_atpoint_crop_urban_nononveg(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks URBAN for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = False
        single_point.dom_pft = [17]
        single_point.plat = [34.05]
        single_point.plon = [118.25]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.zeros((1, 1, 3))

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_URBAN"].data, expected_out)

    def test_modify_surfdata_atpoint_crop_urban_include_nonveg(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks URBAN for one pft
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = True
        single_point.dom_pft = [17]
        single_point.plat = [34.05]
        single_point.plon = [118.25]

        # -- change it to something known
        self.ds_test["PCT_URBAN"][:, :, :] = 1
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        expected_out = np.ones((1, 1, 3))

        # self.assertEqual(ds_out['PCT_NAT_PFT'].data[:,:,5], 100)
        np.testing.assert_array_equal(ds_out["PCT_URBAN"].data, expected_out)

    def test_modify_surfdata_atpoint_crop_lake_include_nonveg(self):
        """
        Test modify_surfdata_atpoint for crop cases
        Checks PCT_LAKE for one pft to make sure it is not zerod-out
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.include_nonveg = True
        single_point.dom_pft = [17]

        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertNotEqual(ds_out["PCT_LAKE"].data[:, :], 0)

    def test_modify_surfdata_atpoint_crop_nopft_zero_nonveg(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NATVEG for no pft and zero nonveg
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = None
        single_point.include_nonveg = False
        self.ds_test["PCT_CROP"].values = [[40]]
        self.ds_test["PCT_LAKE"].values = [[10]]
        self.ds_test["PCT_WETLAND"].values = [[10]]
        self.ds_test["PCT_NATVEG"].values = [[40]]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_NATVEG"].data[:, :], 50)

    def test_modify_surfdata_atpoint_crop_nopft_include_nonveg(self):
        """
        Test modify_surfdata_atpoint
        Checks PCT_NATVEG for no pft and include nonveg
        """
        single_point = SinglePointCase(
            plat=self.plat,
            plon=self.plon,
            site_name=self.site_name,
            create_domain=self.create_domain,
            create_surfdata=self.create_surfdata,
            create_landuse=self.create_landuse,
            create_datm=self.create_datm,
            create_user_mods=self.create_user_mods,
            dom_pft=self.dom_pft,
            evenly_split_cropland=self.evenly_split_cropland,
            pct_pft=self.pct_pft,
            num_pft=self.num_pft,
            cth=self.cth,
            cbh=self.cbh,
            include_nonveg=self.include_nonveg,
            uni_snow=self.uni_snow,
            cap_saturation=self.cap_saturation,
            out_dir=self.out_dir,
            overwrite=self.overwrite,
        )
        single_point.dom_pft = None
        single_point.include_nonveg = True
        self.ds_test["PCT_CROP"].values = [[40]]
        self.ds_test["PCT_LAKE"].values = [[10]]
        self.ds_test["PCT_WETLAND"].values = [[10]]
        self.ds_test["PCT_NATVEG"].values = [[40]]
        ds_out = single_point.modify_surfdata_atpoint(self.ds_test)

        self.assertEqual(ds_out["PCT_NATVEG"].data[:, :], 40)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
