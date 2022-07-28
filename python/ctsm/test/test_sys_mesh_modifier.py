#!/usr/bin/env python3

"""System tests for mesh_mask_modifier"""

import os
import re

import unittest
import tempfile
import shutil

import xarray as xr

from ctsm.path_utils import path_to_ctsm_root
from ctsm import unit_testing
from ctsm.modify_mesh_mask.mesh_mask_modifier import mesh_mask_modifier

# Allow test names that pylint doesn't like; otherwise hard to make them
# readable
# pylint: disable=invalid-name


class TestSysMeshMaskModifier(unittest.TestCase):
    """System tests for mesh_mask_modifier"""

    def setUp(self):
        """
        Obtain path to the existing:
        - modify_template.cfg file
        - /testinputs directory and fsurdat_in, located in /testinputs
        Make /_tempdir for use by these tests.
        Obtain path and names for the files being created in /_tempdir:
        - modify_mesh_mask.cfg
        - mesh_mask_out.nc
        Generate mesh_mask_in file using nco/esmf commands on fsurdat_in.
        Generate landmask_file using nco commands on fsurdat_in.
        """
        self._cfg_template_path = os.path.join(
            path_to_ctsm_root(), "tools/modify_mesh_mask/modify_template.cfg"
        )
        testinputs_path = os.path.join(path_to_ctsm_root(), "python/ctsm/test/testinputs")
        self._fsurdat_in = os.path.join(
            testinputs_path,
            "surfdata_5x5_amazon_16pfts_Irrig_CMIP6_simyr2000_c171214.nc",
        )
        self._tempdir = tempfile.mkdtemp()
        self._cfg_file_path = os.path.join(self._tempdir, "modify_mesh_mask.cfg")
        self._mesh_mask_out = os.path.join(self._tempdir, "mesh_mask_out.nc")
        #TODO Generate mesh_mask_in file using nco/esmf commands on fsurdat_in
        #TODO Generate landmask_file using nco commands on fsurdat_in

    def tearDown(self):
        """
        Remove temporary directory
        """
        shutil.rmtree(self._tempdir, ignore_errors=True)

    def test_allInfo(self):
        """
        This test specifies all the information that one may specify
        Create .cfg file, run the tool, compare mesh_mask_in to mesh_mask_out
        """

        self._create_config_file()

        # run the mesh_mask_modifier tool
        mesh_mask_modifier(self._cfg_file_path)
        # the critical piece of this test is that the above command
        # doesn't generate errors; however, we also do some assertions below

        mesh_mask_in_data = xr.open_dataset(self._mesh_mask_in)
        mesh_mask_out_data = xr.open_dataset(self._mesh_mask_out)
        # assert that mesh_mask_out equals mesh_mask_in
        # TODO Won't work because the Mask variable will now = zeros, not ones
        self.assertTrue(mesh_mask_out_data.equals(mesh_mask_in_data))

    def _create_config_file(self):
        """
        Open the new and the template .cfg files
        Loop line by line through the template .cfg file
        When string matches, replace that line's content
        """
        with open(self._cfg_file_path, "w", encoding="utf-8") as cfg_out:
            with open(self._cfg_template_path, "r", encoding="utf-8") as cfg_in:
                for line in cfg_in:
                    if re.match(r" *mesh_mask_in *=", line):
                        line = f"mesh_mask_in = {self._mesh_mask_in}"
                    elif re.match(r" *mesh_mask_out *=", line):
                        line = f"mesh_mask_out = {self._mesh_mask_out}"
                    elif re.match(r" *landmask_file *=", line):
                        line = f"landmask_file = {self._landmask_file}"
                    elif re.match(r" *lat_dimname *=", line):
                        line = f"lat_dimname = {self._lat_dimname}"
                    elif re.match(r" *lon_dimname *=", line):
                        line = f"lon_dimname = {self._lon_dimname}"
                    elif re.match(r" *lat_varname *=", line):
                        line = f"lat_varname = {self._lat_varname}"
                    elif re.match(r" *lon_varname *=", line):
                        line = f"lon_varname = {self._lon_varname}"
                    cfg_out.write(line)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
