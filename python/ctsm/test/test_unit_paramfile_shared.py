#!/usr/bin/env python3

"""Unit tests for paramfile_shared"""

import unittest

from ctsm import unit_testing

from ctsm.param_utils import paramfile_shared as ps

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name


class TestUnitGetSelectedPftIndices(unittest.TestCase):
    """Unit tests of get_selected_pft_indices"""

    def test_get_selected_pft_indices_1strselected_onlyinlist(self):
        """Check get_selected_pft_indices() given the only one in the list, as a string"""
        selected_pfts = "rice"
        pft_names = ["rice"]
        result = ps.get_selected_pft_indices(selected_pfts=selected_pfts, pft_names=pft_names)
        self.assertListEqual(result, [0])

    def test_get_selected_pft_indices_1selected_onlyinlist(self):
        """Check get_selected_pft_indices() given the only one in the list, as a list"""
        selected_pfts = ["rice"]
        pft_names = ["rice"]
        result = ps.get_selected_pft_indices(selected_pfts=selected_pfts, pft_names=pft_names)
        self.assertListEqual(result, [0])

    def test_get_selected_pft_indices_2selected_sameorder(self):
        """Check get_selected_pft_indices() given 2 selected in the same order as the list"""
        pft_names = ["rice", "irrigated_rice"]
        result = ps.get_selected_pft_indices(selected_pfts=pft_names, pft_names=pft_names)
        self.assertListEqual(result, [0, 1])

    def test_get_selected_pft_indices_2selected_difforder(self):
        """Check get_selected_pft_indices() given 2 selected NOT in the same order as the list"""
        pft_names = ["rice", "irrigated_rice"]
        result = ps.get_selected_pft_indices(
            selected_pfts=list(reversed(pft_names)), pft_names=pft_names
        )
        self.assertListEqual(result, [1, 0])

    def test_get_selected_pft_indices_missing_valueerror(self):
        """Check get_selected_pft_indices() given selected pft NOT in the list"""
        selected_pfts = ["wheat"]
        pft_names = ["rice", "irrigated_rice"]
        with self.assertRaises(ValueError):
            ps.get_selected_pft_indices(selected_pfts=selected_pfts, pft_names=pft_names)


if __name__ == "__main__":
    unit_testing.setup_for_tests()
    unittest.main()
