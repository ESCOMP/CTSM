"""System tests of functions in namelist_utils module (anything touching filesystem)"""

# pylint: disable=too-few-public-methods,protected-access

from pathlib import Path
import xml.etree.ElementTree as ET
from unittest.mock import patch

import pytest

from ctsm.no_nans_in_inputs.namelist_utils import (
    _extract_file_path_list_from_usernl,
    _extract_file_path_set_from_usernl,
    _extract_file_paths_from_xml,
    extract_file_paths_from_file,
    find_user_nl_files,
    _update_usernl_file,
    _update_xml_file,
    update_text_file_referencing_netcdf,
)
from ctsm.no_nans_in_inputs import namelist_utils

# Test constants
TEST_PARAM_CLM60 = "lnd/clm2/paramdata/ctsm60_params.c260108.nc"
TEST_PARAM_CLM50 = "lnd/clm2/paramdata/clm50_params.c250311.nc"
TEST_PARAM_CLM45 = "lnd/clm2/paramdata/clm45_params.c250311.nc"
TEST_PATH_PARAM = "lnd/clm2/paramdata/test_params.nc"
TEST_PATH_SURF = "lnd/clm2/surfdata/test_surf.nc"
TEST_PATH_INIT = "lnd/clm2/initdata/test_init.nc"
TEST_PATH_OTHER = "share/meshes/test_mesh.nc"
TEST_PHYS_CLM60 = "clm6_0"
TEST_PHYS_CLM50 = "clm5_0"
TEST_PHYS_CLM45 = "clm4_5"


class TestUpdateXmlFile:
    """Test the _update_xml_file function."""

    @pytest.mark.parametrize("fn_to_test", [_update_xml_file, update_text_file_referencing_netcdf])
    def test_update_xml_path(self, mock_xml_file_path, fn_to_test):
        """Test updating a file path in XML."""
        # Create XML content in the mocked path
        xml_content = """<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>lnd/clm2/paramdata/test_params.nc</paramfile>
    <surfdata>lnd/clm2/surfdata/test_surf.nc</surfdata>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        # Update one of the paths
        old_path = "lnd/clm2/paramdata/test_params.nc"
        new_path = "lnd/clm2/paramdata/test_params.no_nan_fill.nc"

        # Spy on _check_xml_file: Mock it so that it behaves as usual but we can check call count
        with patch(
            "ctsm.no_nans_in_inputs.namelist_utils._check_xml_file",
            wraps=namelist_utils._check_xml_file,
        ) as mock_check_xml:
            fn_to_test(mock_xml_file_path, old_path, new_path)

        # Read and verify the updated XML
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()
        paramfile = root.find("paramfile")
        assert paramfile is not None
        assert paramfile.text == new_path

        # Other paths should be unchanged
        surfdata = root.find("surfdata")
        assert surfdata.text == "lnd/clm2/surfdata/test_surf.nc"

        # We should have checked the file
        assert mock_check_xml.call_count == 1

    def test_update_xml_path_not_found(self, mock_xml_file_path):
        """Test that updating non-existent path raises ValueError."""
        xml_content = """<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>lnd/clm2/paramdata/test_params.nc</paramfile>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        with pytest.raises(ValueError, match="not found"):
            _update_xml_file(mock_xml_file_path, "nonexistent/path.nc", "new/path.nc")

    def test_update_xml_with_multiple_same_tag(self, mock_xml_file_path):
        """Test updating path when multiple elements have the same tag name."""
        # Simulate the real XML structure with multiple paramfile elements
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="{TEST_PHYS_CLM60}">{TEST_PARAM_CLM60}</paramfile>
    <paramfile phys="{TEST_PHYS_CLM50}">{TEST_PARAM_CLM50}</paramfile>
    <paramfile phys="{TEST_PHYS_CLM45}">{TEST_PARAM_CLM45}</paramfile>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        # Update one specific path
        old_path = TEST_PARAM_CLM60
        new_path = TEST_PARAM_CLM60.replace(".nc", ".no_nan_fill.nc")

        _update_xml_file(mock_xml_file_path, old_path, new_path)

        # Read and verify the updated XML
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()

        # Find all paramfile elements
        paramfiles = root.findall("paramfile")
        assert len(paramfiles) == 3

        # Check that only the clm6_0 one was updated
        clm60_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM60][0]
        assert clm60_param.text == new_path

        # Check that others are unchanged
        clm50_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM50][0]
        assert clm50_param.text == TEST_PARAM_CLM50

        clm45_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM45][0]
        assert clm45_param.text == TEST_PARAM_CLM45

    def test_update_xml_replaces_within_element_text(self, mock_xml_file_path):
        """Test that all occurrences within a single element's text are replaced."""
        # Create XML with a path appearing multiple times in one element's text
        test_path = "lnd/clm2/test/file.nc"
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <multi_path>{test_path} {test_path}</multi_path>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        new_path = test_path.replace(".nc", ".no_nan_fill.nc")
        _update_xml_file(mock_xml_file_path, test_path, new_path)

        # Verify both occurrences were replaced
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()
        multi_path = root.find("multi_path")
        assert multi_path is not None
        assert multi_path.text == f"{new_path} {new_path}"
        # Verify old path is completely gone
        assert test_path not in multi_path.text

    def test_update_xml_replaces_across_different_tags(self, mock_xml_file_path):
        """Test that same path in different element types are all replaced."""
        # Create XML with same path in multiple different element types
        test_path = "lnd/clm2/test/shared_file.nc"
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{test_path}</paramfile>
    <surfdata>{test_path}</surfdata>
    <initdata>{test_path}</initdata>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        new_path = test_path.replace(".nc", ".no_nan_fill.nc")
        _update_xml_file(mock_xml_file_path, test_path, new_path)

        # Verify all three elements were updated
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()

        paramfile = root.find("paramfile")
        assert paramfile is not None
        assert paramfile.text == new_path

        surfdata = root.find("surfdata")
        assert surfdata is not None
        assert surfdata.text == new_path

        initdata = root.find("initdata")
        assert initdata is not None
        assert initdata.text == new_path

    def test_update_xml_replaces_across_same_tag_different_attrs(self, mock_xml_file_path):
        """Test that same path in elements with same tag but different attributes are all replaced."""
        # Simulate scenario where two paramfile elements with different attributes point to same file
        # (like if lines 625 and 626 in the real XML both had the same path)
        test_path = "lnd/clm2/paramdata/shared_params.nc"
        xml_content = f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile phys="{TEST_PHYS_CLM60}">{test_path}</paramfile>
    <paramfile phys="{TEST_PHYS_CLM50}">{test_path}</paramfile>
    <surfdata>lnd/clm2/surfdata/different_file.nc</surfdata>
</namelist_defaults>
"""
        with open(mock_xml_file_path, "w", encoding="utf-8") as f:
            f.write(xml_content)

        new_path = test_path.replace(".nc", ".no_nan_fill.nc")
        _update_xml_file(mock_xml_file_path, test_path, new_path)

        # Verify both paramfile elements were updated
        tree = ET.parse(mock_xml_file_path)
        root = tree.getroot()

        paramfiles = root.findall("paramfile")
        assert len(paramfiles) == 2

        # Both should have the new path
        clm60_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM60][0]
        assert clm60_param.text == new_path

        clm50_param = [p for p in paramfiles if p.get("phys") == TEST_PHYS_CLM50][0]
        assert clm50_param.text == new_path

        # Surfdata should be unchanged
        surfdata = root.find("surfdata")
        assert surfdata.text == "lnd/clm2/surfdata/different_file.nc"


class TestUpdateUsernlFile:
    """Test the _update_usernl_file function."""

    def _do_test(
        self,
        mock_usernl_file_path: Path,
        nc_path_in: str,
        fn_to_test: callable = _update_usernl_file,
    ):

        # Get line with nc_path before substitution
        line_before = None
        with open(mock_usernl_file_path, "r", encoding="utf8") as f:
            n_line = -1
            for line in f.readlines():
                n_line += 1
                if nc_path_in in line:
                    line_before = line
                    break
        if line_before is None:
            raise RuntimeError(f"No line found containing {nc_path_in=}")

        # Substitute
        nc_path_out = "abc123.nc"
        # Spy on _check_usernl_file: Mock it so that it behaves as usual but we can check call count
        with patch(
            "ctsm.no_nans_in_inputs.namelist_utils._check_usernl_file",
            wraps=namelist_utils._check_usernl_file,
        ) as mock_check_usernl:
            fn_to_test(mock_usernl_file_path, nc_path_in, nc_path_out)

        # Get line after substitution
        line_after = None
        with open(mock_usernl_file_path, "r", encoding="utf8") as f:
            n = -1
            for line_after in f.readlines():
                n += 1
                if n == n_line:
                    break
        if line_after is None:
            raise RuntimeError(f"No lines read from {mock_usernl_file_path=}")

        # Make sure *something* happened on the line
        assert line_before != line_after

        # Make sure the line looks like what we expect: The filename was replaced but nothing else
        # was touched.
        expected = line_before.replace(nc_path_in, nc_path_out)
        assert line_after == expected

        # We should have checked the file
        assert mock_check_usernl.call_count == 1

    @pytest.mark.parametrize(
        "fn_to_test", [_update_usernl_file, update_text_file_referencing_netcdf]
    )
    def test_update_usernl_file_relpath(self, create_mock_user_nl_file, fn_to_test):
        """Test _update_usernl_file with rel_path line of our test user_nl_ file"""
        mock_usernl_file_path, nc_paths = create_mock_user_nl_file()
        nc_path_in = nc_paths.rel_path
        self._do_test(mock_usernl_file_path, nc_path_in, fn_to_test)

    def test_update_usernl_file_abspath(self, create_mock_user_nl_file):
        """Test _update_usernl_file with abs_path line of our test user_nl_ file"""
        mock_usernl_file_path, nc_paths = create_mock_user_nl_file()
        nc_path_in = nc_paths.abs_path
        self._do_test(mock_usernl_file_path, nc_path_in)

    def test_update_usernl_file_abs_path_dinlocroot(self, create_mock_user_nl_file):
        """Test _update_usernl_file with abs_path_dinlocroot line of our test user_nl_ file"""
        mock_usernl_file_path, nc_paths = create_mock_user_nl_file()
        nc_path_in = nc_paths.abs_path_dinlocroot
        self._do_test(mock_usernl_file_path, nc_path_in)

    def test_update_usernl_file_abs_path_dinlocrootcurly(self, create_mock_user_nl_file):
        """Test _update_usernl_file with abs_path_dinlocrootcurly line of our test user_nl_ file"""
        mock_usernl_file_path, nc_paths = create_mock_user_nl_file()
        nc_path_in = nc_paths.abs_path_dinlocrootcurly
        self._do_test(mock_usernl_file_path, nc_path_in)


class TestExtractFilePathListFromUserNl:
    """Test the _extract_file_path_list_from_usernl function."""

    def test_extracts_multiple_paths(self, create_mock_user_nl_file):
        """
        Test extracting multiple file paths from user_nl file, both directly and via
        extract_file_paths_from_file.
        """

        mock_usernl_file_path, nc_paths = create_mock_user_nl_file()

        result = _extract_file_path_list_from_usernl(mock_usernl_file_path)
        assert nc_paths.rel_path in result
        assert nc_paths.abs_path in result
        assert nc_paths.abs_path_dinlocroot in result
        assert nc_paths.abs_path_dinlocrootcurly in result
        assert len(result) == 4


class TestExtractFilePathSetFromUserNl:
    """Test the _extract_file_path_set_from_usernl function."""

    @pytest.mark.parametrize(
        "func_to_test, exact",
        [
            (_extract_file_path_set_from_usernl, False),
            (extract_file_paths_from_file, False),
            (_extract_file_path_set_from_usernl, True),
            (extract_file_paths_from_file, True),
        ],
    )
    def test_extracts_multiple_paths(self, create_mock_user_nl_file, func_to_test, exact):
        """
        Test extracting multiple file paths from user_nl file, both directly and via
        extract_file_paths_from_file.
        """

        mock_usernl_file_path, nc_paths = create_mock_user_nl_file()

        result = func_to_test(mock_usernl_file_path, exact)
        set_with_exact_false = {nc_paths.rel_path, nc_paths.abs_path}
        if exact:
            assert result == (
                set_with_exact_false
                | {nc_paths.abs_path_dinlocroot, nc_paths.abs_path_dinlocrootcurly}
            )
        else:
            assert result == set_with_exact_false


class TestExtractFilePathsFromXml:
    """Test the _extract_file_paths_from_xml function."""

    @pytest.mark.parametrize(
        "func_to_test", [_extract_file_paths_from_xml, extract_file_paths_from_file]
    )
    def test_extracts_multiple_paths(self, create_mock_xml_file, func_to_test):
        """
        Test extracting multiple file paths from XML, both directly and via
        extract_file_paths_from_file.
        """
        xml_path = create_mock_xml_file(
            f"""<?xml version="1.0"?>
<namelist_defaults>
    <paramfile>{TEST_PATH_PARAM}</paramfile>
    <surfdata>{TEST_PATH_SURF}</surfdata>
    <initdata>{TEST_PATH_INIT}</initdata>
</namelist_defaults>
"""
        )
        result = func_to_test(xml_path)
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
        result = _extract_file_paths_from_xml(xml_path)
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
        result = _extract_file_paths_from_xml(xml_path)
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
        result = _extract_file_paths_from_xml(xml_path)
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
        result = _extract_file_paths_from_xml(xml_path)
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
        result = _extract_file_paths_from_xml(xml_path)
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
        result = _extract_file_paths_from_xml(xml_path)
        assert result == set()

    def test_file_not_found_exits(self, capsys):
        """Test that a missing XML file causes SystemExit."""
        nonexistent_path = "/nonexistent/file.xml"
        with pytest.raises(SystemExit):
            _extract_file_paths_from_xml(nonexistent_path)
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
        result = _extract_file_paths_from_xml(xml_path)
        assert result == {TEST_PATH_PARAM, TEST_PATH_SURF}


class TestFindUserNlFiles:
    """Tests of find_user_nl_files()"""

    def test_find_user_nl_files(self, tmp_path):
        """Test find_user_nl_files()"""
        found_toplevel = f"{tmp_path}/user_nl_clm"
        Path(found_toplevel).touch()

        # Create files and dirs
        found_secondlevel = f"{tmp_path}/some_dir/user_nl_something"
        Path(found_secondlevel).parent.mkdir()
        Path(found_secondlevel).touch()
        notfound = f"{tmp_path}/some_dir/different_user_nl_confusing"
        Path(notfound).touch()

        # Get and check results
        results = find_user_nl_files(str(tmp_path))
        assert len(results) == 2
        assert found_toplevel in results
        assert found_secondlevel in results
        assert notfound not in results


class TestCheckXmlFile:
    """Test _check_xml_file()"""

    @pytest.fixture
    def run_check(self, tmp_path):
        """
        Helper fixture that:
        - Creates a temp XML file (optionally malformed)
        - Patches ET.parse to raise a provided exception
        - Patches os.remove
        - Executes _check_xml_file
        - Returns (exception, mock_remove, stderr_output)
        """

        def _run(exception, *, malformed=False, capture_output=False):
            xml_file = tmp_path / "file.xml"
            if malformed:
                xml_file.write_text("<root>")  # malformed XML

            with (
                patch("ctsm.no_nans_in_inputs.namelist_utils.ET.parse", side_effect=exception),
                patch("ctsm.no_nans_in_inputs.namelist_utils.os.remove") as mock_remove,
            ):
                if capture_output:
                    # capsys must be accessed inside test, so we return control
                    with pytest.raises(type(exception)):
                        namelist_utils._check_xml_file(str(xml_file))
                else:
                    with pytest.raises(type(exception)):
                        namelist_utils._check_xml_file(str(xml_file))

            return str(xml_file), mock_remove

        return _run

    def test_valid_xml(self, tmp_path):
        """Test _check_xml_file() given a valid string to write"""
        xml_file = tmp_path / "valid.xml"
        xml_file.write_text("<root></root>")

        namelist_utils._check_xml_file(str(xml_file))
        assert xml_file.exists()

    @pytest.mark.parametrize(
        "exception, expected_stderr, malformed",
        [
            (IOError("boom"), "could not be opened", False),
            (ET.ParseError("bad xml"), "not well-formed", True),
            (RuntimeError("boom"), None, False),
        ],
        ids=["ioerror", "parse_error", "generic_exception"],
    )
    def test_exception_paths(
        self,
        run_check,
        capsys,
        exception,
        expected_stderr,
        malformed,
    ):  # pylint: disable=too-many-arguments,too-many-positional-arguments
        """Test _check_xml_file() given an invalid string to write"""
        xml_path, mock_remove = run_check(
            exception,
            malformed=malformed,
        )

        mock_remove.assert_called_once_with(xml_path)

        if expected_stderr is not None:
            captured = capsys.readouterr()
            assert expected_stderr in captured.err

            if malformed:
                assert "<root>" in captured.err
