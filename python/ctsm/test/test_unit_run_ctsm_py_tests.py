"""Unit tests for run_ctsm_py_tests module"""

import pytest
from ctsm.run_ctsm_py_tests import get_pytest_help_item


@pytest.fixture(name="sample_pytest_help")
def fixture_sample_pytest_help():
    """Sample pytest help text for testing"""
    return """usage: pytest [options] [file_or_dir] [file_or_dir] [...]

positional arguments:
  file_or_dir

general:
  -k EXPRESSION         Only run tests which match the given substring expression. An expression is a Python evaluable expression where all names are substring-matched against test names and their
                        parent classes. Example: -k 'test_method or test_other' matches all test functions and classes whose name contains 'test_method' or 'test_other', while -k 'not test_method'
                        matches those that don't contain 'test_method' in their names. -k 'not test_method and not test_other' will eliminate the matches. Additionally keywords are matched to classes
                        and functions containing extra names in their 'extra_keyword_matches' set, as well as functions which have names assigned directly to them. The matching is case-insensitive.
  -m MARKEXPR           Only run tests matching given mark expression. For example: -m 'mark1 and not mark2'.
  --markers             show markers (builtin, plugin and per-project ones).
  -x, --exitfirst       Exit instantly on first error or failed test
  --fixtures, --funcargs
                        Show available fixtures, sorted by plugin appearance (fixtures with leading '_' are only shown with '-v')
  --fixtures-per-test   Show fixtures per test
  --pdbcls=modulename:classname
                        Specify a custom interactive Python debugger for use with --pdb.For example: --pdbcls=IPython.terminal.debugger:TerminalPdb
  --pdb                 Start the interactive Python debugger on errors or KeyboardInterrupt
  --trace               Immediately break when running each test
  --capture=method      Per-test capturing method: one of fd|sys|no|tee-sys
  -s                    Shortcut for --capture=no
  --runxfail            Report the results of xfail tests as if they were not marked
  --lf, --last-failed   Rerun only the tests that failed at the last run (or all if none failed)
  --ff, --failed-first  Run all tests, but run the last failures first. This may re-order tests and thus lead to repeated fixture setup/teardown.
  --nf, --new-first     Run tests from new files first, then the rest of the tests sorted by file mtime
  --cache-show=[CACHESHOW]
                        Show cache contents, don't perform collection or tests. Optional argument: glob (default: '*').
  --cache-clear         Remove all cache contents at start of test run
  --lfnf, --last-failed-no-failures={all,none}
                        With ``--lf``, determines whether to execute tests when there are no previously (known) failures or when no cached ``lastfailed`` data was found. ``all`` (the default) runs the
                        full test suite again. ``none`` just emits a message about no known failures and exits successfully.
  --sw, --stepwise      Exit on test failure and continue from last failing test next time
  --sw-skip, --stepwise-skip
                        Ignore the first failing test but stop on the next failing test. Implicitly enables --stepwise.
"""


class TestGetPytestHelpItem:
    """Tests for get_pytest_help_item function"""

    def test_extract_short_option(self, sample_pytest_help):
        """Test extracting a short option with description"""
        result = get_pytest_help_item(sample_pytest_help, "-k")
        assert result.startswith("-k EXPRESSION")
        assert "Only run tests which match the given substring expression" in result
        assert "case-insensitive" in result

    def test_extract_long_option(self, sample_pytest_help):
        """Test extracting a long option with description"""
        result = get_pytest_help_item(sample_pytest_help, "--markers")
        assert result.startswith("--markers")
        assert "show markers (builtin, plugin and per-project ones)" in result

    def test_extract_option_with_short_and_long(self, sample_pytest_help):
        """Test extracting an option that has both short and long forms"""
        result = get_pytest_help_item(sample_pytest_help, "-x")
        assert result.startswith("-x, --exitfirst")
        assert "Exit instantly on first error or failed test" in result
        print(result)

    def test_extract_option_with_argument(self, sample_pytest_help):
        """Test extracting an option that takes an argument"""
        result = get_pytest_help_item(sample_pytest_help, "--pdbcls")
        assert result.startswith("--pdbcls=modulename:classname")
        assert "custom interactive Python debugger" in result

    def test_extract_multiline_description(self, sample_pytest_help):
        """Test that multiline descriptions are collapsed into single line"""
        result = get_pytest_help_item(sample_pytest_help, "--fixtures")
        assert "\n" not in result
        assert "Show available fixtures" in result
        assert "sorted by plugin appearance" in result

    def test_extract_option_with_optional_argument(self, sample_pytest_help):
        """Test extracting an option with optional argument"""
        result = get_pytest_help_item(sample_pytest_help, "--cache-show")
        assert result.startswith("--cache-show=[CACHESHOW]")
        assert "Show cache contents" in result
        assert "glob (default: '*')" in result

    def test_extract_option_with_choices(self, sample_pytest_help):
        """Test extracting an option with choices"""
        result = get_pytest_help_item(sample_pytest_help, "--lfnf")
        assert result.startswith("--lfnf, --last-failed-no-failures={all,none}")
        assert "With ``--lf``" in result

    def test_option_not_found_raises_error(self, sample_pytest_help):
        """Test that RuntimeError is raised when option is not found"""
        with pytest.raises(RuntimeError, match="Failed to get pytest help for --nonexistent"):
            get_pytest_help_item(sample_pytest_help, "--nonexistent")

    def test_only_matches_line_start(self, sample_pytest_help):
        """Test that option must be at start of line (after whitespace)"""
        # "pdb" appears in both "--pdb" and "--pdbcls", but searching for "pdb"
        # should not match "--pdbcls" since "pdb" is not at the start
        result = get_pytest_help_item(sample_pytest_help, "--pdb")
        assert result.startswith("--pdb")
        assert "Start the interactive Python debugger" in result
        # Should not contain pdbcls content
        assert "modulename:classname" not in result

    def test_handles_indented_options(self):
        """Test that function handles options with leading whitespace"""
        help_text = """
general:
  --foo                 This is foo option
                        with multiple lines
  --bar                 This is bar option
"""
        result = get_pytest_help_item(help_text, "--foo")
        assert result.startswith("--foo")
        assert "This is foo option" in result
        assert "with multiple lines" in result

    def test_stops_at_blank_line(self):
        """Test that description extraction stops at blank line"""
        help_text = """
  --option              First line
                        Second line

  --next                Next option
"""
        result = get_pytest_help_item(help_text, "--option")
        assert "First line" in result
        assert "Second line" in result
        assert "Next option" not in result

    def test_stops_at_less_indented_line(self):
        """Test that description extraction stops when indentation decreases"""
        help_text = """
  --option              First line
                        Second line
  --next                Next option
"""
        result = get_pytest_help_item(help_text, "--option")
        assert "First line" in result
        assert "Second line" in result
        assert "Next option" not in result

    def test_colon_separator_in_output(self, sample_pytest_help):
        """Test that output contains colon separator between header and description"""
        result = get_pytest_help_item(sample_pytest_help, "-s")
        assert ":" in result
        # The format includes the full header line, then colon, then description
        assert result.startswith("-s")
        assert "Shortcut for --capture=no" in result

    def test_empty_help_text(self):
        """Test behavior with empty help text"""
        with pytest.raises(RuntimeError, match="Failed to get pytest help for --option"):
            get_pytest_help_item("", "--option")

    def test_option_with_no_description(self):
        """Test option that has no description lines"""
        help_text = """
  --option

  --next                Next option
"""
        result = get_pytest_help_item(help_text, "--option")
        # Should return just the header with empty description
        assert result.startswith("--option:")
