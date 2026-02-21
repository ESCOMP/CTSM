#!/usr/bin/env python3
"""
Unit tests for get_replacement_fill_values.py script.
"""

import numpy as np
import pytest
import xarray as xr

from ctsm.no_nans_in_inputs.constants import (
    USER_REQ_DELETE,
    USER_REQ_QUIT,
    USER_REQ_SKIP_FILE,
    USER_REQ_SKIP_VAR,
)
from ctsm.no_nans_in_inputs.get_replacement_fill_values import (
    FillValueConfig,
    VarContext,
    get_fill_value_from_user,
    var_data_has_nan,
)


# Test constants used in multiple tests
TEST_VAR_NAME = "test_var"
DEFAULT_VAR_CONTEXT = VarContext(var_name=TEST_VAR_NAME, target_type=float)
INT_VAR_CONTEXT = VarContext(var_name=TEST_VAR_NAME, target_type=int)


class TestGetFillValueFromUser:
    """Test the get_fill_value_from_user function."""

    def test_numeric_input(self, monkeypatch):
        """Test that a numeric input is converted to the target type."""
        expected = 3.14
        monkeypatch.setattr("builtins.input", lambda _: str(expected))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == expected
        assert isinstance(result, float)

    def test_integer_input(self, monkeypatch):
        """Test that an integer input is converted correctly."""
        expected = 42
        monkeypatch.setattr("builtins.input", lambda _: str(expected))
        result = get_fill_value_from_user(INT_VAR_CONTEXT, FillValueConfig())
        assert result == expected
        assert isinstance(result, int)

    def test_delete_command(self, monkeypatch):
        """Test that 'delete' returns USER_REQ_DELETE."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_DELETE)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig(allow_delete=True))
        assert result == USER_REQ_DELETE

    def test_delete_case_insensitive(self, monkeypatch):
        """Test that 'DELETE' is recognized case-insensitively."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_DELETE.upper())
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig(allow_delete=True))
        assert result == USER_REQ_DELETE

    def test_delete_not_allowed(self, monkeypatch):
        """Test that 'delete' is rejected when allow_delete is False, then accepts a number."""
        fallback_value = -1e20
        inputs = iter([USER_REQ_DELETE, str(fallback_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig(allow_delete=False))
        assert result == fallback_value

    def test_skip_variable(self, monkeypatch):
        """Test that 'skip' raises ValueError with SKIP_VARIABLE."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_SKIP_VAR)
        with pytest.raises(ValueError, match="SKIP_VARIABLE"):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_skip_file(self, monkeypatch):
        """Test that 'skipfile' raises ValueError with SKIP_FILE."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_SKIP_FILE)
        with pytest.raises(ValueError, match="SKIP_FILE"):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_quit(self, monkeypatch):
        """Test that 'quit' raises KeyboardInterrupt."""
        monkeypatch.setattr("builtins.input", lambda _: USER_REQ_QUIT)
        with pytest.raises(KeyboardInterrupt):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_default_value_on_empty_input(self, monkeypatch):
        """Test that empty input uses the default value."""
        default = 1.5e36
        monkeypatch.setattr("builtins.input", lambda _: "")
        result = get_fill_value_from_user(
            DEFAULT_VAR_CONTEXT, FillValueConfig(default_value=default)
        )
        assert result == default

    def test_default_delete_on_empty_input(self, monkeypatch):
        """Test that empty input uses delete as default."""
        monkeypatch.setattr("builtins.input", lambda _: "")
        config = FillValueConfig(default_value=USER_REQ_DELETE)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, config)
        assert result == USER_REQ_DELETE

    def test_auto_delete_when_default_is_delete(self):
        """Test that delete_if_none_filled auto-deletes without prompting."""
        config = FillValueConfig(default_value=USER_REQ_DELETE, delete_if_none_filled=True)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, config)
        assert result == USER_REQ_DELETE

    def test_auto_delete_only_when_default_is_delete(self, monkeypatch):
        """Test that delete_if_none_filled does NOT auto-delete when default is not delete."""
        entered_value = 0.0
        default = -888.0
        monkeypatch.setattr("builtins.input", lambda _: str(entered_value))
        config = FillValueConfig(default_value=default, delete_if_none_filled=True)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, config)
        assert result == entered_value

    def test_nan_input_rejected(self, monkeypatch):
        """Test that NaN input is rejected, then accepts a valid number."""
        valid_value = 42.5
        inputs = iter(["nan", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_invalid_input_then_valid(self, monkeypatch):
        """Test that invalid input is rejected, then accepts a valid number."""
        valid_value = -777.0
        inputs = iter(["not_a_number", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_empty_input_no_default_then_valid(self, monkeypatch):
        """Test that empty input with no default shows help, then accepts a value."""
        valid_value = 100.0
        inputs = iter(["", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_whitespace_input_treated_as_empty(self, monkeypatch):
        """Test that whitespace-only input is treated as empty."""
        valid_value = -555.0
        inputs = iter(["   ", str(valid_value)])
        monkeypatch.setattr("builtins.input", lambda _: next(inputs))
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_ctrl_c_twice_raises(self, monkeypatch):
        """Test that pressing Ctrl-C twice raises KeyboardInterrupt."""

        def raise_keyboard_interrupt(_):
            raise KeyboardInterrupt

        monkeypatch.setattr("builtins.input", raise_keyboard_interrupt)
        with pytest.raises(KeyboardInterrupt):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())

    def test_ctrl_c_then_valid_input(self, monkeypatch):
        """Test that Ctrl-C followed by valid input works."""
        valid_value = 123.0
        # call_count is defined in the outer (test method) scope
        call_count = 0

        def input_with_ctrl_c(_):
            # nonlocal lets us modify call_count from the enclosing scope
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                # First call: simulate Ctrl-C. get_fill_value_from_user will
                # catch this, show ncdump output, and loop back to prompt again.
                raise KeyboardInterrupt
            # Second call: return valid input, which should be accepted
            return str(valid_value)

        monkeypatch.setattr("builtins.input", input_with_ctrl_c)
        result = get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())
        assert result == valid_value

    def test_ctrl_c_then_quit(self, monkeypatch):
        """Test that Ctrl-C followed by 'quit' raises KeyboardInterrupt."""
        call_count = 0

        def input_with_ctrl_c_then_quit(_):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise KeyboardInterrupt
            return USER_REQ_QUIT

        monkeypatch.setattr("builtins.input", input_with_ctrl_c_then_quit)
        with pytest.raises(KeyboardInterrupt):
            get_fill_value_from_user(DEFAULT_VAR_CONTEXT, FillValueConfig())


class TestVarDataHasNan:
    """Test the var_data_has_nan function."""

    def test_true_when_nan_present(self):
        """Test detection of NaN values."""
        da = xr.DataArray([[1.0, 2.0], [np.nan, 4.0]])
        assert var_data_has_nan(da)

    def test_false_when_no_nan(self):
        """Test correct response for arrays without NaN values."""
        da = xr.DataArray([1.0, 2.0, 3.0, 4.0])
        assert not var_data_has_nan(da)

    def test_true_when_all_nan(self):
        """Test that arrays with all NaN values are detected."""
        da = xr.DataArray([np.nan, np.nan, np.nan])
        assert var_data_has_nan(da)

    def test_false_for_empty_array(self):
        """Test that empty arrays return False (no NaNs)."""
        da = xr.DataArray([])
        assert not var_data_has_nan(da)
