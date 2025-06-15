"""
General-purpose utilities and functions for handling command-line
config files in ctsm python codes.
"""

import logging
import configparser

from ctsm.utils import abort
from ctsm.longitude import Longitude

logger = logging.getLogger(__name__)

# This string is used in the out-of-the-box ctsm.cfg and modify.cfg files
# to denote a value that needs to be filled in
_CONFIG_PLACEHOLDER = "FILL_THIS_IN"
# This string is used in the out-of-the-box ctsm.cfg and modify.cfg files
# to denote a value that can be filled in, but doesn't absolutely need to be
_CONFIG_UNSET = "UNSET"


def check_lon1_lt_lon2(lon1, lon2, lon_type):
    """
    Description
    -----------
    Given two longitudes, check that lon1 is < lon2. Useful for avoiding CTSM Issue #2017, but note
    that to use this function properly for that purpose, you need to have already converted
    longitudes from lon_type 180 to 360.
    """
    # Convert to Longitude class, if needed
    if not isinstance(lon1, Longitude):
        lon1 = Longitude(lon1, lon_type)
    if not isinstance(lon2, Longitude):
        lon2 = Longitude(lon2, lon_type)

    # Convert to type 360, if needed
    lon1 = lon1.get(360)
    lon2 = lon2.get(360)

    if lon1 < lon2:
        return

    msg = f"--lon1 ({lon1}) must be < --lon2 ({lon2})\n"
    msg += "See CTSM issue #2017: https://github.com/ESCOMP/CTSM/issues/2017"
    if lon_type == 180:
        msg = "After converting to --lon-type 360, " + msg
    raise ValueError(msg)


def get_config_value(
    config,
    section,
    item,
    file_path,
    *,
    allowed_values=None,
    default=None,
    is_list=False,
    convert_to_type=None,
    can_be_unset=False,
):
    """Get a given item from a given section of the config object
    Give a helpful error message if we can't find the given section or item
    Note that the file_path argument is only used for the sake of the error message
    If allowed_values is present, it should be a list of strings giving allowed values
    The function _handle_config_value determines what to do if we read:
    - a list or
    - a str that needs to be converted to int / float / bool
    - _CONFIG_UNSET: anything with the value "UNSET" will become "None"
    """
    try:
        val = config.get(section, item)
    except configparser.NoSectionError:
        abort("ERROR: Config file {} must contain section '{}'".format(file_path, section))
    except configparser.NoOptionError:
        abort(
            "ERROR: Config file {} must contain item '{}' in section '{}'".format(
                file_path, item, section
            )
        )

    if val == _CONFIG_PLACEHOLDER:
        abort("Error: {} needs to be specified in config file {}".format(item, file_path))

    val = _handle_config_value(
        var=val,
        default=default,
        item=item,
        is_list=is_list,
        convert_to_type=convert_to_type,
        can_be_unset=can_be_unset,
        allowed_values=allowed_values,
    )
    return val


def get_config_value_or_array(
    config,
    section,
    item,
    convert_to_type=None,
):
    """Get a config value as a single value or as an array if it's expressed as an array
    for cases when you don't know how it's going to be expressed"""
    val = config.get(section, item)
    vallist = val.split()
    if convert_to_type is not None:
        if (
            convert_to_type is not float
            and convert_to_type is not int
            and convert_to_type is not str
        ):
            abort(
                "get_config_value_or_array can only have convert_to_type as float, int or str not "
                + str(convert_to_type)
            )
    is_list = bool(len(vallist) > 1)

    val = _handle_config_value(
        var=val,
        default=None,
        item=item,
        is_list=is_list,
        convert_to_type=convert_to_type,
        can_be_unset=False,
        allowed_values=None,
    )
    return val


def _handle_config_value(
    *, var, default, item, is_list, convert_to_type, can_be_unset, allowed_values
):
    """
    Description
    -----------
    Assign the default value or the user-specified one to var.
    Convert from default type (str) to reqested type (int or float).

    If is_list is True, then default should be a list
    """
    if var == _CONFIG_UNSET:
        if can_be_unset:
            return default  # default may be None
        abort("Must set a value for .cfg file variable: {}".format(item))

    # convert string to list of strings; if there is just one element,
    # we will get a list of size one, which we will convert back to a
    # scalar later if needed
    var = var.split()

    if convert_to_type is bool:
        try:
            var = [_convert_to_bool(v) for v in var]
        except ValueError:
            abort("Non-boolean value found for .cfg file variable: {}".format(item))
    elif convert_to_type is not None:
        try:
            var = [convert_to_type(v) for v in var]
        except ValueError:
            abort("Wrong type for .cfg file variable: {}".format(item))

    if allowed_values is not None:
        for val in var:
            if val not in allowed_values:
                print("val = ", val, " in var not in allowed_values")
                errmsg = (
                    "{} is not an allowed value for {} in .cfg file. "
                    "Check allowed_values".format(val, item)
                )
                abort(errmsg)

    if not is_list:
        if len(var) > 1:
            abort("More than 1 element found for .cfg file variable: {}".format(item))
        var = var[0]

    return var


def _convert_to_bool(var):
    """
    Function for converting different forms of
    boolean strings to boolean value.

    Args:
        var (str): String bool input

    Raises:
        if the argument is not an acceptable boolean string
        (such as yes or no ; true or false ; y or n ; t or f ; 0 or 1).
        ValueError: The string should be one of the mentioned values.

    Returns:
        var_out (bool): Boolean value corresponding to the input.
    """
    if var.lower() in ("yes", "true", "t", "y", "1", "on"):
        var_out = True
    elif var.lower() in ("no", "false", "f", "n", "0", "off"):
        var_out = False
    else:
        raise ValueError("Boolean value expected. [true or false] or [y or n]")

    return var_out
