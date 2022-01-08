"""General-purpose utility functions"""

import logging
import os
import sys
import string
import pdb
import subprocess

from datetime import date
from getpass import getuser
from configparser import NoSectionError, NoOptionError
from ctsm.path_utils import path_to_ctsm_root

logger = logging.getLogger(__name__)

# This string is used in the out-of-the-box ctsm.cfg and modify.cfg files
# to denote a value that needs to be filled in
_CONFIG_PLACEHOLDER = "FILL_THIS_IN"
# This string is used in the out-of-the-box ctsm.cfg and modify.cfg files
# to denote a value that can be filled in, but doesn't absolutely need to be
_CONFIG_UNSET = "UNSET"


def abort(errmsg):
    """Abort the program with the given error message

    No traceback is given, but if the logging level is DEBUG, then we'll enter pdb
    """
    if logger.isEnabledFor(logging.DEBUG):
        pdb.set_trace()

    sys.exit("ERROR: {}".format(errmsg))


def fill_template_file(path_to_template, path_to_final, substitutions):
    """Given a template file (based on python's template strings), write a copy of the
    file with template values filled in.

    Args:
    path_to_template (str): path to the existing template file
    path_to_final (str): path to where the final version will be written
    substitutions (dict): key-value pairs for the template string substitutions
    """

    with open(path_to_template) as template_file:
        template_file_contents = template_file.read()
    template = string.Template(template_file_contents)
    final_file_contents = template.substitute(substitutions)
    with open(path_to_final, "w") as final_file:
        final_file.write(final_file_contents)

def add_tag_to_filename(filename, tag):
    """
    Add a tag and replace timetag of a filename
    Expects file to end with [._]cYYMMDD.nc or [._]YYMMDD.nc
    Add the tag to just before that ending part
    and change the ending part to the current time tag.

    Parameters
    ----------
        filename (str) : file name
        tag (str) : string of a tag to be added to the end of filename

    Raises
    ------
        Error: When it cannot find . and _ in the filename.

    Returns
    ------
        fname_out (str): filename with the tag and date string added

    """
    basename = os.path.basename(filename)
    cend = -10
    if basename[cend] == "c":
        cend = cend - 1
    if (basename[cend] != ".") and (basename[cend] != "_"):
        logger.error("Trouble figuring out where to add tag to filename:" , filename)
        abort()
    today = date.today()
    today_string = today.strftime("%y%m%d")
    fname_out = basename[:cend] + "_" + tag + "_c" + today_string + ".nc"
    return fname_out


def update_metadata(file, title, summary, contact, data_script, description):
    """
    Description
    -----------
    Update netcdf file's metadata

    Arguments
    ---------
    title: No more than short one-sentence explanation.

    summary: No more than two-sentence explanation.

    contact: E.g. CAM bulletin board at https://bb.cgd.ucar.edu

    data_script: Script or instructions used to generate the dataset.

    description: Anything else that's relevant. Capturing the command-line
                 would be good (sys.argv) here or in data_script.
    """

    # update attributes
    today = date.today()
    today_string = today.strftime("%Y-%m-%d")

    # This is the required metadata for inputdata files
    file.attrs["title"] = title
    file.attrs["summary"] = summary
    file.attrs["creator"] = getuser()
    file.attrs["contact"] = contact
    file.attrs["creation_date"] = today_string
    file.attrs["data_script"] = data_script
    file.attrs["description"] = description

    # delete unrelated attributes if they exist
    del_attrs = [
        "source_code",
        "SVN_url",
        "hostname",
        "history",
        "History_Log",
        "Logname",
        "Host",
        "Version",
        "Compiler_Optimized",
    ]
    attr_list = file.attrs

    for attr in del_attrs:
        if attr in attr_list:
            del file.attrs[attr]


def lon_range_0_to_360(lon_in):
    """
    Description
    -----------
    Restrict longitude to 0 to 360 when given as -180 to 180.
    """
    if -180 <= lon_in < 0:
        lon_out = lon_in + 360
        logger.info(
            "Resetting longitude from %s to %s to keep in the range " " 0 to 360",
            str(lon_in),
            str(lon_out),
        )
    elif 0 <= lon_in <= 360 or lon_in is None:
        lon_out = lon_in
    else:
        errmsg = "lon_in needs to be in the range 0 to 360"
        abort(errmsg)

    return lon_out


def get_config_value(
    config,
    section,
    item,
    file_path,
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
    except NoSectionError:
        abort(
            "ERROR: Config file {} must contain section '{}'".format(file_path, section)
        )
    except NoOptionError:
        abort(
            "ERROR: Config file {} must contain item '{}' in section '{}'".format(
                file_path, item, section
            )
        )

    if val == _CONFIG_PLACEHOLDER:
        abort(
            "Error: {} needs to be specified in config file {}".format(item, file_path)
        )

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


def _handle_config_value(
    var, default, item, is_list, convert_to_type, can_be_unset, allowed_values
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


def _convert_to_bool(val):
    """Convert the given value to boolean

    Conversion is as in config files 'getboolean'
    """
    if val.lower() in ["1", "yes", "true", "on"]:
        return True
    if val.lower() in ["0", "no", "false", "off"]:
        return False
    raise ValueError("{} cannot be converted to boolean".format(val))
