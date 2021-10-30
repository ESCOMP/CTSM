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

logger = logging.getLogger(__name__)

# This string is used in the out-of-the-box ctsm.cfg and modify.cfg files
# to denote a value that needs to be filled in
_CONFIG_PLACEHOLDER = 'FILL_THIS_IN'
# This string is used in the out-of-the-box ctsm.cfg and modify.cfg files
# to denote a value that can be filled in, but doesn't absolutely need to be
CONFIG_UNSET = 'UNSET'

def abort(errmsg):
    """Abort the program with the given error message

    No traceback is given, but if the logging level is DEBUG, then we'll enter pdb
    """
    if logger.isEnabledFor(logging.DEBUG):
        pdb.set_trace()

    sys.exit('ERROR: {}'.format(errmsg))

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
    with open(path_to_final, 'w') as final_file:
        final_file.write(final_file_contents)

def get_git_sha():
    """
    Returns Git short SHA for the currect directory.
    """
    return subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).strip().decode()

def add_tag_to_filename(filename, tag):
    """
    Add a tag and replace timetag of a filename
    Expects file to end with [._]cYYMMDD.nc or [._]YYMMDD.nc
    Add the tag to just before that ending part
    and change the ending part to the current time tag
    """

    basename = os.path.basename(filename)
    cend = -10

    if basename[cend] == "c":
        cend = cend - 1
    if ( (basename[cend] != ".") and (basename[cend] != "_") ):
        errmsg = 'Trouble figuring out where to add tag to filename: ' + filename
        abort(errmsg)

    today = date.today()
    today_string = today.strftime("%y%m%d")

    return basename[:cend] + "_" + tag + "_c" + today_string + '.nc'

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

    #update attributes
    today = date.today()
    today_string = today.strftime("%Y-%m-%d")

    # This is the required metadata for inputdata files
    file.attrs['title'] = title
    file.attrs['summary'] = summary
    file.attrs['creator'] = getuser()
    file.attrs['contact'] = contact
    file.attrs['creation_date'] = today_string
    file.attrs['data_script'] = data_script
    file.attrs['description'] = description

    #delete unrelated attributes if they exist
    del_attrs = ['source_code', 'SVN_url', 'hostname', 'history'
                 'History_Log', 'Logname', 'Host', 'Version',
                 'Compiler_Optimized']
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
        message = 'INFO: Resetting longitude from ' + str(lon_in) + \
                  ' to ' + str(lon_out) + ' to keep in the range 0 to 360'
        print(message)  # TODO Use logging to print this
    elif 0 <= lon_in <= 360 or lon_in is None:
        lon_out = lon_in
    else:
        errmsg = 'lon_in needs to be in the range 0 to 360'
        abort(errmsg)

    return lon_out

def get_config_value(config, section, item, file_path, allowed_values=None, default=None, is_list=False, convert_to_type=None):
    """Get a given item from a given section of the config object
    Give a helpful error message if we can't find the given section or item
    Note that the file_path argument is only used for the sake of the error message
    If allowed_values is present, it should be a list of strings giving allowed values
    """
    try:
        val = config.get(section, item)
    except NoSectionError:
        abort("ERROR: Config file {} must contain section '{}'".format(file_path, section))
    except NoOptionError:
        abort("ERROR: Config file {} must contain item '{}' in section '{}'".format(
            file_path, item, section))

    if val == _CONFIG_PLACEHOLDER:
        abort("Error: {} needs to be specified in config file {}".format(item, file_path))

    if allowed_values is not None:
        if val not in allowed_values:
            abort("Error: {} is not an allowed value for {} in config file {}\n"
                  "Allowed values: {}".format(val, item, file_path, allowed_values))

    val = _handle_config_value(var=val, default=default,
                               is_list=is_list,
                               convert_to_type=convert_to_type)

    return val

def _handle_config_value(var, default, is_list, convert_to_type):
    """
    Description
    -----------
    Assign the default value or the user-specified one to var.
    Convert from default type (str) to reqested type (int or float).
    """
    if var == CONFIG_UNSET:
        var = default  # default may be None
    elif is_list:
        var = list(var.split())  # convert string to list of strings
        var = list(map(convert_to_type, var))  # convert all elements
    elif convert_to_type is not None:
        var = convert_to_type(var)
    else:
        var = var

    return var
