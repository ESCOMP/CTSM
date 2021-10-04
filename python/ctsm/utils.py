"""General-purpose utility functions"""

import logging
import os
import sys
import string
import pdb
import subprocess

from datetime import date
from getpass import getuser

logger = logging.getLogger(__name__)

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
    lon_1 = -180
    lon_2 = 0
    lon_3 = 180
    lon_4 = 360
    if lon_1 <= lon_in < lon_2:
        lon_out = lon_in + lon_4
        message = 'WARNING: Resetting longitude from ' + str(lon_in) + \
                  ' to ' + str(lon_out) + ' to keep in the range 0 to 360'
        print(message)  # TODO Use logging to print this
    elif lon_2 <= lon_in <= lon_4 or lon_in == -999:
        lon_out = lon_in
    else:
        errmsg = 'lon_in needs to be in the range 0 to 360'
        abort(errmsg)

    return lon_out
