"""General-purpose utility functions"""

import logging
import os
import sys
import string
import re
import pdb

from datetime import date
from getpass import getuser

from ctsm.git_utils import get_ctsm_git_short_hash

logger = logging.getLogger(__name__)


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


def add_tag_to_filename(filename, tag, replace_res=False):
    """
    Add a tag and replace timetag of a filename
    Expects file to end with [._]cYYMMDD.nc or [._]YYMMDD.nc
    or with 4-digit years YYYYMMDD.
    Add the tag to just before that ending part
    and change the ending part to the current time tag.

    if replace_res is True, then replace the resolution
    part of the filename. Expects the file to start with
    [a-z.]_ and then the resolution.

    Parameters
    ----------
        filename (str) : file name
        tag (str) : string of a tag to be added to the end of filename
                    (or to replace the resolution part of the filename)

    Raises
    ------
        Error: When it cannot find . and _ in the filename.
        Error: When it's asked to replace the resolution and
               can't figure out where that is in the filename.

    Returns
    ------
        fname_out (str): filename with the tag and date string added

    """
    basename = os.path.basename(filename)
    cend = -10
    if basename[cend] == "c":
        cend = cend - 1
    if (basename[cend] != ".") and (basename[cend] != "_"):
        # Check if date stirng at end includes a 4 digit year
        cend = -12
        if basename[cend] == "c":
            cend = cend - 1
        if (basename[cend] != ".") and (basename[cend] != "_"):
            err_msg = "Trouble figuring out where to add tag to filename: " + filename
            abort(err_msg)
    today = date.today()
    today_string = today.strftime("%y%m%d")
    if not replace_res:
        fname_out = basename[:cend] + "_" + tag + "_c" + today_string + ".nc"
    else:
        match = re.fullmatch(r"([a-z.]+)_([Cfvnenp0-9x.crunldasA-Z]+)_(.+?)", basename[:cend])
        if match is not None:
            fname_out = (
                match.group(1) + "_" + tag + "_" + match.group(3) + "_c" + today_string + ".nc"
            )
        else:
            abort(
                "Trouble figuring out where to replace the resolution in the filename: " + filename
            )
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


def write_output(file, file_in, file_out, file_type):
    """
    Description
    -----------
    Write output file

    Arguments
    ---------
    file_in:
        (str) User-defined entry of input file
    file_out:
        (str) User-defined entry of output file
    file_type:
        (str) examples: mesh, fsurdat
    """

    # update attributes
    title = "Modified " + file_type + " file"
    summary = "Modified " + file_type + " file"
    contact = "N/A"
    data_script = os.path.abspath(__file__) + " -- " + get_ctsm_git_short_hash()
    description = "Modified this file: " + file_in
    update_metadata(
        file,
        title=title,
        summary=summary,
        contact=contact,
        data_script=data_script,
        description=description,
    )

    # mode 'w' overwrites file if it exists
    file.to_netcdf(path=file_out, mode="w", format="NETCDF3_64BIT")
    logger.info("Successfully created: %s", file_out)
    file.close()
