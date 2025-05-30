"""
Check that all XML files in our repo (except for those in submodules) are well-formed. Error if not.
"""
import sys
import glob
import os
import subprocess
import xml.etree.ElementTree as ET


def get_submodule_paths():
    """
    Get list of submodules
    """
    cmd = "git config --file .gitmodules --get-regexp path | awk '{ print $2 }'"
    result = subprocess.run(cmd, capture_output=True, text=True, shell=True, check=True)
    result_list = result.stdout.split("\n")[:-1]
    result_list = [x for x in result_list if x]
    return result_list


def is_in_submodule(file_path, submodule):
    """
    Return True if file is in given submodule, False otherwise
    """
    return os.path.commonpath([file_path, submodule]) == os.path.commonpath([submodule])


def is_in_any_submodule(file_path, submodule_paths):
    """
    Return True if file is in any submodule, False otherwise
    """
    file_path = os.path.abspath(file_path)
    submodule_paths = map(os.path.abspath, submodule_paths)
    return any(is_in_submodule(file_path, submodule) for submodule in submodule_paths)


def validate_xml(file_path):
    """
    Return True if XML file is well-formed, False otherwise
    """
    try:
        ET.parse(file_path)
    except ET.ParseError:
        print(f"❌ {file_path} is NOT well-formed")
        return False
    print(f"✅ {file_path} is well-formed")
    return True


def main():
    # pylint: disable=missing-function-docstring
    submodule_paths = get_submodule_paths()
    all_valid = True
    for xml_file in glob.glob("**/*.xml", recursive=True):
        if is_in_any_submodule(xml_file, submodule_paths):
            continue
        if not validate_xml(xml_file):
            all_valid = False

    if not all_valid:
        print("\nUse xmllint to show problems in malformed files")

    return all_valid


sys.exit(0 if main() else 1)
