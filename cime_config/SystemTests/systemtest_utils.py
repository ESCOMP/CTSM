"""
Reduce code duplication by putting reused functions here.
"""

import os, subprocess, re, glob
from collections import OrderedDict


def cmds_to_setup_conda(caseroot):
    # Add specific commands needed on different machines to get conda available
    # Use semicolon here since it's OK to fail
    #
    conda_setup_commands = ". " + caseroot + "/.env_mach_specific.sh; "
    # Setting CONDA_PREFIX to empty ensures that this works even if called from
    # a shell with a conda environment activated
    conda_setup_commands += "CONDA_PREFIX=; "
    # Execute the module unload/load when "which conda" fails
    # eg on cheyenne
    try:
        subprocess.run("which conda", shell=True, check=True)
    except subprocess.CalledProcessError:
        # Remove python and add conda to environment for cheyennne
        unload_python_load_conda = "module unload python; module load conda;"
        # Make sure that adding this actually loads conda
        subprocess.run(unload_python_load_conda + "which conda", shell=True, check=True)
        # Save
        conda_setup_commands += " " + unload_python_load_conda

    return conda_setup_commands


def cmds_to_run_via_conda(caseroot, conda_run_call, command):
    # Run in the specified conda environment
    conda_setup_commands = cmds_to_setup_conda(caseroot)
    conda_setup_commands += " " + conda_run_call

    # Finish with Python script call
    command = conda_setup_commands + " " + command
    print(f"command: {command}")

    return command


def run_python_script(caseroot, this_conda_env, command_in, tool_path):

    # First, try with "conda run -n"
    command = cmds_to_run_via_conda(caseroot, f"conda run -n {this_conda_env}", command_in)

    # Run with logfile
    tool_name = os.path.split(tool_path)[-1]
    try:
        with open(tool_name + ".log", "w") as f:
            subprocess.run(
                command, shell=True, check=True, text=True, stdout=f, stderr=subprocess.STDOUT
            )
    except subprocess.CalledProcessError as error:
        # Retry with the original "conda activate" method
        command = cmds_to_run_via_conda(
            caseroot,
            f"conda activate {this_conda_env} && ",
            command_in,
        )
        try:
            with open(tool_name + ".log2", "w") as f:
                subprocess.run(
                    command, shell=True, check=True, text=True, stdout=f, stderr=subprocess.STDOUT
                )
        except subprocess.CalledProcessError as error:
            print("ERROR while getting the conda environment and/or ")
            print(f"running the {tool_name} tool: ")
            print(f"(1) If your {this_conda_env} environment is out of date or you ")
            print(f"have not created the {this_conda_env} environment, yet, you may ")
            print("get past this error by running ./py_env_create ")
            print("in your ctsm directory and trying this test again. ")
            print("(2) If conda is not available, install and load conda, ")
            print("run ./py_env_create, and then try this test again. ")
            print("(3) If (1) and (2) are not the issue, then you may be ")
            print(f"getting an error within {tool_name} itself. ")
            print("Default error message: ")
            print(error.output)
            raise
        except:
            print(f"ERROR trying to run {tool_name}.")
            raise
    except:
        print(f"ERROR trying to run {tool_name}.")
        raise


# Read a user_nl file and return the namelist option if found
def find_user_nl_option(caseroot, component, namelist_option):

    # This is a copy of the CIME _get_list_of_user_nl_files
    # which could be used if this moved into the CIME project
    file_pattern = "user_nl_" + component + "*"
    file_list = glob.glob(os.path.join(caseroot, file_pattern))

    # Check that there is at least one file
    if len(file_list) == 0:
        raise RuntimeError("No user_nl files found for component " + component)

    # Read through the file list and look for a match and return the whole entry
    output = OrderedDict()
    for one_file in file_list:
        with open(one_file, "r") as user_nl_file:
            user_nl_text = user_nl_file.read()
            reg = rf"{namelist_option}.*?(?=,|\n)"
            find_out = re.findall(reg, user_nl_text)
            output[one_file] = find_out
    return output
