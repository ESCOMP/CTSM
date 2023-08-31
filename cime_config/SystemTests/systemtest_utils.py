"""
Reduce code duplication by putting reused functions here.
"""

import os, subprocess


def cmds_to_setup_conda(caseroot):
    # Add specific commands needed on different machines to get conda available
    # Use semicolon here since it's OK to fail
    #
    conda_setup_commands = ". " + caseroot + "/.env_mach_specific.sh; "
    # Execute the module unload/load when "which conda" fails
    # eg on cheyenne
    try:
        subprocess.run("which conda", shell=True, check=True)
    except subprocess.CalledProcessError:
        # Remove python and add conda to environment for cheyennne
        conda_setup_commands += " module unload python; module load conda;"

    return conda_setup_commands


def run_python_script(caseroot, this_conda_env, command, tool_path):

    # Run in the specified conda environment
    conda_setup_commands = cmds_to_setup_conda(caseroot)
    conda_setup_commands += f" conda run -n {this_conda_env}"

    # Finish with Python script call
    command = conda_setup_commands + " " + command
    print(f"command: {command}")

    # Run with logfile
    tool_name = os.path.split(tool_path)[-1]
    try:
        with open(tool_name + ".log", "w") as f:
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
