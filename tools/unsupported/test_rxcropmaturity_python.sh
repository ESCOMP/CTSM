#!/usr/bin/env bash
set -e

# This script is designed to easily test updated versions of the code in python/ctsm/crop_calendars/
# on outputs of the rxcropmaturity test suite. For RXCROPMATURITYSKIPGEN tests (i.e., skipping the
# GDD-generating step, it will rerun check_rxboth_run.py. It will also do this for RXCROPMATURITY
# tests whose call of generate_gdds.py completed successfully. For RXCROPMATURITY tests where that
# failed, it will retry the generate_gdds.py call. Tests are performed with both the ctsm_pylib
# and npl conda environments.
#
# Note that the python/ctsm/crop_calendars/ this test uses will be in the same CTSM directory as
# you used to start the tests.
#
# The script takes one positional input:
#     suite_dir: The directory where your test suite was performed. E.g., $SCRATCH/tests_0123-142858de.
#
# Output will look something like this (✅ for success, 🔴 for failure):
#     ctsm_pylib
#        RXCROPMATURITYINST_Lm61.f10_f10_mg37.IHist check_rxboth_run ✅
#        RXCROPMATURITY_Lm61.f10_f10_mg37.IHist generate_gdds 🔴
#        RXCROPMATURITYSKIPGENINST_Ld1097.f10_f10_mg37.IHist check_rxboth_run ✅
#        RXCROPMATURITYSKIPGEN_Ld1097.f10_f10_mg37.IHist check_rxboth_run 🔴
#     npl
#        RXCROPMATURITYINST_Lm61.f10_f10_mg37.IHist check_rxboth_run ✅
#        RXCROPMATURITY_Lm61.f10_f10_mg37.IHist generate_gdds ✅
#        RXCROPMATURITYSKIPGENINST_Ld1097.f10_f10_mg37.IHist check_rxboth_run ✅
#        RXCROPMATURITYSKIPGEN_Ld1097.f10_f10_mg37.IHist check_rxboth_run ✅
#
# Log files for each will be saved as TEST_SHORTNAME.CONDA_ENV.log.

# Process input arguments
suite_dir="$1"
if [[ ! -d "${suite_dir}" ]]; then
    echo "You must provide suite_dir" >&2
    exit 1
fi

# Where do we save log files?
log_dir="$PWD"

cd "${suite_dir}"
test_dir_list="$(ls -d RXCROPMATURITY*/ | grep -v gddgen)"

for conda_env in ctsm_pylib npl; do
    echo "${conda_env}"
    for d in ${test_dir_list}; do

        # Get test shortname
        t="$(echo $d | grep -oE ".*IHist")"
        echo -n "   $t "

        # Set up
        pushd $d 1>/dev/null 2>&1
        logfile="${log_dir}/${t}.${conda_env}.log"

        # Get python command to test
        script="check_rxboth_run"
        set +e
        cmd="$(grep -h "${script}.py" *int.o* | grep -oE "python3.*" | tail -n 1)"
        set -e
        if [[ "${cmd}" == "" ]]; then
            # check_rxboth_run.py wasn't run. Look for generate_gdds.py command.
            script="generate_gdds"
            set +e
            cmd="$(grep -h "${script}.py" *int.o* | grep -oE "python3.*" | tail -n 1)"
            set -e

            # Neither were found
            if [[ "${cmd}" == "" ]]; then
                echo -e "\n      Command not found" >&2
                popd 1>/dev/null 2>&1
                continue
            fi

            cd ../$t*gddgen
        fi

        # Strip extraneous text off command
        cmd="$(echo ${cmd} | sed "s/ returned non-zero exit status 1.//")"
        cmd="$(echo ${cmd} | sed -E "s/ '$//")"
        cmd="$(echo ${cmd} | sed -E "s/'$//")"
        
        # generate_gdds.py should be run with --no-pickle for troubleshooting
        if [[ "${cmd}" == *"generate_gdds.py "* ]]; then
            cmd="${cmd} --no-pickle"
        fi

        # Run the command
        echo -n "${script} "
        set +e
        conda run -n ${conda_env} ${cmd} 1>"${logfile}" 2>&1
        result=$?
        set -e

        # Print emoji to indicate result
        if [[ ${result} -eq 0 ]]; then
            echo ✅
        else
            echo 🔴
        fi

        # Return to test suite directory
        popd 1>/dev/null 2>&1
    done
done

exit 0
