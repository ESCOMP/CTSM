#!/usr/bin/env bash
set -e

# Process input arguments
suite_dir="$1"
if [[ ! -d "${suite_dir}" ]]; then
    echo "You must provide suite_dir" >&2
    exit 1
fi
shift
script="$1"
if [[ ! -f "${script}" ]]; then
    echo "You must provide script path" >&2
    exit 1
fi
shift
sdates_file="$1"
if [[ ! -f "${script}" ]]; then
    echo "You must provide sdates_file" >&2
    exit 1
fi

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
