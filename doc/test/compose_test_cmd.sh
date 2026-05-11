# This should only be run locally within another shell

if [[ "${cli_tool}" == "" ]]; then
    echo "${msg} (no container)"
else
    cmd="${cmd} -d"
    if [[ "${cli_tool}" != "default" ]]; then
        cmd="${cmd} --container-cli-tool ${cli_tool}"
    fi
    echo "${msg} (container: ${cli_tool})"
fi

echo cmd
