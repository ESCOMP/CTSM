#!/bin/bash
set -e

script="./cs.status"
if [[ ! -e "${script}" ]]; then
    command="ls ./cs.status.[0-9]*"
    printed=0
    while ! ${command} 1>/dev/null 2>&1; do
        if compgen -G "STDERR*" > /dev/null; then
            if [[ $(stat -c %s STDERR* | sort | tail -n 1) -gt 0 ]]; then
                echo "STDERR file(s) found; stopping." >&2
                exit 1
            fi
        elif [[ ${printed} -eq 0 ]]; then
            echo "Waiting for ./cs.status.[0-9]* file(s) to appear..." >&2
            printed=1
        fi
        sleep 1
    done
    script="$(${command})"
fi

echo ${script}
exit 0
