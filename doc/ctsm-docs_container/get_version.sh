#!/usr/bin/env bash
set -e
cd doc/ctsm-docs_container

# Extract version from Dockerfile
version="$(grep "org.opencontainers.image.version" Dockerfile | cut -d'"' -f2 | sort |uniq)"
n_found=$(echo $version | wc -w)

# Error if anything other than exactly one version tag was found
if [[ ${n_found} -gt 1 ]]; then
    echo -e "Multiple version tags found:\n${version}" >&2
    exit 2
elif [[ ${n_found} -lt 1 ]]; then
    echo "Expected 1 but found 0 version tags" >&2
    exit -1
fi

# Error if version doesn't start with v
if [[ "${version}" != "v"* ]]; then
    echo "Version '${version}' doesn't start with v" >&2
    exit 22
fi

echo ${version}

exit 0
