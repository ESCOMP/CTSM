#!/usr/bin/env bash

# Fail on any non-zero exit code
set -e

# Compare docs built with container vs. ctsm_pylib

cli_tool="$1"

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
cd "${SCRIPT_DIR}/.."

rm -rf _publish*

# Build all docs using container
echo "~~~~~ Build all docs using container"
# Also do a custom --conf-py-path
rm -rf _build _publish
d1="$PWD/_publish_container"
./build_docs_to_publish -r _build -d --site-root "$PWD/_publish"
# VERSION LINKS WILL NOT RESOLVE IN _publish_container
cp -a _publish "${d1}"

# Build all docs using ctsm_pylib
echo "~~~~~ Build all docs using ctsm_pylib"
rm -rf _build _publish
d2="$PWD/_publish_nocontainer"
conda run -n ctsm_pylib --no-capture-output ./build_docs_to_publish -r _build --site-root "$PWD/_publish"  --conf-py-path doc-builder/test/conf.py --static-path ../_static --templates-path ../_templates
# VERSION LINKS WILL NOT RESOLVE IN _publish_nocontainer
cp -a _publish "${d2}"

# Make sure container version is identical to no-container version
echo "~~~~~ Make sure container version is identical to no-container version"
diff -qr "${d1}" "${d2}"
echo "Successful: Docs built with container are identical to those built without"

exit 0
