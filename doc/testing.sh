#!/bin/bash
set -e
set -x

rm -rf _publish*

# Build all docs using container
echo "~~~~~ Build all docs using container"
# Also do a custom --conf-py-path
rm -rf _build _publish
d1="$PWD/_publish_container"
./build_docs_to_publish -r _build -d --site-root "$PWD/_publish" --conf-py-path doc-builder/test/conf.py --static-path ../_static --templates-path ../_templates
# VERSION LINKS WILL NOT RESOLVE IN _publish_container
cp -a _publish "${d1}"

# Build all docs using ctsm_pylib
echo "~~~~~ Build all docs using ctsm_pylib"
rm -rf _build _publish
d2="$PWD/_publish_nocontainer"
conda run -n ctsm_pylib ./build_docs_to_publish -r _build --site-root "$PWD/_publish"
# VERSION LINKS WILL NOT RESOLVE IN _publish_nocontainer
cp -a _publish "${d2}"

# Make sure container version is identical to no-container version
echo "~~~~~ Make sure container version is identical to no-container version"
diff -qr "${d1}" "${d2}"

# Check that -r -v works
echo "~~~~~ Check that -r -v works"
# Also do a custom --conf-py-path
rm -rf _build_docker
./build_docs -r _build_docker -v latest -d -c --conf-py-path doc-builder/test/conf.py --static-path ../_static --templates-path ../_templates

# Check that Makefile method works
echo "~~~~~ Check that Makefile method works"
rm -rf _build
make SPHINXOPTS="-W --keep-going" BUILDDIR=${PWD}/_build html

# Check that -b works
echo "~~~~~ Check that -b works"
rm -rf _build_docker
./build_docs -b _build_docker -d -c

# Check that doc-builder tests pass
echo "~~~~~ Check that doc-builder tests pass"
cd doc-builder/test
conda run -n ctsm_pylib make test

exit 0