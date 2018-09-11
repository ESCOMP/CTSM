#!/usr/bin/env bash

set -e
set -x

echo "building lilac"

# build lilac
mkdir -p /lilac/build
cd /lilac/build && cmake ..
make -j 4

echo "done building lilac, time to run the tests..."

# run test suite
ctest