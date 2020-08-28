# Testing the code here

## Running everything

To run all tests (unit tests, system tests and pylint), simply run `make
all` from this directory.

## Unit and system tests

Unit and system tests can be run in one of two ways; these do the same
thing, but support different options:

1. via `make test`

   You can specify a few arguments to this:
   
   - python version: `make python=python3 test`
   - verbose: `make verbose=true test`
   - debug: `make debug=true test`

   Note that unit tests and system tests can be run separately with
   `make utest` or `make stest`, or they can all be run with `make
   test`.

2. via `./run_ctsm_py_tests`

   You can specify various arguments to this; run `./run_ctsm_py_tests
   -h` for details

## pylint

You can run pylint on everything in the ctsm package with `make lint`.

Note: you should expect some errors if using a python2 version of
pylint, but this should be clean if running with a python3 version.
