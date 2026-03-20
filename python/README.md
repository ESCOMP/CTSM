# Testing the code here

## Running everything

To run all tests (unit tests, system tests and pylint), simply run `make
all` from this directory.

When you run `make all`, you need to first execute `module load nco`.

## Python environment

Another way is to use the file conda\_env\_ctsm\_py.txt to setup
a python environment. Comments in the file tell how to do this.

You can also use the script in the top level directory to do
all the conda commands and do this for you.

 ../py_env_create
 conda activate ctsm_pylib

Conda requirements files:

conda_env_ctsm_py.txt --------- Standard conda environment to use for most machines
conda_env_ctsm_py_latest.txt -- Test environment with latest versions that work

## Unit and system tests

Unit and system tests can be run in one of two ways; these do the same
thing, but support different options:

1. via `make test`

   You can specify a few arguments to this:
   
   - python version: `make python=python3.9 test` (defaults to `python3`; you should expect errors if trying to run with python2)
   - verbose: `make verbose=true test`
   - debug: `make debug=true test`

   Note that unit tests and system tests can be run separately with
   `make utest` or `make stest`, or they can all be run with `make
   test`.

   When you run `make test` or `make stest`, you need to first execute
   `module load nco`.

2. via `./run_ctsm_py_tests`

   You can specify various arguments to this; run `./run_ctsm_py_tests
   -h` for details. Please specify either --unit or --sys rather than
   not including any arguments.

   In any configuration where you run the system tests, you need to
   first execute `module load nco`.

## pylint

You can run pylint on everything in the ctsm package with `make lint`.

Note, that the listing of errors with pylint is very specific to the version
of pylint being used. Using a python environment as detailed earlier is important
in order to get a clean run of pylint.

## black

You can run a check for the black formatting with `make black`.
This won't change the code, but check if it would be reformatted
with black.
