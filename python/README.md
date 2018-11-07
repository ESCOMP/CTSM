# Testing the code here

## Unit tests

Unit tests can be run in one of two ways; these do the same thing, but
support different options:

1. via `make test`

   You can specify a few arguments to this:
   
   - python version: `make python=python3 test`
   - verbose: `make verbose=true test`
   - debug: `make debug=true test`
   
2. via `./run_ctsm_py_tests`

   You can specify various arguments to this; run `./run_ctsm_py_tests
   -h` for details
   
## pylint

You can run pylint on everything in the ctsm package with `make lint`
