# Updating `ctsm_pylib`

`ctsm_pylib` should be updated somewhat regularly to avoid using outdated Python versions. Here is a procedure for doing so.

## Update package versions
1. Edit `python/conda_env_ctsm_py.txt` to remove all version specifiers except for Python, which you should change to the version you want to use.
1. Run `conda` with that file.
1. Set `python/conda_env_ctsm_py.txt` to use the versions `conda` chose during its solve step.

## Resolve errors
1. Run `make lint` from `python/`. Resolve any errors (i.e., ignore anything that's just a warning or even less serious). This step will make subsequent steps easier, because `make lint` can catch errors very quickly with good error messages.
1. Run `black` with the version in the updated environment.
1. Run the Python unit tests (in `python/`, do `./run_ctsm_py_tests -u`). Fix any that are failing, then do the same with the Python system tests (`-s` instead of `-u`).
1. Run the `clm_pymods` test suite. Fix any failing tests.

## Update `pylint` settings
1. Generate a fresh `pylintrc` file for the new version of `pylint`: `pylint --generate-rcfile > filename`.
1. Look through the new `pylintrc` file, comparing settings to what we had in our existing one (`python/.pylintrc`). Change any default settings that our old one disagrees with, if you think it makes sense to do so.
1. Replace `python/.pylintrc` with your new file.
1. Run `make lint` from `python/`. Resolve all complaints. (`black` can help throughout this process.)

## Finish up
1. Repeat testing (`black`, `pylint`, Python system and unit tests, and `clm_pymods` suite) until all errors and complaints are resolved.
1. Update `.github/workflows/black.yml` to use the version of `black` in the updated `ctsm_pylib`.
