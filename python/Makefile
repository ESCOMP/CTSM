# Makefile for running tests on the python code here

# These variables can be overridden from the command-line
python = not-set
verbose = not-set
debug = not-set

ifneq ($(python), not-set)
    PYTHON=$(python)
else
    PYTHON=python3
endif

ifneq ($(debug), not-set)
    TEST_ARGS+=--debug
endif
ifneq ($(verbose), not-set)
    TEST_ARGS+=--verbose
endif

PYLINT=pylint
PYLINT_ARGS=-j 4 --rcfile=ctsm/.pylintrc
PYLINT_SRC = \
	ctsm

all: test lint black
	@echo
	@echo
	@echo "Successfully ran all standard tests"

test: utest stest

.PHONY: utest
utest: FORCE
	$(PYTHON) ./run_ctsm_py_tests $(TEST_ARGS) --unit

.PHONY: stest
stest: FORCE
	$(PYTHON) ./run_ctsm_py_tests $(TEST_ARGS) --sys

.PHONY: lint
lint: FORCE
	$(PYLINT) $(PYLINT_ARGS) $(PYLINT_SRC)

.PHONY: black
# Run the black check on all of the python files here and undeneath.
# Use the black configure file to explicitly set a few things and specifiy the exact files.
black: FORCE
	black --check --config pyproject.toml .

.PHONY: run_black
# Run black on all of the python files here and undeneath.
# Use the black configure file to explicitly set a few things and specifiy the exact files.
run_black: FORCE
	black --config pyproject.toml .

.PHONY: clean
clean: FORCE
	find . -name '*.pyc' -exec rm {} \;

FORCE:

