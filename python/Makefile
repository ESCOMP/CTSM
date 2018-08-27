# Makefile for running tests on the python code here

# These variables can be overridden from the command-line
python = not-set
verbose = not-set
debug = not-set

ifneq ($(python), not-set)
PYTHON=$(python)
else
PYTHON=python
endif

# common args for running tests
TEST_ARGS=-m unittest discover

ifeq ($(debug), not-set)
    ifeq ($(verbose), not-set)
        # summary only output
        TEST_ARGS+=--buffer
    else
        # show individual test summary
        TEST_ARGS+=--buffer --verbose
    endif
else
    # show detailed test output
    TEST_ARGS+=--verbose
endif

PYLINT=pylint
PYLINT_ARGS=-j 2
PYLINT_SRC = \
	ctsm

.PHONY: test
test: FORCE
	$(PYTHON) $(TEST_ARGS)

.PHONY: lint
lint: FORCE
	$(PYLINT) $(PYLINT_ARGS) $(PYLINT_SRC)

.PHONY: clean
clean: FORCE
	find . -name '*.pyc' -exec rm {} \;

FORCE:

