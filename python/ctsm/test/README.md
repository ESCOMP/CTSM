CTSM-specific boiler-plate needed for most unit test modules:

(1) If cime stuff is invoked by these unit tests (directly or
    indirectly, then: the first ctsm import statement near the top of
    the module should be:

from ctsm import add_cime_to_path # pylint: disable=unused-import

(2) Import the ctsm-specific unit_testing module:

from ctsm import unit_testing

(3) Allow names that pylint doesn't like:

# Allow names that pylint doesn't like, because otherwise I find it hard
# to make readable unit test names
# pylint: disable=invalid-name

(4) Have a 'main' block at the bottom:

if __name__ == '__main__':
    unit_testing.setup_for_tests()
    unittest.main()
