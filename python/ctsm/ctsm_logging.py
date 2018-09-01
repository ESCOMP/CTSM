"""Utilities to facilitate logging

A guide to logging in ctsm python scripts:

- At the top of each module, you should have:
  logger = logging.getLogger(__name__)

- Logging should be done via that logger, NOT via logging.[whatever]

- If you want to allow the user to control logging via command-line arguments, you should:

  (1) At the very start of a script / application, call setup_logging_pre_config(). (We
      need to initialize logging to avoid errors from logging calls made very early in the
      script.)

  (2) When setting up the argument parser, call add_logging_args(parser)

  (3) After parsing arguments, call process_logging_args(args)

- If you don't want to allow the user to control logging via command-line arguments, then
  simply:

  (1) At the very start of a script / application, call setup_logging() with the desired
      arguments

- In unit tests, to avoid messages about loggers not being set up, you should call
  setup_logging_for_tests (this is typically done via unit_testing.setup_for_tests)
"""

import logging

logger = logging.getLogger(__name__)

def setup_logging_pre_config():
    """Setup logging for a script / application

    This function should be called at the very start of a script / application where you
    intend to allow the user to control logging preferences via command-line arguments.

    This sets initial options that may be changed later by process_logging_args.
    """
    setup_logging(level=logging.WARNING)

def setup_logging_for_tests(enable_critical=False):
    """Setup logging as appropriate for unit tests"""
    setup_logging(level=logging.CRITICAL)
    if not enable_critical:
        logging.disable(logging.CRITICAL)

def setup_logging(level=logging.WARNING):
    """Setup logging for a script / application

    This function should be called at the very start of a script / application where you
    do NOT intend to allow the user to control logging preferences via command-line
    arguments, so that all of the final logging options are set here.
    """
    logging.basicConfig(format='%(levelname)s: %(message)s', level=level)

def add_logging_args(parser):
    """Add common logging-related options to the argument parser"""

    logging_level = parser.add_mutually_exclusive_group()

    logging_level.add_argument('-v', '--verbose', action='store_true',
                               help='Output extra logging info')

    logging_level.add_argument('--debug', action='store_true',
                               help='Output even more logging info for debugging')

def process_logging_args(args):
    """Configure logging based on the logging-related args added by add_logging_args"""
    root_logger = logging.getLogger()

    if args.debug:
        root_logger.setLevel(logging.DEBUG)
    elif args.verbose:
        root_logger.setLevel(logging.INFO)
    else:
        root_logger.setLevel(logging.WARNING)
