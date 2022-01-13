.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This documents various conventions for the user interface of python scripts. The focus here is on user-visible aspects, as opposed to implementation details.

====================================================
 Separation of front-end script from implementation
====================================================

Most python code resides in the ``python/ctsm`` directory, and modules there have a ``.py`` extension. However, scripts that are run from the command-line have small wrapper files that reside in appropriate places throughout the repository (e.g., in the ``tools`` directory). These wrapper files should *not* have a ``.py`` extension, and should contain as little python code as possible (since these files aren't checked by pylint and cannot easily be unit tested): typically they contain just enough code to setup the python path, load a python module and then call the main function from that module. See examples throughout CTSM for details.

Rationale: Modules meant to be imported (i.e., everything under ``python/ctsm``) should have a ``.py`` extension. However, it is valuable to keep the extension off of the scripts that users run for a few reasons:
1. A user shouldn't need to know what language a script is in when all they want to do is to run the script
2. We want to avoid the need for retraining users, updating documentation, etc. if we change the implementation language
3. Standard Unix utilities rarely if ever require specifying the language extension when running them; the same is true for the core scripts in CIME (``create_newcase``, ``xmlchange``, etc.); we would like to stay consistent with these other scripts and utilities

Counter-arguments: Arguments for keeping a ``.py`` extension even on the files meant to be run by users are the following:
1. It's obvious from looking at a directory listing what language a file is in
2. This extension may be needed on Windows systems
3. The file extension is needed to support linting, unit testing, etc.

Although there are good arguments on both sides, the stability of the user experience is supported by excluding the language extension. We choose to prioritize the user experience over developer convenience. (The possible inability of running on Windows systems does not currently feel like an issue, because users don't typically try to run these scripts on Windows systems as far as we know.)

========================================
 Conventions for command-line arguments
========================================

1. Longer options:

Options that are more than a single character should be formatted as ``--some-variable`` -- i.e., with a double leading hyphen and a hyphen as word separator.

   * Rationale: This follows the `POSIX conventions <https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html>`_.

2. Standard options: verbose, silent, help and debug:

Scripts should support ``--verbose`` / ``-v`` and ``--debug`` options. See comments at the top of ``ctsm_logging.py`` for details.
Also the options silent and help are recommended as well.

3. Value flags:

For arguments that need a value: If there is no obvious default for the value, then the argument should be required rather than being given a default that the user will often need to override. A general rule of thumb is: if, say, 75% -- 90% of users will want a particular value, then it can be okay to have this value as the default; otherwise, it should be a required argument. This holds even for arguments that are only needed in certain circumstances: For example, if argument ``--fooval FOOVAL`` needs to be given whenever ``--bar`` is given, but there is no obvious default for ``FOOVAL``, then the argument parsing should be coded so that ``--fooval`` does not have a default, and there is error checking code that ensures that ``fooval`` has been specified (not ``None``) if ``--bar`` was given.

   * Rationale: Giving an argument a default value means that we can't catch and report when a user forgets to set it but should have set it. Leaving off the default allows for this error checking.
   * Example: To set the start year for DATM forcing: ``--datm-syr 1850``

4. Switch flags:

For logical flags, use a flag without an argument -- ``--feature`` for the case where something is off by default and you want to turn it on, or ``--no-feature`` for the case where something is on by default and you want to turn it off -- instead of something like ``--feature true`` / ``--feature false``. Prefer phrasing these as positives whenever possible (e.g., ``--allow-multiple-pfts`` instead of ``--no-single-pft``). To the extent possible, try to clearly document what the behavior is if you add the flag, and what the behavior is if you don't add the flag.

   * Rationale:
      * This use of ``--feature`` / ``--no-feature`` is more common behavior for Unix tools.
      * For something that is either on or off by default, you can see what the default operation is and what you can change just by looking at the available flag names, without needing to read through all of the documentation of default values.
      * Whenever possible avoid the use of --no-feature arguments as they often are confusing. Try to express the same thing using positive wording.
   * Examples: To turn prognostic crop on (or off): ``--crop`` or ``--no-crop``

=========
 Logging
=========

We try to follow the guide at the top of `Python's logging howto <https://docs.python.org/3/howto/logging.html>`_. In particular, print statements should be used for "console output for ordinary usage of a command line script or program"; ``logger.info`` or ``logger.debug`` should be used to "report events that occur during normal operation of a program (e.g. for status monitoring or fault investigation)", etc. The distinction between when to use print and when to use logging can admittedly be a bit subjective, as it comes down to the question of whether the given output is part of the fundamental operation of the script – i.e., part of what the script is designed to do is to give this output. For example, ``run_sys_tests`` prints a variety of information when it starts, particularly concerning the git and manage_externals status of the current repository. The rationale for using ``print`` statements for this is that Bill designed ``run_sys_tests`` to replace some of the repetitive items that he did whenever running the system tests. One of these items was running ``git status`` and ``./manage_externals/checkout_externals -S -v`` to check that the repository is in a clean state. Thus, in this case, Bill's view is that the output from these commands is part of the fundamental purpose of ``run_sys_tests``: it is something he always wants to see, and he feels that it is important for anyone running the system tests to review, and thus ``print`` statements are appropriate here.

Near the top of each python module where logging is used, there should be a line, ``logger = logging.getLogger(__name__)``. Then logging statements should be done using statements like ``logger.info(...)``, *not* ``logging.info(...)``: this allows more contextual information in logging output.
