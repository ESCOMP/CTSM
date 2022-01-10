.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This documents our conventions for command-line arguments, primarily as they pertain to python scripts. These primarily focus on user-visible aspects, as opposed to implementation details.

========================================
 Conventions for command-line arguments
========================================

1. Options that are more than a single character should be formatted as ``--some-variable`` -- i.e., with a double leading hyphen and a hyphen as word separator.

   * Rationale: This follows the `POSIX conventions <https://www.gnu.org/software/libc/manual/html_node/Argument-Syntax.html>`_.

2. Scripts should support ``--verbose`` / ``-v`` and ``--debug`` options. See comments at the top of ``ctsm_logging.py`` for details.

3. For arguments that need a value: If there is no obvious default for the value, then the argument should be required rather than being given a default that the user will often need to override. A general rule of thumb is: if, say, 75% -- 90% of users will want a particular value, then it can be okay to have this value as the default; otherwise, it should be a required argument. This holds even for arguments that are only needed in certain circumstances: For example, if argument ``--fooval FOOVAL`` needs to be given whenever ``--bar`` is given, but there is no obvious default for ``FOOVAL``, then the argument parsing should be coded so that ``--fooval`` does not have a default, and there is error checking code that ensures that ``fooval`` has been specified (not ``None``) if ``--bar`` was given.

   * Rationale: Giving an argument a default value means that we can't catch and report when a user forgets to set it but should have set it. Leaving off the default allows for this error checking.

4. For logical flags, use a flag without an argument -- ``--feature`` for the case where something is off by default and you want to turn it on, or ``--no-feature`` for the case where something is on by default and you want to turn it off -- instead of something like ``--feature true`` / ``--feature false``. Prefer phrasing these as positives whenever possible (e.g., ``--allow-multiple-pfts`` instead of ``--no-single-pft``). To the extent possible, try to clearly document what the behavior is if you add the flag, and what the behavior is if you don't add the flag.

   * Rationale:

      * This use of ``--feature`` / ``--no-feature`` is more common behavior for Unix tools.
      * For something that is either on or off by default, you can see what the default operation is and what you can change just by looking at the available flag names, without needing to read through all of the documentation of default values.
