.. sectnum::

.. contents::

==================================
 Overview of this design document
==================================

This document describes various design considerations for writing pFUnit-based unit tests.

=================
 Handling errors
=================

When building the unit tests, we replace the standard ``shr_abort_mod`` with a version that throws a pFUnit exception instead of actually aborting. This is important for expected-error testing, where we want to check that the code aborts when it should: by throwing a pFUnit exception, we can check for this exception in the test. If we actually aborted, we wouldn't be able to catch and check this, and all remaining tests in the given module would be skipped.

This behavior can sometimes be a problem, though, because Fortran doesn't have real exception handling. So when a pFUnit exception is "thrown", the code continues along as if nothing abnormal happened. This can lead to later errors, such as accessing unallocated variables, array bounds exceptions, floating point exceptions, etc. This means that you may not be able to do some expected error testing that you'd like to do, or in order to do this expected error testing, you'll need to add some code to prevent errors after calling ``endrun``. (This code after ``endrun`` won't be executed in the production case, but would be executed in unit testing.)

In unit test-specific code, if we detect an unexpected error - which in this case would typically indicate an issue in the unit test code, rather than an issue in the production code - we use different techniques in different places. Where possible, we prefer calling ``shr_sys_abort`` if an unexpected error occurs in unit test-specific code. However, there are some situations where this is problematic, and we instead use ``stop 1``. (Note that ``stop 1`` causes the program to abort with a non-zero error code, which pFUnit interprets as a test failure, as is desired in these cases; in contrast, ``stop`` with no argument, an argument of 0, or a string argument causes the program to abort with a zero error code, which pFUnit interprets as a pass, which is often not what you want.) Some situations where we use ``stop 1`` instead of ``shr_sys_abort`` are:

1. Unit test-specific code where throwing a pFUnit exception and continuing will lead to cryptic errors later (e.g., array bounds exceptions, etc.).

2. Code called using pFUnit's ``EXTRA_FINALIZE`` technique - i.e., code called at the end of all tests in a given pFUnit executable. Throwing a pFUnit exception in these cases doesn't lead to a FAIL result, so we use a ``stop 1`` to ensure that we get a FAIL result if there's a problem in this finalization.

(See https://github.com/ESCOMP/CTSM/pull/3017/commits/c4377dda9b45a9649f53cd0981657ee6eb268e8b for an illustration of some of these cases.)

======================================
 ESMF initialization and finalization
======================================

Some unit tests require ESMF to be initialized. It is an error to re-initialize ESMF after finalizing it, so we ensure that initialization and finalization are only done once.

To ensure that ESMF initialization is only done once, any call to ``ESMF_Initialize`` is within a conditional on the return value of ``ESMF_IsInitialized``. (ESMF appears to allow multiple initializations, with subsequent initializations not doing anything, but this behavior may not be guaranteed in the future, so to be safe, we ensure that initialization is only done once in unit testing.)

Handling ESMF finalization is trickier, because we only want to do it at the end of a given pFUnit executable. (In practice, it might be okay to skip ESMF finalization, but this might cause issues like un-flushed log files.) To accomplish this, we use the ``EXTRA_FINALIZE`` argument to ``add_pfunit_ctest`` in the ``CMakeLists.txt`` file for any unit test that might initialize ESMF (as indicated by including ``esmf`` in the ``LINK_LIBRARIES``. (e.g., see https://github.com/ESCOMP/CTSM/pull/3017/commits/55b373b184ac56c4c1f9b837f8556cfc9ef33339.) This argument specifies a subroutine that pFUnit calls after the last test in a given executable. The associated ``EXTRA_USE`` argument gives the module name in which this subroutine is defined. Note that the ``ESMF_Finalize`` call is only done if ``ESMF_IsInitialized`` returns ``.true.``, so it's safe to use this in situations where ESMF may or may not have been initialized.

We could handle ESMF initialization in the same way as finalization. But this might lead to unnecessary ``ESMF_Initialize`` calls. So, currently, we invoke ``ESMF_Initialize`` from the unit tests that need it, so that it's only done when truly needed by a unit test.