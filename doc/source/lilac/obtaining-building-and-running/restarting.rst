.. highlight:: shell

.. _restarting:

=====================================
 Continuing a run from restart files
=====================================

All of the information that CTSM and LILAC need to continue a run from restart files is
given in the restart files themselves. No namelist changes need to be made (other than
whatever is needed in the host atmosphere model), but the ``starttype_in`` argument to the
``lilac_init2`` subroutine call from the atmosphere model will need to be changed to
"continue" rather than "startup".

CTSM and LILAC use ``rpointer`` files to indicate the specific restart files that should
be read. These files, ``rpointer.lnd`` and ``rpointer.lilac``, are one-line text files
that simply specify the name of the respective restart files. When restart files are
written (according to the ``write_restarts_now`` argument to the ``lilac_run``
subroutine), these ``rpointer`` files are updated to point to the latest set of restarts.

If you want to restart from the latest set of restart files, the ``rpointer`` files should
already be set up to facilitate this. However, if you want to restart from an earlier set
of restarts, you can simply edit ``rpointer.lnd`` and ``rpointer.lilac`` to point to the
appropriate restart files.

.. important::

   Be sure that the ``rpointer.lnd`` and ``rpointer.lilac`` files point to restart files
   from the same time as each other, and from the same time as the atmosphere model's
   restart time.
