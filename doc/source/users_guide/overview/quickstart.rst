.. include:: ../substitutions.rst

.. _quickstart:

============
 Quickstart
============

Running the CLM requires a suite of UNIX utilities and programs and you should make sure you have all of these available before trying to go forward with using it. If you are missing one of these you should contact the systems administrator for the machine you wish to run on and make sure they are installed.

List of utilities required for CESM in the `Software/OS Prerequisites <https://www2.cesm.ucar.edu/models/cesm1.2/cesm/doc/usersguide/x32.html#software_system_prerequisites>`_ section of the CESM User's Guide.

- UNIX bash shell (for some of the CLM tools scripts)
- NCL (for some of the offline tools for creating/modifying CLM input datasets; see :ref:`using-ncl` for more information)
- Python

Before working with |version| read the `CESM QuickStart Guide <https://escomp.github.io/CESM/versions/cesm2.2/html/>`_. Once you are familiar with how to setup cases for any type of simulation with CESM you will want to direct your attention to the specifics of using CLM.

For some of the details of setting up cases for |version| read the README and text files available from the ``$CTSMROOT/doc`` directory (see the "CLM Web pages" section for a link to the list of these files). Here are the important ones that you should be familiar with:

- :ref:`readme` describing the directory structure.
- The IMPORTANT_NOTES file talks about important things for users to know about using the model scientifically. It content is given in the next chapter on :ref:`scientific-validiation`.
- The ChangeLog/ChangeSum talk about advances in different versions of CLM. The content of these files is largely explained in the previous chapter on :ref:`what-is-new-with-clm`.
- The release-clm5.0.ChangeLog gives the specific changes that have gone on the release-clm5.0 branch. clm3_0_ChangeLog, clm4_0_ChangeLog, clm4_5_ChangeLog gives the changes that culimated in that given version of the CLM.

Note other directories have README files that explain different components and tools used when running CLM and are useful in understanding how those parts of the model work and should be consulted when using tools in those directories. For more details on configuring and customizing a case with CLM see :ref:`customizing_section`.

The Quickstart.GUIDE (which can be found in ``$CTSMROOT/doc``) is repeated here.

.. include:: ../../../Quickstart.GUIDE
   :literal:
