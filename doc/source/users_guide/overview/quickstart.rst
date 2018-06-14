.. _quickstart:

============
 Quickstart
============

Running the CLM requires a suite of UNIX utilities and programs and you should make sure you have all of these available before trying to go forward with using it. 
If you are missing one of these you should contact the systems administrator for the machine you wish to run on and make sure they are installed.

List of utilities required for CESM in the "CESM1.2.0 Software/Operating System Prerequisites" section in `http://www.cesm.ucar.edu/models/cesm1.2//cesm/doc/usersguide/book1.html <CLM-URL>`_
- UNIX bash shell (for some of the CLM tools scripts)
- NCL (for some of the offline tools for creating/modifying CLM input datasets see `Chapter 2 <CLM-URL>`_ for more information on NCL)
- Python (optional, needed for PTCLM)
- xsltproc, docbook and docbook utilities (optional, needed to build the Users-Guide)

Before working with CLM4.5 read the QuickStart Guide in the `CESM1.2.0 Scripts User's Guide <CLM-URL>`_. Once you are familiar with how to setup cases for any type of simulation with CESM you will want to direct your attention to the specifics of using CLM.

For some of the details of setting up cases for CLM4.5 read the README and text files available from the "models/lnd/clm/doc" directory (see the "CLM Web pages" section for a link to the list of these files). Here are the important ones that you should be familiar with.

1. `README file <CLM-URL>`_ describing the directory structure.

2. `Quickstart.userdatasets <CLM-URL>`_ file describing how to use your own datasets in the model (also see `the Section called Creating your own single-point/regional surface datasets in Chapter 5 <CLM-URL>`_).

3. `models/lnd/clm/doc/KnownBugs <CLM-URL>`_ file describing known problems in CLM4.5 (that we expect to eventually fix).

4. `models/lnd/clm/doc/KnownLimitationss <CLM-URL>`_ file describing known limitations in CLM4.5 and workarounds that we do NOT expect to fix.

The IMPORTANT_NOTES file talks about important things for users to know about using the model scientifically. It content is given in the next chapter on `"What is scientifically validated and functional in CLM4.5 in CESM1.2.0?" <CLM-URL>`_.

The ChangeLog/ChangeSum talk about advances in different versions of CLM. The content of these files is largely explained in the previous chapter on `"What is new with CLM4.5 in CESM1.2.0 since previous public releases?" <CLM-URL>`_.

Note other directories have README files that explain different components and tools used when running CLM and are useful in understanding how those parts of the model work and should be consulted when using tools in those directories. For more details on configuring and customizing a case with CLM see `Chapter 1 <CLM-URL>`_.

The Quickstart.GUIDE (which can be found in ``models/lnd/clm/doc``) is repeated here.

.. include:: ../../clm5.0/doc/Quickstart.GUIDE
   :literal:
