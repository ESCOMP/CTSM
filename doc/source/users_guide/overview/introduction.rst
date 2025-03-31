.. include:: ../substitutions.rst

.. _introduction:

**User's Guide to version |version| of the Community Land Model (CLM)**

**Authors: Benjamin Andre, Erik Kluzek, William Sacks**

The National Center for Atmospheric Research (NCAR) is operated by the nonprofit University Corporation for Atmospheric Research (UCAR) under the sponsorship of the National Science Foundation. Any opinions, findings, conclusions, or recommendations expressed in this publication are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

National Center for Atmospheric Research
P. O. Box 3000, Boulder, Colorado 80307-3000

.. _rst_Users_Guide_Introduction:

==============
Introduction
==============

The Community Land Model (|release| in |cesmrelease|) is the latest in a series of global land models developed by the CESM Land Model Working Group (LMWG) and maintained at the National Center for Atmospheric Research (NCAR). This guide is intended to instruct both the novice and experienced user on running CLM. This guide pertains to the latest version |version| in |cesmrelease| available for download from the public release subversion repository as a part of |cesmrelease|. Documentation may be different if you are using an older version, you should either use the documentation for that release version, update to the latest version, or use the documentation inside your own source tree. There is information in the ChangeLog file and in :ref:`what-is-new-with-clm` regarding the changes from previous versions of CESM.

.. note:: This release of |version| in |cesmrelease| includes BOTH CLM4.0 physics and CLM4.5 physics used in previous releases as well as the updated |version| physics. CLM allow you to trigger between the three physics modes. Most often when we refer to CLM4.0 we are referring to the CLM4.0 physics in |version| in |cesmrelease| rather than to a specific version of CLM4.0 (where we would give the exact version). And when we refer to CLM4.5 we are referring to the CLM4.5 physics in |version| in |cesmrelease| rather than to a specific version of CLM4.5. Likewise, when referring to |version| we are referring to the |version| physics in |version| in |cesmrelease|.

The novice user should read the :ref:`overview_section` chapter in detail before beginning work, while the expert user should read :ref:`what-is-new-with-clm` and :ref:`quickstart`, and then use the more detailed sections as reference. Before novice users go onto more technical problems covered in :ref:`customizing_section`, :ref:`using-clm-tools-section`, :ref:`adding-new-resolutions-section`, :ref:`running-special-cases-section` or :ref:`running-single-points`, they should know the material covered in :ref:`overview_section` and be able to replicate some of the examples given there.

All users should read :ref:`how-to-use-this-document` and :ref:`getting-help` to understand the document conventions and the various ways of getting help on using |version|. Users should also read the :ref:`scientific-validiation` section to see if their planned use of the model is something that has been scientifically validated and well tested. Users that are NOT using NCAR machines or our list of well tested machines should also read the What are the UNIX utilities required to use |version|? section to make sure they have all the required UNIX utilities on the system they want to do their work.

Developers that are making changes to CLM either for their own development or for development that they hope will eventually become a part of the main CLM should read :ref:`testing_section`. We have a suite of test scripts that automatically test many different model configurations and namelist options, as well as ensuring things like restarts are bit-for-bit and the like. It's helpful to use these scripts to ensure your changes are working correctly. As well as being a required part of the process to bring in new code developments. And it's far easier to use the automated scripts rather than having to figure out, what to test, how to do it, and then finally do it by hand. If you are using non supported machines you may also want to use the test scripts to make sure your machine is working correctly.

.. _what-is-new-with-clm:

============================
 What is New with |version|
============================

`What's new with |version| science <https://escomp.github.io/ctsm-docs/doc/build/html/tech_note/Introduction/CLM50_Tech_Note_Introduction.html#|version|/>`_ gives a synopsis of the changes to CLM since the CLM4.5 release. More details are given in the CLM ChangeLog files:

- `CLM 3.0 ChangeLog file <https://github.com/ESCOMP/CTSM/blob/master/doc/clm3_0_ChangeLog>`_
- `CLM 4.0 ChangeLog file <https://github.com/ESCOMP/CTSM/blob/master/doc/clm4_0_ChangeLog>`_
- `CLM 4.5 ChangeLog file <https://github.com/ESCOMP/CTSM/blob/master/doc/clm4_5_ChangeLog>`_
- `CLM 5.0 ChangeLog file <https://github.com/ESCOMP/CTSM/blob/master/doc/clm5_0_ChangeLog>`_
- `CTSM 1.0 ChangeLog file <https://github.com/ESCOMP/CTSM/blob/master/doc/ctsm1_0_ChangeLog>`_
- `Latest ChangeLog file <https://github.com/ESCOMP/CTSM/blob/master/doc/ChangeLog>`_

Previous release pages give similar list of changes for previous versions of the model.

.. _users-guide-overview:

==========================
 Overview of User's Guide
==========================

In this introduction we first give a simple guide to understand the document conventions in :ref:`how-to-use-this-document`. The next section, :ref:`what-is-new-with-clm`, gives references to describe the differences between |version| in |cesmrelease| and previous CESM releases both from a scientific as well as a software engineering point of view. For information on previous releases of |version| before |version| in |cesmrelease| see the CESM1.2.2 documentation. The :ref:`quickstart` section is for users that are already experts in using CLM and gives a quickstart guide to the bare details on how to use |version|. It also lists the UNIX utilities required to use |version| and is important if you are running on non-NCAR machines, generic local machines, or machines NOT as well tested by us at NCAR. The :ref:`scientific-validiation` section tells you about what has been extensively tested and scientifically validated (and maybe more importantly) what has NOT. Next we have :ref:`best-practices-for-usage` to detail some of the best practices for using |version| for science. The last introductory section is :ref:`getting-help`, which lists different resources for getting help with |version| and |cesmrelease|.

:ref:`customizing_section` goes into detail on how to setup and run simulations with |version| and especially how to customize cases. Details of cesm_setup modes and build-namelist options as well as namelist options are given in this chapter.

:ref:`using-clm-tools-section` gives instructions on the CLM tools for either CLM4.5 or |version| physics for creating input datasets for use by CLM, for the expert user. There's an overview of what each tool does, and some general notes on how to build the FORTRAN tools. Then each tool is described in detail along with different ways in which the tool might be used. A final section on how to customize datasets for observational sites for very savvy expert users is given as the last section of this chapter.

As a followup to the tools chapter, :ref:`adding-new-resolutions-section` tells how to add files to the XML database for build-namelist to use. This is important if you want to use the XML database to automatically select user-created input files that you have created when you setup new cases with CLM (CLM4.0, CLM4.5 and |version| physics).

In :ref:`running-special-cases-section`, again for the expert user, we give details on how to do some particularly difficult special cases. For example, we give the protocol for spinning up the |version|-BGC and CLMCN models as well as CLM with dynamic vegetation active (CNDV). We give instructions to do a spinup case from a previous case with Coupler history output for atmospheric forcing. We also give instructions on running both the prognostic crop and irrigation models. Lastly we tell the user how to use the DATM model to send historical CO2 data to CLM.

:ref:`running-single-points` outlines how to do single-point or regional simulations using |version|. This is useful to either compare |version| simulations with point observational stations, such as tower sites (which might include your own atmospheric forcing), or to do quick simulations with CLM for example to test a new parameterization. There are several different ways given on how to perform single-point simulations which range from simple sampling of existing inputs to more complex where you create all your own datasets, tying into :ref:`using-clm-tools-section` and also :ref:`adding-new-resolutions-section` to add the files into the build-namelist XML database.

There is also :ref:`pts_mode`, which is useful for running single points as part of the Single Column Atmospheric Model (SCAM).

:ref:`troubleshooting-index` gives some guidance on trouble-shooting problems when using |version|. It doesn't cover all possible problems with CLM, but gives you some guidelines for things that can be done for some common problems.

:ref:`testing_section` goes over the automated testing scripts for validating that the CLM is working correctly. The test scripts run many different configurations and options with CLM4.0 physics as well and |version| physics making sure that they work, as well as doing automated testing to verify restarts are working correctly, and testing at many different resolutions. In general this is an activity important only for a developer of |version|, but could also be used by users who are doing extensive code modifications and want to ensure that the model continues to work correctly.

In the appendices we talk about some issues that are useful for advanced users and developers of |version|.

Finally on Github we give `instructions <https://github.com/ESCOMP/CTSM/wiki/Directions-for-editing-CLM-documentation-on-github-and-sphinx>`_ on how to build the documentation associated with |version| (i.e. how to build this document). This document is included in every CLM distribution and can be built so that you can view a local copy rather than having to go to the CESM website. This also could be useful for developers who need to update the documentation due to changes they have made.

.. _readme:

================================
README file describing |version|
================================

The README (which can be found in ``$CTSMROOT/doc``) is repeated here.

.. include:: ../../../../README
   :literal:

.. _best-practices-for-usage:

================
 Best Practices
================

- |version| includes BOTH the old CLM4.0, CLM4.5 physics AND the new |version| physics and you can toggle between those three. The "standard" practice for CLM4.0 is to run with CN on, and with Qian atmospheric forcing. While the "standard" practice for CLM4.5 is to run with BGC on, and CRUNCEP atmospheric forcing. And finally the "standard" practice for |version| is to run with BGC and Prognostic Crop on, with the MOSART model for river routing, as well as the CISM ice sheet model, and using GSWP3 atmospheric forcing. "BGC" is the new |version| biogeochemistry and include CENTURY-like pools, vertical resolved carbon, as well as Nitrification and de-Nitrification (see :ref:`acronyms-and-terms`).

- When running with CLMCN (either CLM4.0 or |version| physics) or |version|-BGC, it is critical to begin with initial conditions that are provided with the release or to spin the model up following the CN spinup procedure before conducting scientific runs (see :ref:`spinning-up-clm-bgc` or :ref:`spinning-up-sp`). Simulations without a proper spinup will effectively be starting from an unvegetated world. See :ref:`setting-initial-conditions` for information on how to provide initial conditions for your simulation.

- Initial condition files are provided for CLM4.0-CN as before, for fully coupled BCN and offline ICN cases for 1850 and 2000 at finite volume grids: 1deg (0.9x1.25), 2deg (1.9x2.5), and T31 resolutions. We also have interpolated initial conditions for BCN for 1850 and 2000 for two finite volume grids: 10x15, 4x5 and two HOMME grids (ne30np4 and ne120np4). There's also an initial condition file for ICN with the prognostic crop model for 2000 at 2deg resolution, and one with CLMSP for 2000 at 2deg resolution. We also have initial conditions for offline CNDV for 1850. The 1850 initial condition files are in 'reasonable' equilibrium. The 2000 initial condition files represent the model state for the year 2000, and have been taken from transient simulations. Therefore, by design the year 2000 initial condition files do not represent an equilibrium state. Note also that spinning the 2000 initial conditions out to equilibrium will not reflect the best estimate of the real carbon/nitrogen state for the year 2000.

- Initial condition files are also provided for |version| for several configurations and resolutions. For CLM4.5-SP and CLM4.5-BGC with both CRUNCEP and GSWP3 forcing we have initial conditions at 1deg resolution for 1850. For |version|-SP and |version|-BGC-Crop with both CRUNCEP and GSWP3 forcing we have initial conditions at 1deg resolution for 1850. Normally, these files are interpolated to any other resolution that you run at.

- Users can interpolate initial condition files at different resolutions at startup of a CLM4.5 or |version| simulation. And the file created can be stored for later use. Interpolated initial condition files may no longer be in 'reasonable' equilibrium.

- In |version| for both |version|-CN, |version|-BGC, and |version|-BGC-Crop the new fire model requires lightning frequency data, and human population density (both are read inside of CLM). By default we have provided a climatology dataset for lightning frequency and a dataset with coverage from 1850 to 2014 for population density. Both of these datasets are interpolated from the native resolution of the datasets to the resolution you are running the model on. If you are running with an atmosphere model or forcing that is significantly different than present day -- the lightning frequency may NOT appropriately correspond to your atmosphere forcing and fire initiation would be inappropriate.

- Aerosol deposition is a required field to both CLM4.0, CLM4.5 and |version| physics, sent from the atmosphere model. Simulations without aerosol deposition will exhibit unreasonably high snow albedos. The model sends aerosol deposition from the atmospheric model (either CAM or DATM). When running with prescribed aerosol the atmosphere model will interpolate the aerosols from 2-degree resolution to the resolution the atmosphere model is running at.

.. _ctsm_vs_cesm_checkout:

=============================
A CTSM versus a CESM checkout
=============================

The directory structure for |version| is different depending on if it's checked out from |release| or |cesmrelease|. If |version| is checked out from |ctsm_gh| the CLM source code is directly under the top level directory. If |cesmrelease| is checkout out from |cesm_gh| then the CLM source directories are under ``components/clm`` from the top-level directory. We will refer to this directory for the CLM source directories in the User's Guide as ``$CTSMROOT``.

.. _how-to-use-this-document:

========================================================
How To Use This Document
========================================================

Links to descriptions and definitions have been provided in the code below. We use the same conventions used in the CESM documentation as outlined below.

::

   Throughout the document this style is used to indicate shell commands and options, fragments of code, namelist variables, etc. Where examples from an interactive shell session are presented, lines starting with > indicate the shell prompt. A backslash "\" at the end of a line means the line continues onto the next one (as it does in standard UNIX shell). Note that $EDITOR" is used to refer to the text editor of your choice. $EDITOR is a standard UNIX environment variable and should be set on most UNIX systems. Comment lines are signaled with a "#" sign, which is the standard UNIX comment sign as well. $CSMDATA is used to denote the path to the inputdata directory for your CESM data.

   > This is a shell prompt with commands \
   that continues to the following line.
   > $EDITOR filename # means you are using a text editor to edit "filename"
   # This is a comment line

   $CTSMROOT means the path to the root of the CTSM model
