.. include:: ../substitutions.rst

.. _creating-maps-for-mksurfdata:

*********************************************
Creating input for surface dataset generation
*********************************************

Generating ESMF mesh files
==================================

The ``mksurfdata_esmf`` tool requires ESMF mesh files to describe the input and output grids used for generating fsurdat and landuse files. CLM? provides the ``mesh_maker`` tool for generating such files.

ESMF mesh files for all the standard model resolutions and the raw surface datasets already exist and the files are in the XML database. Hence, you may skip this step -- UNLESS YOU ARE CREATING YOUR OWN GRIDS.

