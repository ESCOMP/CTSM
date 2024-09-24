.. include:: ../substitutions.rst

.. _running-prognostic-crop-model:

===================================
 Running the prognostic crop model
===================================

The prognostic crop model is setup to work with CLM4.0, CLM4.5 or |version| with either BGC or CN (with or without DV). In order to use the initial condition file, we need to set the ``RUN_TYPE`` to startup rather than ``hybrid`` since the compset for f19 sets up to use an initial condition file without crop active. To activate the crop model you can choose a compset that has "Crop" in the name such as "I1850Clm50BgcCropCru" or simply add ``-crop`` to ``CLM_BLDNML_OPTS`` (or for CLM4.0 add ``-crop on`` to ``CLM_CONFIG_OPTS``).

Example: Crop Simulation
------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case CROP -res f19_g17_gl4 -compset I1850Clm50BgcCropCru
   > cd CROP

   > ./case.setup

   # Now build and run normally
   > ./case.build
   > ./case.submit
