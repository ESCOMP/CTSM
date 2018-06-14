.. running-prognostic-crop-model:

===================================
 Running the prognostic crop model
===================================

The prognostic crop model is setup to work with +|version|-BGC or CLM4.0-CN (with or without DV) for present day conditions and we have surface and initial condition datasets at f19 resolution. 
In order to use the initial condition file, we need to set the ``RUN_TYPE`` to startup rather than ``hybrid`` since the compset for f19 sets up to use an initial condition file without crop active. 
To activate the crop model you can choose a compset that has "CROP" in the name such as "ICRUCLM45BGCCROP" or simply add "-crop on" to ``CLM_CONFIG_OPTS``.

Example: Crop Simulation
------------------------------------
::

   > cd scripts
   > ./create_newcase -case CROP -res f19_g17_gl4 -compset I1850Clm50BgcCropCru 
   > cd CROP

   # Append "-crop on" to CLM_CONFIG_OPTS in env_build.xml (you could also use an editor)
   > ./xmlchange CLM_CONFIG_OPTS="-crop on" -append

   # Change to startup type so uses spunup initial conditions file for crop if it exists
   # By default the model will do a hybrid startup with an initial condition file
   # incompatible with the crop surface dataset.

   > ./xmlchange RUN_TYPE=startup
   > ./case.setup

   # Now build and run normally
   > ./case.build
   > ./case.submit
