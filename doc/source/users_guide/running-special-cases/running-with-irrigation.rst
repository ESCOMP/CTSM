.. running-with-irrigation:

===================================
 Running with irrigation
===================================

In CLM4.5 irrigation can ONLY be used WITH crop. 
To turn on irrigation in CLM4.5 we simply add "-irrig on" to ``CLM_BLDNML_OPTS``. 
Just as in the crop example we also change ``RUN_TYPE`` to ``startup`` so that we don't use an initial condition file that is incompatible with irrigation.

Example: Irrigation Simulation
------------------------------------------
::

   # Note here we do a CLMSP simulation as that is what has been validated
   > cd scripts
   > ./create_newcase -case IRRIG -res f19_g16 -compset I 
   > cd IRRIG

   # Append "-irrig" to CLM_BLDNML_OPTS in env_run.xml (you could also use an editor)
   > ./xmlchange CLM_BLDNML_OPTS="-irrig" -append

   # Change to startup type so uses spunup initial conditions file for irrigation if it exists
   # By default the model will do a hybrid startup with an initial condition file
   # incompatible with the irrigation surface dataset.
   > ./xmlchange RUN_TYPE=startup
   > ./case.setup

   # Now build and run normally
   > ./case.build
   > ./case.submit


