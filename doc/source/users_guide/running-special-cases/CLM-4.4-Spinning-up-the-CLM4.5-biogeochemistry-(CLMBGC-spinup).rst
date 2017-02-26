.. _spinning-up-clm45-bgc:

======================
 Spinup of CLM4.5-BGC
======================

To get the CLM4.5-BGC model to a steady state, you first run it from arbitrary initial conditions using the "accelerated decomposition spinup" (-bgc_spinup on in CLM **configure**) mode for 1000 simulation years. 
After this you branch from this mode in the "final spinup" (-bgc_spinup off in CLM **configure**), and run for (at least 200+ simulation years).

**1. 45_AD_SPINUP**
     For the first step of running 1000+ years in "-bgc_spinup on" mode, you will setup a case, and then edit the values in env_build.xml and env_run.xml so that the right configuration is turned on and the simulation is setup to run for the required length of simulation time. So do the following:
   
Example:: AD_SPINUP Simulation for CLM4.5-BGC
--------------------------------------------------------
::

   > cd scripts
   > ./create_newcase -case BGC_spinup -res f19_g16 -compset I1850CRUCLM45BGC -mach yellowstone_intel
   > cd BGC_spinup
   # Append "-spinup on" to CLM_BLDNML_OPTS
   > ./xmlchange CLM_BLDNML_OPTS="-bgc_spinup on" -append
   # The following sets CLM_FORCE_COLDSTART to "on", and run-type to startup (you could also use an editor)
   > ./xmlchange CLM_FORCE_COLDSTART=on,RUN_TYPE=startup
   # Make the output history files only annual, by adding the following to the user_nl_clm namelist
   > echo 'hist_nhtfrq = -8760' >> user_nl_clm
   # Now setup
   > ./cesm_setup -case
   # Now build
   > ./BGC_spinup.build
   # The following sets RESUBMIT to 30 times in env_run.xml (you could also use an editor)
   # The following sets STOP_DATE,STOP_N and STOP_OPTION to Jan/1/1001, 20, "nyears" in env_run.xml (you could also use an       editor)
   > ./xmlchange RESUBMIT=20,STOP_N=50,STOP_OPTION=nyears,STOP_DATE=10010101
   # Now run normally
   > ./BGC_spinup.submit

.. note:: This same procedure works for CLM4.5-CN as well, you can typically shorten the spinup time from 1000 years to 600 though.

Afterwards save the last restart file from this simulation to use in the next step.

**2. Final spinup for CLM4.5-BGC**
     Next save the last restart file from this step and use it as the "finidat" file to use for one more spinup for at least 200+ years in normal mode. So do the following:

Example: Final CLMBGC Spinup Simulation for CLM4.5-BGC
------------------------------------------------------------------
::

   > cd scripts
   > ./create_newcase -case BGC_finalspinup -res f19_g16 -compset I1850CRUCLM45BGC -mach yellowstone_intel
   > cd BGC_finalspinup
   # Now, Copy the last CLM restart file from the earlier case into your run directory
   > cp /ptmp/$LOGIN/archive/BGC_spinup/rest/BGC_spinup.clm*.r*.1002-01-01-00000.nc \
   /glade/scratch/$LOGIN/CN_finalspinup
   # Set the runtype to startup
   > ./xmlchange RUN_TYPE=startup
   # And copy the rpointer files for datm and drv from the earlier case
   > cp /glade/scratch/$LOGIN/archive/BGC_spinup/rest/rpointer.atm /ptmp/$LOGIN/CN_finalspinup
   > cp /glade/scratch/$LOGIN/archive/BGC_spinup/rest/rpointer.drv /ptmp/$LOGIN/CN_finalspinup
   # Set the finidat file to the last restart file saved in previous step
   > echo ' finidat = "BGC_spinup.clm2.r.1002-01-01-00000.nc"' > user_nl_clm
   # Now setup
   > ./cesm_setup
   > Now build
   > ./BGC_finalspinup.build
   # The following sets RESUBMIT to 4 times in env_run.xml (you could also use an editor)
   # The following sets STOP_N and STOP_OPTION to 50 and "nyears" in env_run.xml (you could also use an editor)
   > ./xmlchange RESUBMIT=4,STOP_OPTION=nyears,STOP_N=50
   > Now run as normal
   > ./BGC_finalspinup.submit

To assess if the model is spunup plot trends of CLMBGC variables of interest. If you see a trend, you may need to run the simulation longer. Finally save the restart file from the end of this simulation to use as an "finidat" file for future simulations.

.. note:: This same final spinup procedure works for CLM4.5-CN as well, you can typically shorten the spinup time from 200 years to 50 though.


   
