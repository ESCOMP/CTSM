.. include:: ../substitutions.rst

.. _spinning-up-clm-bgc:

=============================
 Spinup of |version|-BGC-Crop
=============================

To get the |version|-BGC model to a steady state, you start it from arbitrary initial conditions using the "accelerated decomposition spinup" (``CLM_ACCELERATED_SPINUP on`` in CLM ``env_run.xml``, see example below) mode for 300-400 simulation years. :numref:`Figure BGC AD spinup plot for 1850` shows spinup behavior for an 1850 BGC accelerated decomposition (AD) case using CRUJRA atmospheric forcing. Generally, the criterion that less than 3% of the land surface be in total ecosystem carbon disequilibrium takes the longest to satisfy due to slow soil carbon (TOTSOMC) turnover times in the Arctic.

.. _Figure BGC AD spinup plot for 1850:

.. figure:: AD_spinup_placeholder.png

 BGC AD spinup plot for a year 1850 case with CRUJRA atmospheric forcing. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

After this you continue in "SASU" mode (``CLM_ACCELERATED_SPINUP sasu`` in CLM ``env_run.xml``, see example below), and run for 300-350 simulation years. :numref:`Figure BGC SASU spinup plot for 1850` shows spinup behavior for an 1850 BGC SASU case using CRUJRA atmospheric forcing. The criterion that less than 3% of the land surface be in total ecosystem carbon disequilibrium takes the longest to satisfy and need not be met for this step.

.. _Figure BGC SASU spinup plot for 1850:

.. figure:: SASU_spinup_placeholder.png

 BGC SASU spinup plot for a year 1850 case with CRUJRA atmospheric forcing and initialization from the end of the BGC AD spinup case. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

After this you continue in standard mode for about 160 years. We refer to this as the post-SASU or pSASU phase (``CLM_ACCELERATED_SPINUP off`` in CLM ``env_run.xml``, see example below). :numref:`Figure BGC pSASU spinup plot for 1850` shows spinup behavior for an 1850 BGC pSASU case using CRUJRA atmospheric forcing. As before, the criterion that less than 3% of the land surface be in total ecosystem carbon disequilibrium takes the longest to satisfy.

.. _Figure BGC pSASU spinup plot for 1850:

.. figure:: pSASU_spinup_placeholder.png

 BGC pSASU spinup plot for a year 1850 case with CRUJRA atmospheric forcing and initialization from the end of the BGC SASU spinup case. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

As an alternative to spinning up, one may start from a default initial file that is setup as part of the selected compset. When the simulation's spatial resolution is identical to the initial file's resolution, it may still take 10 or more years for variables such as TLAI (total leaf area index), GPP (gross primary production), and TWS (total water storage) to reach a new equilibrium state due to the different atmospheric forcing. Similarly, it may take 10 or more years for these variables to reach a new equilibrium when switching atmospheric CO2 from 1850 to a present-day value.

Example: AD_spinup Simulation for |version|-BGC
--------------------------------------------------------
For the first step of running in ``CLM_ACCELERATED_SPINUP on`` mode, you will setup a case, and then edit the values in env_build.xml and env_run.xml so that the right configuration is turned on and the simulation is setup to run for the required length of simulation time. Try the following:

::

   > cd cime/scripts
   > ./create_newcase -case AD_spinup -res f19_g17 -compset I1850Clm60BgcCrop
   > cd AD_spinup
   # Change accelerated spinup mode
   > ./xmlchange CLM_ACCELERATED_SPINUP="on"
   # Now setup
   > ./case.setup
   # Now build
   > ./case.build
   # The following sets RESUBMIT to 7 times and
   # STOP_DATE, STOP_N, and STOP_OPTION to Jan/1/0401, 50, "nyears" in env_run.xml (you could also modify these with an editor)
   > ./xmlchange RESUBMIT=7,STOP_N=50,STOP_OPTION=nyears,STOP_DATE=04010101
   # Now run normally
   > ./case.submit

Save the last restart file from this simulation to use in the next step.

.. _eg-sasu-spinup:

Example: SASU_spinup Simulation for |version|-BGC
------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case SASU_spinup -res f19_g17 -compset I1850Clm60BgcCrop
   > cd SASU_spinup
   # Copy the last restart files from the AD_spinup case into your run directory
   > cp /scratch/archive/AD_spinup/rest/0401-01-01-00000/* /scratch/SASU_spinup/run
   # Runtype should already be startup, but this will ensure it
   > ./xmlchange RUN_TYPE=startup
   # Set finidat to the restart file copied in the previous step
   > echo ' finidat = "AD_spinup.clm2.r.0401-01-01-00000.nc"' > user_nl_clm
   # Now setup
   > ./case.setup
   > Now build
   > ./case.build
   # The following sets RESUBMIT to 6 times and
   # STOP_N to 50 and STOP_OPTION to "nyears" in env_run.xml (you could also modify these with an editor)
   > ./xmlchange RESUBMIT=6,STOP_OPTION=nyears,STOP_N=50
   > Now run as normal
   > ./case.submit

Save the last restart file from this step and use it as the ``finidat`` file for the pSASU step. Save the restart file from the end of the pSASU simulation to use as a "finidat" file for future simulations.


Example: pSASU_spinup Simulation for |version|-BGC
------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case pSASU_spinup -res f19_g17 -compset I1850Clm60BgcCrop
   > cd pSASU_spinup
   # Copy the last restart files from the SASU_spinup case into your run directory
   > cp /scratch/archive/SASU_spinup/rest/0351-01-01-00000/* /scratch/pSASU_spinup/run
   # Runtype should already be startup, but this will ensure it
   > ./xmlchange RUN_TYPE=startup
   # Set finidat to the restart file copied in the previous step
   > echo ' finidat = "SASU_spinup.clm2.r.0351-01-01-00000.nc"' > user_nl_clm
   # Now setup
   > ./case.setup
   > Now build
   > ./case.build
   # The following sets RESUBMIT to 3 times and
   # STOP_N to 50 and STOP_OPTION to "nyears" in env_run.xml (you could also modify these with an editor)
   > ./xmlchange RESUBMIT=3,STOP_OPTION=nyears,STOP_N=50
   > Now run as normal
   > ./case.submit

