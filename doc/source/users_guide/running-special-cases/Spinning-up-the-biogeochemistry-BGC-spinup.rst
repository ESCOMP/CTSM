.. include:: ../substitutions.rst

.. _spinning-up-clm-bgc:

=============================
 Spinup of CLM-BGC-Crop
=============================

To get the CLM-BGC-Crop model to a steady state, you start it from arbitrary initial conditions using the "accelerated decomposition spinup" (``CLM_ACCELERATED_SPINUP on`` in CLM `env_run.xml`, see example below) mode for 300-400 simulation years. :numref:`Figure BGC-Crop AD spinup plot for 1850` shows spinup behavior for an 1850 BGC-Crop accelerated decomposition (AD) case using CRUJRA atmospheric forcing. Generally, the criterion that less than 3% of the land surface be in total ecosystem carbon disequilibrium takes the longest to satisfy due to slow soil carbon (TOTSOMC) turnover times in the Arctic.

.. _Figure BGC-Crop AD spinup plot for 1850:

.. figure:: ctsm5.4.CMIP7_ciso_ctsm5.3.075_f09_124_AD_Spinup-0.png

 BGC-Crop AD spinup plot for a year 1850 case with CRUJRA atmospheric forcing. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

After this you continue in "SASU" mode (``CLM_ACCELERATED_SPINUP sasu`` in CLM `env_run.xml`, see example below), and run for 300-350 simulation years. :numref:`Figure BGC-Crop SASU spinup plot for 1850` shows spinup behavior for an 1850 BGC-Crop SASU case using CRUJRA atmospheric forcing. The criterion that less than 3% of the land surface be in total ecosystem carbon disequilibrium takes the longest to satisfy and need not be met for this step.

.. _Figure BGC-Crop SASU spinup plot for 1850:

.. figure:: ctsm5.4.CMIP7_ciso_ctsm5.3.075_f09_124_SASU_Spinup-0.png

 BGC-Crop SASU spinup plot for a year 1850 case with CRUJRA atmospheric forcing and initialization from the end of the BGC-Crop AD spinup case. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

After this you continue in standard mode for 200 years. We refer to this phase as post-SASU, pSASU, or normal mode (``CLM_ACCELERATED_SPINUP off`` in CLM `env_run.xml`, see example below). :numref:`Figure BGC-Crop normal mode plot for 1850` shows spinup behavior for an 1850 BGC-Crop normal mode case using CRUJRA atmospheric forcing. As before, the criterion that less than 3% of the land surface be in total ecosystem carbon disequilibrium takes the longest to satisfy.

.. _Figure BGC-Crop normal mode plot for 1850:

.. figure:: ctsm5.4.CMIP7_ciso_ctsm5.3.075_f09_124_pSASU_Spinup-0.png

 BGC-Crop normal mode plot for a year 1850 case with CRUJRA atmospheric forcing and initialization from the end of the BGC-Crop SASU spinup case. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

As an alternative to spinning up, one may start from a default initial file that is setup as part of the selected compset. When the simulation's spatial resolution is identical to the initial file's resolution, it may still take 10 or more years for variables such as TLAI (total leaf area index), GPP (gross primary production), and TWS (total water storage) to reach a new equilibrium state due to the different atmospheric forcing. Similarly, it may take 10 or more years for these variables to reach a new equilibrium when switching atmospheric CO2 from 1850 to a present-day value.

Example: AD_spinup Simulation for CLM-BGC-Crop
--------------------------------------------------------
For the first step of running in ``CLM_ACCELERATED_SPINUP on`` mode, you will setup a case, and then edit the values in env_build.xml and env_run.xml so that the right configuration is turned on and the simulation is setup to run for the required length of simulation time. Try the following:

::

   > cd cime/scripts
   > ./create_newcase -case AD_spinup -res f19_g17 -compset I1850Clm60BgcCrop --run-unsupported
   > cd AD_spinup
   # Change accelerated spinup mode
   > ./xmlchange CLM_ACCELERATED_SPINUP="on"
   # Now setup and build
   > ./case.setup
   > ./case.build
   # The following sets RESUBMIT to 7 times and
   # STOP_DATE, STOP_N, and STOP_OPTION to Jan/1/0401, 50, "nyears" in env_run.xml (you could also modify these with an editor)
   > ./xmlchange RESUBMIT=7,STOP_N=50,STOP_OPTION=nyears,STOP_DATE=04010101
   # The following makes sure that we run with MOSART off
   > ./xmlchange MOSART_MODE=NULL
   # Now run
   > ./case.submit

While this simulation progresses, use SpinupStability_BGC_v11.ncl to assess whether the simulation is approaching equilibrium. When the simulation ends, save the last restart file for use in the SASU_spinup step.

Using the SpinupStability.ncl scripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In CLM's /tools/contrib directory there are three versions of this .ncl script:

- SpinupStability_BGC_v11.ncl for Bgc and BgcCrop compsets run on 2D lat/lon grids.
- SpinupStability_BGC_v12_SE.ncl for Bgc, BgcCrop, or Fates compsets run on certain spectral element grids (currently ne120, ne30, ne16).
- SpinupStability_SP_v10.ncl for Sp compsets run on 2D lat/lon grids. See section :numref:`spinning-up-sp` for helpful pointers about this script that may also apply to the BGC-Crop versions.

To run one of these scripts on derecho, one loads ncl (``module load ncl``) and submits with ``ncl SpinupStability_BGC_v11.ncl``, for example. Before running one needs to confirm a few easy settings appearing near the top of each script.

One of the settings that may not be intuitive at first glance is ``annual_hist``. By default the phases AD_spinup and SASU_spinup generate annual history, so set this to "True", while normal mode generates monthly history, so set this to "False".

.. _eg-sasu-spinup:

Example: SASU_spinup Simulation for CLM-BGC-Crop
------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case SASU_spinup -res f19_g17 -compset I1850Clm60BgcCrop --run-unsupported
   > cd SASU_spinup
   # Change accelerated spinup mode
   > ./xmlchange CLM_ACCELERATED_SPINUP="sasu"
   # Now setup
   > ./case.setup
   # Copy the last restart files from the AD_spinup case into your run directory
   # On NSF-NCAR's derecho computer, cd to /glade/derecho/scratch/$USER
   > cp archive/AD_spinup/rest/0401-01-01-00000/* SASU_spinup/run
   # Runtype should already be startup, but this will ensure it
   > ./xmlchange RUN_TYPE=startup
   # Set finidat to the restart file copied in the previous step
   > echo ' finidat = "AD_spinup.clm2.r.0401-01-01-00000.nc"' > user_nl_clm
   # Now build
   > ./case.build
   # The following sets RESUBMIT to 6 times and
   # STOP_N to 50 and STOP_OPTION to "nyears" in env_run.xml (you could also modify these with an editor)
   > ./xmlchange RESUBMIT=6,STOP_OPTION=nyears,STOP_N=50
   # The following makes sure that we run with MOSART off
   > ./xmlchange MOSART_MODE=NULL
   # Now run
   > ./case.submit

Save the last restart file from this step and use it as the ``finidat`` file for the normal mode simulation. Save the restart file from the end of the normal mode simulation to use as a "finidat" file for future simulations.


Example: Normal mode simulation for CLM-BGC-Crop
--------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case pSASU_spinup -res f19_g17 -compset I1850Clm60BgcCrop --run-unsupported
   > cd pSASU_spinup
   # Change accelerated spinup mode
   > ./xmlchange CLM_ACCELERATED_SPINUP="off"
   # Now setup
   > ./case.setup
   # Copy the last restart files from the SASU_spinup case into your run directory
   # On NSF-NCAR's derecho computer, cd to /glade/derecho/scratch/$USER
   > cp archive/SASU_spinup/rest/0351-01-01-00000/* pSASU_spinup/run
   # Runtype should already be startup, but this will ensure it
   > ./xmlchange RUN_TYPE=startup
   # Set finidat to the restart file copied in the previous step
   > echo ' finidat = "SASU_spinup.clm2.r.0351-01-01-00000.nc"' > user_nl_clm
   # Now build
   > ./case.build
   # The following sets RESUBMIT to 3 times and
   # STOP_N to 50 and STOP_OPTION to "nyears" in env_run.xml (you could also modify these with an editor)
   > ./xmlchange RESUBMIT=3,STOP_OPTION=nyears,STOP_N=50
   # The following makes sure that we run with MOSART off
   > ./xmlchange MOSART_MODE=NULL
   # Now run
   > ./case.submit

