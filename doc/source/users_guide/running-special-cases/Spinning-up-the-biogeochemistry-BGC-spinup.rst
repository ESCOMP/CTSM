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

 BGC pSASU spinup plot for a year 1850 case with CRUJRA atmospheric forcing and initialization from the end of the BGC AD spinup case. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

As an alternative to spinning up, you may start from a default initial file that is setup as part of the selected compset. :numref:`Figure BGC initialized spinup plot for 1850` shows spinup behavior for an 1850 pSASU BGC case that loops over one year of coupler history output for atmospheric forcing (generated from the fully coupled model), initialized with a BGC initial file generated from a CRUJRA atmospheric forcing case. Note that it takes about 10 years for variables such as TLAI (total leaf area index), GPP (gross primary production), and TWS (total water storage) to reach a specified equilibrium state (denoted by the dotted lines) due to the different atmospheric forcing.

.. _Figure BGC initialized spinup plot for 1850:

.. figure:: initialized_spinup_placeholder.png

 BGC initialized spinup plot for year 1850. Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

:numref:`Figure BGC initialized spinup plot for 2000 CO2` shows spinup behavior for the same case but also changes CO2 to present-day conditions (379 ppmv). Again, it takes about 10 years to reach equilibrium for TLAI, GPP, and TWS.

.. _Figure BGC initialized spinup plot for 2000 CO2:

.. figure:: initialized_2000_placeholder.png

 BGC initialized spinup plot for year 2000 CO2.  Variables examined are TOTECOSYSC (total ecosystem carbon), TOTSOMC (total soil organic matter carbon), TOTVEGC (total vegetation carbon), TLAI (total leaf area index), GPP (gross primary production) and TWS (total water storage). Generated using .../tools/contrib/SpinupStability_BGC_v11.ncl.

If you use the default initial file and you signficantly change model behavior or atmospheric forcing, and you are concerned about the carbon equilibrium (e.g., TOTECOSYSC, TOTSOMC, TOTVEGC), particularly at high latitudes, then we recommend you put the model back into AD mode to reach a new equilibrium. In this configuration, this will also automatically reseed "dead" plant functional types in the initial file with a bit of leaf carbon to give those plant functional types another chance to grow under the new atmospheric forcing or model conditions.

**1. |version| accelerated-decomposition (AD) spinup**
     For the first step of running 200+ years in ``-bgc_spinup on`` mode, you will setup a case, and then edit the values in env_build.xml and env_run.xml so that the right configuration is turned on and the simulation is setup to run for the required length of simulation time. So do the following:

Example: AD_SPINUP Simulation for |version|-BGC
--------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case BGC_spinup -res f19_g17_gl4 -compset I1850Clm50BgcCropCru
   > cd BGC_spinup
   # Change accelerated spinup mode
   > ./xmlchange CLM_ACCELERATED_SPINUP="on"
   # Now setup
   > ./case.setup -case
   # Now build
   > ./case.build
   # The following sets RESUBMIT to 3 times in env_run.xml (you could also use an editor)
   # The following sets STOP_DATE,STOP_N and STOP_OPTION to Jan/1/0201, 20, "nyears" in env_run.xml (you could also use an editor)
   > ./xmlchange RESUBMIT=3,STOP_N=50,STOP_OPTION=nyears,STOP_DATE=02010101
   # Now run normally
   > ./case.submit

.. note:: This same procedure works for |version|-CN as well.

Afterwards save the last restart file from this simulation to use in the next step.

**2. Final spinup for |version|-BGC**
     Next save the last restart file from this step and use it as the ``finidat`` file to use for one more spinup for at least 400+ years in normal mode. So do the following:

.. _eg-final-clmbgc-spinup:

Example: Final CLMBGC Spinup Simulation for |version|-BGC
------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case BGC_finalspinup -res f19_g17_gl4 -compset I1850Clm50BgcCropCru
   > cd BGC_finalspinup
   # Now, Copy the last CLM restart file from the earlier case into your run directory
   > cp /ptmp/$LOGIN/archive/BGC_spinup/rest/BGC_spinup.clm*.r*.0201-01-01-00000.nc \
   /glade/scratch/$LOGIN/CN_finalspinup/run
   # Set the runtype to startup
   > ./xmlchange RUN_TYPE=startup
   # And copy the rpointer files for datm and drv from the earlier case
   > cp /glade/scratch/$LOGIN/archive/BGC_spinup/rest/rpointer.atm /glade/scratch/$LOGIN/CN_finalspinup/run
   # Set the finidat file to the last restart file saved in previous step
   > echo ' finidat = "BGC_spinup.clm2.r.0201-01-01-00000.nc"' > user_nl_clm
   # Now setup
   > ./case.setup
   > Now build
   > ./case.build
   # The following sets RESUBMIT to 7 times in env_run.xml (you could also use an editor)
   # The following sets STOP_N and STOP_OPTION to 50 and "nyears" in env_run.xml (you could also use an editor)
   > ./xmlchange RESUBMIT=7,STOP_OPTION=nyears,STOP_N=50
   > Now run as normal
   > ./case.submit

To assess if the model is spunup, plot trends for CLMBGC variables of interest using .../tools/contrib/SpinupStability.ncl. If you don't meet the equilibrium criteria, you may need to run the simulation longer. Finally save the restart file from the end of this simulation to use as an "finidat" file for future simulations.

.. note:: This same final spinup procedure works for |version|-CN as well.

