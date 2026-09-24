.. include:: ../substitutions.rst

.. _running-with-previous-simulation-forcing:

=============================================================
 Running with atmospheric forcing from a previous simulation
=============================================================

Another way that you might want to spinup the model is to run your own simulation for a relatively short period (either a B, E, or F compset) and then use it as forcing for your "I" case later. By only running 20 to 50 years for the fully coupled case, you'll save a substantial amount of computer time rather than running the entire spinup period with a fully coupled model.

The first thing we need to do is to run a fully coupled case and save the atmospheric coupling fields on a three hourly basis. In this example, we will run on derecho and archive the data to a local disk that we can then use in the next simulation.

Example: Fully Coupled Simulation to Create Data to Force Next Example Simulation
----------------------------------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case myB1850 -res f09_g17_gl4 -compset B1850
   > cd myB1850
   > ./case.setup
   # Set the followng auxiliary history settings to true in your user_nl_cpl file
   > cat << EOF > user_nl_cpl
   histaux_atm2med_file1_enabled = .true.
   histaux_atm2med_file2_enabled = .true.
   histaux_atm2med_file3_enabled = .true.
   histaux_atm2med_file4_enabled = .true.
   histaux_atm2med_file5_enabled = .true.
   histaux_atm2med_file5_history_n = 1
   histaux_atm2med_file5_history_option = 'ndays'
   histaux_atm2med_file5_ntperfile = 1
   EOF
   # Now build
   > ./case.build
   # The following sets the archival disk space (you could also use an editor)
   > ./xmlchange DOUT_S_ROOT='/glade/home/$USER/$CASE'
   # Make sure files are archived to disk, but NOT to long term storage
   # (you could also use an editor)
   > ./xmlchange DOUT_S=TRUE,DOUT_L_MS=FALSE
   # Set the run length to run a total of 20 years (you could also use an editor)
   > ./xmlchange RESUBMIT=9,STOP_OPTION=nyears,STOP_N=2
   # Now run as normal
   > ./case.submit

Now we run an I compset forced with the data from the previous simulation using the ``CPLHIST-CESM3`` option to DATM_MODE. See :ref:`cplhistforcing` for more information.

.. _eg-sim-data-from-prev-sim:

Example: Simulation Forced with Data from the Previous Simulation
------------------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case frcwmyB1850 -res f09_f09_mt233 -compset I1850Clm60BgcCropSpinup
   > cd frcWmyB1850
   # By default this compset will use the CPLHIST-CESM3 and give you the default data from a CESM3 spinup simulation.

   # HOWEVER, if you want to use your case you would do something like this:
   # The following sets the directory and casename to point to for atm forcing (you could also use an editor)
   > ./xmlchange DATM_CPLHIST_DIR='$CIME_OUTPUT_ROOT/archive/$DATM_CPLHIST_CASE/cpl/hist'
   > ./xmlchange DATM_CPLHIST_CASE="myB1850"
   # The following sets the align year and years to run over for atm forcing
   #  (you could also use an editor)
   > ./xmlchange DATM_YR_ALIGN="1",DATM_YR_START=1,DATM_YR_END=20
   > ./case.setup
   # Now build and run as normal
   > ./case.build
   > ./case.submit
