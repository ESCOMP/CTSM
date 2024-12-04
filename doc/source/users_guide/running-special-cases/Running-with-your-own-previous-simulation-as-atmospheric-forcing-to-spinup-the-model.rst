.. include:: ../substitutions.rst

.. _running-with-previous-simulation-forcing:

=============================================================
 Running with atmospheric forcing from a previous simulation
=============================================================

Another way that you might want to spinup the model is to run your own simulation for a relatively short period (either a B, E, or F compset) and then use it as forcing for your "I" case later. By only running 20 to 50 years for the fully coupled case, you'll save a substantial amount of computer time rather than running the entire spinup period with a fully coupled model.

The first thing we need to do is to run a fully coupled case and save the atmospheric coupling fields on a three hourly basis. In this example, we will run on cheyenne and archive the data to a local disk that we can then use in the next simulation.

Example: Fully Coupled Simulation to Create Data to Force Next Example Simulation
----------------------------------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case myB1850 -res f09_g17_gl4 -compset B1850
   > cd myB1850
   > ./case.setup
   # Set histaux_a2x3hr to .true. in your user_nl_cpl output from the atmosphere model
   # will be saved 3 hourly
   echo "histaux_a2x3hr=.true." >> user_nl_cpl
   # edit the driver code in order to save the correct list of fields (see note below)
   > cp ../../models/drv/driver/ccsm_comp_mod.F90 SourceMods/src.cpl
   > $EDITOR SourceMods/src.cpl
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

Now we run an I compset forced with the data from the previous simulation using the ``CPLHISTForcing`` option to DATM_MODE. See :ref:`cplhistforcing` for more information.

.. _eg-sim-data-from-prev-sim:

Example: Simulation Forced with Data from the Previous Simulation
------------------------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case frcwmyB1850 -res f09_g17_gl4 -compset I1850Clm50BgcSpinup
   > cd frcWmyB1850
   # The following sets the casename to point to for atm forcing (you could also use an editor)
   > ./xmlchange DATM_CPLHIST_CASE="myB1850"
   # The following sets the align year and years to run over for atm forcing
   #  (you could also use an editor)
   > ./xmlchange DATM_YR_ALIGN="1",DATM_YR_START=1,DATM_YR_END=20
   # Set the strm_datdir in the namelist_defaults_datm.xml
   # file to the archival path of the case above in the form of: /glade/home/achive/$USER/$DATM_CPLHIST_CASE/cpl/hist
   # NOTE: THIS WILL CHANGE THE PATH FOR ALL I1850Clm50BgcSpinup COMPSET CASES MADE AFTER THIS!
   > $EDITOR ../../models/atm/datm/bld/namelist_files/namelist_defaults_datm.xml
   > ./case.setup
   # Now build and run as normal
   > ./case.build
   > ./case.submit

.. note:: We did this by editing the "namelist_defaults_datm.xml" which will change the settings for ALL future ``I1850Clm50BgcSpinup`` cases you run. You could also do this by editing the path in the resulting streams text files in the CaseDocs directory, and then create a "user\_" streams file with the correct path. This would change the streams file JUST for this case. The steps do it this way are:

::

   > ./preview_namelists
   > cp CaseDocs/datm.streams.txt.CPLHIST3HrWx.Precip            user_datm.streams.txt.CPLHIST3HrWx.Precip
   > cp CaseDocs/datm.streams.txt.CPLHIST3HrWx.Solar             user_datm.streams.txt.CPLHIST3HrWx.Solar
   > cp CaseDocs/datm.streams.txt.CPLHIST3HrWx.nonSolarNonPrecip user_datm.streams.txt.CPLHIST3HrWx.nonSolarNonPrecip
   # Change the <fieldInfo> field <filePath> to point to the correct directory i.e.: /glade/home/achive/$USER/$DATM_CPLHIST_CASE/cpl/hist
   > $EDITOR user_datm.streams.txt.CPLHIST3HrWx.*
   > ./preview_namelists
   # Then make sure the CaseDocs/datm.streams.txt.CPLHIST3HrWx.* files have the correct path
