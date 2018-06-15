.. _ptclm-examples:

.. include:: ../substitutions.rst

=========================
 Examples of using PTCLM
=========================

Now let's give a few more complex examples using some of the options we have discussed above.

In this first example, we'll demonstrate using a supported single point dataset, which then requires using the "nopointdata". We'll also demonstrate the compset option, "stdurbpt" and "caseidprefix" options.

Example: Running PTCLM for the Mexicocity supported single point dataset
------------------------------------------------------------------------
::

   > cd scripts/ccsm_utils/Tools/lnd/clm/PTCLM
   > ./PTCLM.py -m cheyenne_intel -s 1x1_mexicocityMEX -d $CSMDATA --nopointdata \
   --stdurbpt -c ICRUCLM45 --caseidprefix `pwd`/myPTCLMcases/site
   > cd myPTCLMcases/site_1x1_mexicocityMEX_I
   > ./cesm_setup
   # Now build and run normally
   > ./site_1x1_mexicocityMEX_I.build
   # Here we show running interactively
   > ./site_1x1_mexicocityMEX_I.run

Now, let's demonstrate using a different group list, doing a spinup, running with Qian global forcing data, but using tower years to set the years to run over. This uses the options: sitegroupname, useQIAN, and QIANtower_years.

Example: Running PTCLM for a spinup simulation with Qian data for tower years.
------------------------------------------------------------------------------
::
   
   > cd scripts/ccsm_utils/Tools/lnd/clm/PTCLM
   > ./PTCLM.py -m cheyenne_intel -s US-Ha1 -d $CSMDATA --sitegroupname AmeriFlux --useQIAN --QIAN_tower_yrs
   > cd ../../../../../US-Ha1_ICRUCLM45BGC_QIAN
   > ./cesm_setup
   # Now build and run normally
   > ./US-Ha1_ICRUCLM45BGC_QIAN.build
   # Here we show running interactively
   > ./US-Ha1_ICRUCLM45BGC_QIAN.run
   ```

Finally, let's demonstrate using a generic machine (which then requires the scratchroot option), using the global grid for PFT and soil types, and setting the run length to two months.

Example: Running PTCLM on a user-defined machine with global PFT and soil types dataset
---------------------------------------------------------------------------------------
::

   > cd scripts/ccsm_utils/Tools/lnd/clm/PTCLM
   # Note, see the the Section called Converting AmeriFlux Data for use by PTCLM with permission information
   # to use the US-UMB data.
   > ./PTCLM.py -m userdefined_intel -s US-UMB -d $CSMDATA --pftgrid --soilgrid \
   --scratchroot $HOME --run_n 2 --run_units nmonths
   > cd ../../../../../US-UMB_ICRUCLM45BGC
   # If userdefined is NOT set up for you Uncomment the following and set OS, NTASKS_PER_NODE, TMPDIR
   # > ./xmlchange OS=$OS,MAX_TASKS_PER_NODE=$NTASKS_PER_NODE,MPILIB=mpi-serial
   # > ./xmlchange RUNDIR=$TMPDIR/$USER/\$CASE/run,DIN_LOC_ROOT=$CSMDATA,COMPILER=intel
   # > ./xmlchange EXEROOT=$TMPDIR/$USER/\$CASE
   > ./cesm_setup
   # Now build
   > ./US-UMB_ICRUCLM45BGC.userdefined_intel.build
   # To get the files from the svn server...
   # First list the files from the streams text file
   > ../ccsm_utils/Tools/listfilesin_streams \
   -t $HOME/US-UMB_ICRUCLM45BGC/run/clm1PT.1x1pt_US-UMB.stream.txt -l \
   > Buildconf/datm.input_data_list
   # And now run the script to export data to your machine
   > ../ccsm_utils/Tools/check_input_data -export
   # Here we show running interactively
   > ./US-UMB_ICRUCLM45BGC.userdefined_intel.run

.. warning: Because of Bug 1364, when running this case as above we get a floating point error after reaching time-step 124 for the example exactly as above. Other machines or compilers probably won't have this problem. See the `$CTSMROOT/doc/KnownBugs <CLM-URL>`_ file for more information on this problem.
