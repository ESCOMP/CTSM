.. _running-with-moar-data:

========================
 Running with MOAR data
========================

Because it takes so long to spinup the CN model (as we just saw previously), if you are doing fully coupled simulations with active atmosphere and ocean, you will want to do the spinup portion of this "offline". 
So instead of doing expensive fully coupled simulations for the spinup duration, you run CLM in a very cheap "I" compset using atmospheric forcing from a shorter fully coupled simulation (or a simulation run previously by someone else).

In this example we will use the ``I1850Clm50BgcSpinup compset`` to setup CLM to run with atmospheric forcing from a previous fully coupled simulation with data that is already stored on disk on cheyenne. 
There are several simulations that have high frequency data for which we can do this. You can also do this on a machine other than cheyenne, but would need to download the data from the Earth System Grid and change the datapath similar to `Example 4-11 <CLM-URL>`_.

Example: Simulation with MOAR Data on cheyenne
-------------------------------------------------------------
::

   > cd scripts
   > ./create_newcase -case MOARforce1850 -res f19_g17_gl4 -compset I1850Clm50BgcSpinup -mach cheyenne_intel
   > cd MOARforce1850
   # The following sets the casename to point to for atm forcing (you could also use an editor)
   > ./xmlchange DATM_CPL_CASE=b40.1850.track1.1deg.006a
   # The following sets the align year and years to run over for atm forcing 
   #  (you could also use an editor)
   > ./xmlchange DATM_CPL_YR_ALIGN=1,DATM_CPL_YR_START=960,DATM_CPL_YR_END=1030
   > ./cesm_setup
   # Now build and run as normal
   > ./MOARforce1850.build
   > ./MOARforce1850.submit


.. warning: Because of bug 1339 (see the `$CTSMROOT/doc/KnownBugs <CLM-URL>`_ file on this) you can't run with 83 or more years of forcing. If you do need to run with more years of forcing, you'll need to address the issue as outlined in the `$CTSMROOT/doc/KnownBugs <CLM-URL>`_ file.
