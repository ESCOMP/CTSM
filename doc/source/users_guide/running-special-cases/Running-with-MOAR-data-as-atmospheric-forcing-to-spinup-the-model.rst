.. include:: ../substitutions.rst

.. _running-with-moar-data:

========================
 Running with MOAR data
========================

.. warning::
    These instructions are outdated and will not work. This page will be either updated or removed as part of the CTSM6/CLM3 release.

Because it takes so long to spinup the CN model (as we just saw previously), if you are doing fully coupled simulations with active atmosphere and ocean, you will want to do the spinup portion of this "offline". So instead of doing expensive fully coupled simulations for the spinup duration, you run CLM in a very cheap "I" compset using atmospheric forcing from a shorter fully coupled simulation (or a simulation run previously by someone else).

In this example we will use the ``I1850Clm50BgcSpinup compset`` to setup CLM to run with atmospheric forcing from a previous fully coupled simulation with data that is already stored on disk on Cheyenne. There are several simulations that have high frequency data for which we can do this. You can also do this on a machine other than Cheyenne, but would need to download the data from the Earth System Grid and change the datapath similar to Example :numref:`eg-sim-data-from-prev-sim`.

Example: Simulation with MOAR Data on derecho
-------------------------------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case MOARforce1850 -res f19_g17_gl4 -compset I1850Clm50BgcSpinup
   > cd MOARforce1850
   # The following sets the casename to point to for atm forcing (you could also use an editor)
   > ./xmlchange DATM_CPLHIST_CASE=b40.1850.track1.1deg.006a
   # The following sets the align year and years to run over for atm forcing
   #  (you could also use an editor)
   > ./xmlchange DATM_YR_ALIGN=1,DATM_YR_START=960,DATM_YR_END=1030
   > ./case.setup
   # Now build and run as normal
   > ./case.build
   > ./case.submit

