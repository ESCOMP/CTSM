.. _running-with-excess-ground-ice:

.. include:: ../substitutions.rst

===================================
 Running with excess ground ice
===================================


Excess ground ice can be toggled with ``use_excess_ice`` namelist option. By default this option is ``.false.``. When
``use_excess_ice`` is true, CTSM needs initial excess ice amount within soil layers to initialize. A second namelist option ``use_excess_ice_streams`` exists to control this process (``.false.`` is default). If ``.true.`` and ``use_excess_ice`` is ``.true.``,
initial conditions will be read from a data-stream file (default is based :ref:`on IPA map from 1997 <Brownetal1997>`).
This is useful, since in this way, a run with excess ground ice can be started from a restart or initial dataset, that does not include excess ground ice.
If the run is a continue-run, excess ice variables will **always** be expected on a restart file.

.. note:: Excess ice amount provided by the stream file is expressed in excess ice concentration (%) and does not have a vertical distribution. Each soil layer beneath 0.5 m (or maximum active layer depth from the previous year if it is greater) down to bedrock will receive the same concentration, but the ice mass will be scaled by the soil layer depth. Both naturally vegetated and crop columns get excess ice.


Since presence of excess ice within the soil significantly alters heat diffusion within it, when starting from initial conditions where excess ice was not present, an additional spinup is required.
Usually such spinup takes 100-150 years (depending on your climate) to completely equilibrate soil temperatures.



Example: Crop Simulation
------------------------------------
::

   > cd cime/scripts
   > ./create_newcase -case I1850Clm50BgcCrop_with_exice -res f19_g17_gl4 -compset I1850Clm50BgcCrop
   > cd I1850Clm50BgcCrop_with_exice

   > ./case.setup

   # turn on excess ice and its "stream" initialization
   > echo "use_excess_ice=.true." >> user_nl_clm
   > echo "use_excess_ice_streams=.true." >> user_nl_clm

   # Now build and run normally
   > ./case.build
   > ./case.submit
