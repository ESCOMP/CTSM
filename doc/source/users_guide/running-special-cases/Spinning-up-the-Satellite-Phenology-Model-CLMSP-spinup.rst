.. include:: ../substitutions.rst

.. _spinning-up-sp:

===========================================
 Spinning up the Satellite Phenology Model
===========================================

To spin-up the CLMSP model you merely need to run CLMSP for about 50 simulation years starting from arbitrary initial conditions. You then use the final restart file for initial conditions in other simulations. Because this is a straight forward operation we will NOT give the details on how to do that here, but leave it as an exercise for the reader. See the Example :numref:`eg-final-clmbgc-spinup` as an example of doing this as the last step for CLMCN.

You can also start from a default initial file that is setup as part of the selected compset. :numref:`Figure SP spinup plot for 1850` shows spinup behavior for an 1850 SP case that loops over one year of coupler history output for atmospheric forcing (generated from the fully coupled model), initialized with an initial file generated from a GSWP3 atmospheric forcing case. Note that it takes less than 10 years for state variables such as FSH (sensible heat flux), EFLX_LH_TOT (latent heat flux), GPP (photosynthesis), H2OSOI (soil water), and TSOI (soil temperature) to reach a specified equilibrium state (denoted by the dotted lines) due to the different atmospheric forcing. TWS (total water storage) may take a bit longer.

.. _Figure SP spinup plot for 1850:

.. figure:: image1.png

 SP spinup plot for year 1850. Variables examined are FSH (sensible heat flux), EFLX_LH_TOT (latent heat flux), GPP (photosynthesis), TWS (total water storage), H2OSOI (volumetric soil water in layer 8) and TSOI (soil temperature in layer 10). Generated using .../tools/contrib/SpinupStability_SP.ncl.

:numref:`Figure SP spinup plot for 2000 CO2` shows spinup behavior for the same case but also changes CO2 to present-day conditions (379ppmv). Again, it takes about 10 years to reach equilibrium.

.. _Figure SP spinup plot for 2000 CO2:

.. figure:: image2.png

 SP spinup plot for year 2000 CO2. Variables examined are FSH (sensible heat flux), EFLX_LH_TOT (latent heat flux), GPP (photosynthesis), TWS (total water storage), H2OSOI (volumetric soil water in layer 8) and TSOI (soil temperature in layer 10). Generated using .../tools/contrib/SpinupStability_SP.ncl.
