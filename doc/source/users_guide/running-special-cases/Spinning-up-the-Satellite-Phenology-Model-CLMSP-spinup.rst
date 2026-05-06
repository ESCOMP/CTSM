.. include:: ../substitutions.rst

.. _spinning-up-sp:

===========================================
 Spinning up the Satellite Phenology Model
===========================================

To spin-up the CLM60SP model you generally need to run CLM60SP for a few decades starting from arbitrary initial conditions, the main goal being to ensure that the turbulent fluxes and soil water and temperature are stable (minimal trends). You then use the final restart file resulting from that simulation as initial conditions in other simulations.  Alternatively, you can also start from an initial file that is typically already provided for you as part of the selected compset. Generally, this will allow for shorter spinup times if your simulation configuration is similar to the one run to generate the default initial file.

The following steps illustrate how to setup and run a 50 year CLM60SP spinup from arbitrary initial conditions using the ``I2000Clm60SpCrujra`` compset and ``f09_t232`` spatial resolution.  From a checkout of the CLM code (choose your own case name):
::

   cd cime/scripts
   ./create_newcase --case Clm60SP_ctsm54030_1deg_CRUJRA2024_arbi_2000 --compset I2000Clm60SpCrujra --res f09_t232 --run-unsupported --project XX
   cd Clm60SP_ctsm54030_1deg_CRUJRA2024_arbi_2000/
   ./case.setup
   ./xmlchange CLM_FORCE_COLDSTART=on
   ./xmlchange RUN_STARTDATE=0001-01-01
   ./xmlchange DATM_YR_START=1991
   ./xmlchange DATM_YR_END=2000
   ./xmlchange DATM_YR_ALIGN=1
   ./xmlchange STOP_OPTION=nyears
   ./xmlchange STOP_N=50
   ./case.build
   ./case.submit

Setting ``CLM_FORCE_COLDSTART=on`` forces the model to use arbitrary initial conditions (see Section :numref:`Initialization` for a description of these initial conditions).  This will result in ``finidat=' '`` in the ``lnd_in`` namelist. Spinups are generally started at year 1 (``RUN_STARTDATE=0001-01-01``). Here we've chosen to loop over years 1991-2000 (``DATM_YR_START=1991``, ``DATM_YR_END=2000``) of the atmospheric forcing (10 years total), align model year 1 (``DATM_YR_ALIGN=1``) with the first year of atmospheric forcing, and run for 50 years (``STOP_OPTION=nyears``, ``STOP_N=50``).

The spinup stability script available in the CLM checkout at ``tools/contrib/SpinupStability_SP_v10.ncl`` can be used to assess the stability or equilibrium of key model variables. Key settings in that script for this example simulation are
::

   caseid = "Clm60SP_ctsm54030_1deg_CRUJRA2024_fini_2000"
   subper = 10

The ``subper`` setting tells the script how many years of atmospheric forcing were repeated.

:numref:`Figure CLM60SP spinup plot for arbitrary initial conditions` shows spinup behavior for this simulation. Variables are plotted every 10 years, hence five points are plotted.  These include FSH (sensible heat flux), EFLX_LH_TOT (latent heat flux), FPSN (photosynthesis), H2OSOI (soil water at layer 8 which is about 1 meter), TSOI (soil temperature at layer 10 which is about 3 meters), and TWS (total water storage). The speed at which these variables reach a specified equilibrium state (denoted by falling within the dotted lines) varies by variable, TWS generally takes the longest to equilibrium.  The plot in the lower left denotes the percent of land area that is not in TWS equilibrium. The equilibrium thresholds are fairly arbitrary for the SP configuration and can be chosen by the user.  The current settings are
::

  glob_thresh_fsh = 0.02         ; global threshold for FSH equilibrium (delta W m-2 / yr)
  glob_thresh_lh = 0.02          ; global threshold for EFLX_LH_TOT equilibrium (delta W m-2 / yr)
  glob_thresh_gpp = 0.02         ; global threshold for FPSN equilibrium (delta PgC / yr)
  glob_thresh_tws = 0.001        ; global threshold for TWS equilibrium (delta m / yr)
  glob_thresh_h2osoi = 0.01      ; global threshold for H2OSOI equilibrium (delta mm mm-3 / yr)
  glob_thresh_tsoi = 0.02        ; global threshold for TSOI equilibrium (delta K / yr)
  glob_thresh_area = 3.0         ; global threshold percent area with TWS disequilibrium gt 0.01 m

.. _Figure CLM60SP spinup plot for arbitrary initial conditions:

.. figure:: Clm60SP_ctsm54030_1deg_CRUJRA2024_fini_2000_SP_Spinup.png

 SP spinup plot for arbitrary initial conditions. Variables examined are FSH (sensible heat flux), EFLX_LH_TOT (latent heat flux), GPP (photosynthesis), TWS (total water storage), H2OSOI (volumetric soil water in layer 8) and TSOI (soil temperature in layer 10). Generated using ``tools/contrib/SpinupStability_SP_v10.ncl``.

You can also start from a default initial file that is provided as part of the selected compset. The following steps illustrate how to setup and run a 50 year CLM60SP spinup from default initial conditions again using the ``I2000Clm60SpCrujra`` compset and ``f09_t232`` spatial resolution.  From a checkout of the CLM code (choose your own case name):
::

   cd cime/scripts
   ./create_newcase --case Clm60SP_ctsm54030_1deg_CRUJRA2024_fini_2000 --compset I2000Clm60SpCrujra --res f09_t232 --run-unsupported --project XX
   cd Clm60SP_ctsm54030_1deg_CRUJRA2024_fini_2000/
   ./case.setup
   echo "use_init_interp = .true" >> user_nl_clm
   ./xmlchange RUN_STARTDATE=0001-01-01
   ./xmlchange DATM_YR_START=1991
   ./xmlchange DATM_YR_END=2000
   ./xmlchange DATM_YR_ALIGN=1
   ./xmlchange STOP_OPTION=nyears
   ./xmlchange STOP_N=50
   ./case.build
   ./case.submit

The difference from the previous simulation is that we don't set ``CLM_FORCE_COLDSTART=on`` so that the model uses the default provided initial conditions.  In this case, setting ``use_init_interp = .true`` is required because the model configuration used is slightly different from that used to generate the initial file.

:numref:`Figure CLM60SP spinup plot for default initial conditions` shows spinup behavior for this simulation. Here we can see that equilbrium is reached much sooner because the default initial file is from a spinup where the model configuration was very similar to this one.

.. _Figure CLM60SP spinup plot for default initial conditions:

.. figure:: Clm60SP_ctsm54030_1deg_CRUJRA2024_fini_2000_SP_Spinup.png

 SP spinup plot for default initial conditions. Variables examined are FSH (sensible heat flux), EFLX_LH_TOT (latent heat flux), GPP (photosynthesis), TWS (total water storage), H2OSOI (volumetric soil water in layer 8) and TSOI (soil temperature in layer 10). Generated using ``tools/contrib/SpinupStability_SP_v10.ncl``.
