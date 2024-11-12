.. include:: ../substitutions.rst

.. _running-with-anomaly-forcing:

==============================
 Running with anomaly forcing
==============================
Because performing fully coupled climate simulations is computationally expensive, an alternate method of running land-only simulations forced by future climate projections was developed for CTSM called 'anomaly forcing'.  The anomaly forcing method uses a previously completed fully coupled simulation to create monthly anomalies, relative to the present day, of near-surface atmospheric states and fluxes.  These anomalies, representing the evolution of future climate projections, are applied to a repeating cycle of present day atmospheric forcing data, either as an additive (for states) or multiplicative (for fluxes) quantity.  Thus, high-frequency variability is obtained from the present day atmospheric forcing data, while the long-term evolution of the climate is determined by the anomaly forcing dataset.

To enable anomaly forcing in a CTSM simulation, the following namelist variable can be added to the user\_nl\_datm file:

  anomaly\_forcing = 'Anomaly.Forcing.Precip','Anomaly.Forcing.Temperature','Anomaly.Forcing.Pressure','Anomaly.Forcing.Humidity','Anomaly.Forcing.Uwind','Anomaly.Forcing.Vwind','Anomaly.Forcing.Shortwave','Anomaly.Forcing.Longwave'

Any combination or subset of forcing variables can be used, e.g. to modify only a single atmospheric forcing variable, one could use:

  anomaly\_forcing = 'Anomaly.Forcing.Temperature'

which will only adjust the temperature (TBOT).

After the namelist has been created, the run directory will be populated with files such as these:

  datm.streams.txt.Anomaly.Forcing.Temperature

which will contain the location of the default anomaly forcing datasets.  To use alternative data, copy these files to the case directory with the 'user\_' prefix, and modify the 'user\_*' files accordingly, e.g.:

  user\_datm.streams.txt.Anomaly.Forcing.Temperature

    For example, one could use the user\_datm.streams.txt.Anomaly.Forcing.* files to point to these SSP-derived anomaly forcing datasets:

    /glade/p/cesmdata/cseg/inputdata/atm/datm7/anomaly\_forcing/CMIP6-SSP3-7.0

    af.huss.cesm2.SSP3-7.0.2015-2100\_c20200329.nc
    af.pr.cesm2.SSP3-7.0.2015-2100\_c20200329.nc
    af.ps.cesm2.SSP3-7.0.2015-2100\_c20200329.nc
    af.rlds.cesm2.SSP3-7.0.2015-2100\_c20200329.nc
    af.rsds.cesm2.SSP3-7.0.2015-2100\_c20200329.nc
    af.tas.cesm2.SSP3-7.0.2015-2100\_c20200329.nc

Users may wish to also update files such as the landuse\_timeseries and aerosol and Ndepostion files to correspond to the appropriate SSP.

For single point simulations, the global anomaly forcing files can be used, but the map_algo namelist variable should be appended with nearest neighbor values for each of the anomaly forcing fields, e.g.

    mapalgo = 'nn','nn','nn','nn','nn','nn','nn','nn','nn','nn','nn','nn','nn' (the number of 'nn' values will depend on the number of original streams plus the number of anomaly forcing streams)

The cycling of the present-day (base) climate is controlled through the DATM\_YR\_START and DATM\_YR\_END variables in env\_run.xml.
