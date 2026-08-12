.. include:: ../substitutions.rst

.. _running-with-anomaly-forcing:

==============================
 Running with anomaly forcing
==============================
Because performing fully coupled climate simulations is computationally expensive, an alternate method of running land-only simulations forced by future climate projections was developed for CTSM called "anomaly forcing." The anomaly forcing method uses a previously-completed, fully-coupled simulation to create monthly anomalies, relative to the present day, of near-surface atmospheric states and fluxes. These anomalies, representing the evolution of future climate projections, are applied to a repeating cycle of present day atmospheric forcing data, either as an additive (for states) or multiplicative (for fluxes) quantity. Thus, high-frequency variability is obtained from the present day atmospheric forcing data, while the long-term evolution of the climate is determined by the anomaly forcing dataset.

Anomaly climate forcings are automatically enabled for ``ISSP`` compsets (e.g., ``ISSP585``). After the namelist has been created, ``CaseDocs/datm.streams.xml`` will have an entry like this pointing to the anomaly forcing file being used:

::

  <stream_info name="Anomaly.Forcing.cmip6.ssp585">
     ...
     <datafiles>
        <file>/glade/campaign/cesm/cesmdata/inputdata/atm/datm7/anomaly_forcing/CMIP6-SSP5-8.5/af.allvars.CESM.SSP5-8.5.2015-2100_c20220628.nc</file>
     </datafiles>
     <datavars>
        <var>huss  Sa_shum_af</var>
        <var>pr    Faxa_prec_af</var>
        <var>ps    Sa_pbot_af</var>
        <var>rlds  Faxa_lwdn_af</var>
        <var>rsds  Faxa_swdn_af</var>
        <var>tas   Sa_tbot_af</var>
        <var>uas   Sa_u_af</var>
        <var>vas   Sa_v_af</var>
     </datavars>
     ...
  </stream_info>

To use alternative data, add a ``user_nl_datm_streams`` namelist file to your case with contents like so:

::

  Anomaly.Forcing.cmip6.ssp585:datafiles = /path/to/your/datafile
  Anomaly.Forcing.cmip6.ssp585:meshfile = /path/to/meshfile/for/your/datafile

  ! List of Data types to use
  ! Remove the variables you do NOT want to include in the anomaly forcing:
  !     pr is precipitation
  !     tas is temperature
  !     huss is humidity
  !     uas and vas are U and V winds
  !     rsds is solare
  !     rlds is LW down
  Anomaly.Forcing.cmip6.ssp585:datavars = pr    Faxa_prec_af, \
                                          tas   Sa_tbot_af, \
                                          ps    Sa_pbot_af, \
                                          huss  Sa_shum_af, \
                                          uas   Sa_u_af, \
                                          vas   Sa_v_af, \
                                          rsds  Faxa_swdn_af, \
                                          rlds  Faxa_lwdn_af

To instead disable anomaly forcing in ``ISSP`` compsets, the following can be added to the ``user_nl_datm`` file:

::

  anomaly_forcing = 'none'

Note that other inputs are also set automatically for ``ISSP`` compsets, including CO2 (``co2tseries``), ozone (``preso3``), N deposition (``presndep``), and aerosols (``presaero``).

The first and last years over which the present-day (base) climate should cycle are set through the ``DATM_YR_START`` and ``DATM_YR_END`` XML variables.
