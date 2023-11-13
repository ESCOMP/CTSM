.. _rst_Land-only Mode:

Land-Only Mode
================

In land-only mode (uncoupled to an atmospheric model), the atmospheric forcing required by CLM (:numref:`Table Atmospheric input to land model`) is supplied by observed datasets. The standard forcing provided with the model is a 110-year (1901-2010) dataset provided by the Global Soil Wetness Project (GSWP3; NEED A REFERENCE). The GSWP3 dataset has a spatial resolution of 0.5° X 0.5° and a temporal resolution of three hours.

An alternative forcing dataset is also available, CRUNCEP, a 110-year (1901-2010) dataset (CRUNCEP; :ref:`Viovy 2011 <Viovy2011>`) that is a combination of two existing datasets; the CRU TS3.2 0.5° X 0.5° monthly data covering the period 1901 to 2002 (:ref:`Mitchell and Jones 2005 <MitchellJones2005>`) and the NCEP reanalysis 2.5° X 2.5° 6-hourly data covering the period 1948 to 2010. The CRUNCEP dataset has been used to force CLM for studies of vegetation growth, evapotranspiration, and gross primary production (:ref:`Mao et al. 2012 <Maoetal2012>`, :ref:`Mao et al. 2013 <Maoetal2013>`, :ref:`Shi et al. 2013 <Shietal2013>`) and for the TRENDY (trends in net land-atmosphere carbon exchange over the period 1980-2010) project (:ref:`Piao et al. 2012 <Piaoetal2012>`). Version 7 is available here (:ref:`Viovy 2011 <Viovy2011>`).

Here, the GSWP3 dataset, which does not include data for particular fields over oceans, lakes, and Antarctica is modified. This missing data is filled with :ref:`Qian et al. (2006) <Qianetal2006>` data from 1948 that is interpolated by the data atmosphere model to the 0.5° GSWP3 grid. This allows the model to be run over Antarctica and ensures data is available along coastlines regardless of model resolution.

The forcing data is ingested into a data atmosphere model in three "streams"; precipitation (:math:`P`) (mm s\ :sup:`-1`), solar radiation (:math:`S_{atm}` ) (W m\ :sup:`-2`), and four other fields [atmospheric pressure :math:`P_{atm}` (Pa), atmospheric specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`), atmospheric temperature :math:`T_{atm}` (K), and atmospheric wind :math:`W_{atm}` (m s\ :sup:`-1`)]. These are separate streams because they are handled differently according to the type of field. In the GSWP3 dataset, the precipitation stream is provided at three hour intervals and the data atmosphere model prescribes the same precipitation rate for each model time step within the three hour period. The four fields that are grouped together in another stream (pressure, humidity, temperature, and wind) are provided at three hour intervals and the data atmosphere model linearly interpolates these fields to the time step of the model.

The total solar radiation is also provided at three hour intervals. The data is fit to the model time step using a diurnal function that depends on the cosine of the solar zenith angle :math:`\mu` to provide a smoother diurnal cycle of solar radiation and to ensure that all of the solar radiation supplied by the three-hourly forcing data is actually used. The solar radiation at model time step :math:`t_{M}` is

.. math::
   :label: 31.1

   \begin{array}{lr}
   S_{atm} \left(t_{M} \right)=\frac{\frac{\Delta t_{FD} }{\Delta t_{M} } S_{atm} \left(t_{FD} \right)\mu \left(t_{M} \right)}{\sum _{i=1}^{\frac{\Delta t_{FD} }{\Delta t_{M} } }\mu \left(t_{M_{i} } \right) } & \qquad {\rm for\; }\mu \left(t_{M} \right)>0.001 \\
   S_{atm} \left(t_{M} \right)=0 & \qquad {\rm for\; }\mu \left(t_{M} \right)\le 0.001
   \end{array}

where :math:`\Delta t_{FD}` is the time step of the forcing data (3 hours :math:`\times` 3600 seconds hour\ :sup:`-1` = 10800 seconds), :math:`\Delta t_{M}` is the model time step (seconds), :math:`S_{atm} \left(t_{FD} \right)` is the three-hourly solar radiation from the forcing data (W m\ :sup:`-2`), and :math:`\mu \left(t_{M} \right)` is the cosine of the solar zenith angle at model time step :math:`t_{M}` (section :numref:`Solar Zenith Angle`). The term in the denominator of equation :eq:`31.1` is the sum of the cosine of the solar zenith angle for each model time step falling within the three hour period. For numerical purposes, :math:`\mu \left(t_{M_{i} } \right)\ge 0.001`.

The total incident solar radiation :math:`S_{atm}` at the model time step :math:`t_{M}` is then split into near-infrared and visible radiation and partitioned into direct and diffuse according to factors derived from one year's worth of hourly CAM output from CAM version cam3\_5\_55 as

.. math::
   :label: 31.2

   S_{atm} \, \downarrow _{vis}^{\mu } =R_{vis} \left(\alpha S_{atm} \right)

.. math::
   :label: 31.3

   S_{atm} \, \downarrow _{nir}^{\mu } =R_{nir} \left[\left(1-\alpha \right)S_{atm} \right]

.. math::
   :label: 31.4

   S_{atm} \, \downarrow _{vis} =\left(1-R_{vis} \right)\left(\alpha S_{atm} \right)

.. math::
   :label: 31.5

   S_{atm} \, \downarrow _{nir} =\left(1-R_{nir} \right)\left[\left(1-\alpha \right)S_{atm} \right].

where :math:`\alpha`, the ratio of visible to total incident solar radiation, is assumed to be

.. math::
   :label: 31.6

   \alpha =\frac{S_{atm} \, \downarrow _{vis}^{\mu } +S_{atm} \, \downarrow _{vis}^{} }{S_{atm} } =0.5.

The ratio of direct to total incident radiation in the visible :math:`R_{vis}` is

.. math::
   :label: 31.7

   R_{vis} =a_{0} +a_{1} \times \alpha S_{atm} +a_{2} \times \left(\alpha S_{atm} \right)^{2} +a_{3} \times \left(\alpha S_{atm} \right)^{3} \qquad 0.01\le R_{vis} \le 0.99

and in the near-infrared :math:`R_{nir}`  is

.. math::
   :label: 31.8

   R_{nir} =b_{0} +b_{1} \times \left(1-\alpha \right)S_{atm} +b_{2} \times \left[\left(1-\alpha \right)S_{atm} \right]^{2} +b_{3} \times \left[\left(1-\alpha \right)S_{atm} \right]^{3} \qquad 0.01\le R_{nir} \le 0.99

where :math:`a_{0} =0.17639,\, a_{1} =0.00380,\, a_{2} =-9.0039\times 10^{-6},\, a_{3} =8.1351\times 10^{-9}` and :math:`b_{0} =0.29548,b_{1} =0.00504,b_{2} =-1.4957\times 10^{-5},b_{3} =1.4881\times 10^{-8}` are coefficients from polynomial fits to the CAM data.

The additional atmospheric forcing variables required by :numref:`Table Atmospheric input to land model` are derived as follows. The atmospheric reference height :math:`z'_{atm}` (m) is set to 30 m. The directional wind components are derived as :math:`u_{atm} =v_{atm} ={W_{atm} \mathord{\left/ {\vphantom {W_{atm} \sqrt{2} }} \right.} \sqrt{2} }`. The potential temperature :math:`\overline{\theta _{atm} }` (K) is set to the atmospheric temperature :math:`T_{atm}`. The atmospheric longwave radiation :math:`L_{atm} \, \downarrow` (W m\ :sup:`-2`) is derived from the atmospheric vapor pressure :math:`e_{atm}` and temperature :math:`T_{atm}` (:ref:`Idso 1981<Idso1981>`) as

.. math::
   :label: 31.9

   L_{atm} \, \downarrow =\left[0.70+5.95\times 10^{-5} \times 0.01e_{atm} \exp \left(\frac{1500}{T_{atm} } \right)\right]\sigma T_{atm}^{4}

where

.. math::
   :label: 31.10

   e_{atm} =\frac{P_{atm} q_{atm} }{0.622+0.378q_{atm} }

and :math:`\sigma` is the Stefan-Boltzmann constant (W m\ :sup:`-2` K\ :sup:`-4`) (:numref:`Table Physical constants`). The fraction of precipitation :math:`P` (mm s\ :sup:`-1`) falling as rain and/or snow is

.. math::
   :label: 31.11

   q_{rain} =P\left(f_{P} \right),

.. math::
   :label: 31.12

   q_{snow} =P\left(1-f_{P} \right)

where

.. math::
   :label: 31.13

   f_{P} =0<0.5\left(T_{atm} -T_{f} \right)<1.

The aerosol deposition rates :math:`D_{sp}` (14 rates as described in :numref:`Table Atmospheric input to land model`) are provided by a time-varying, globally-gridded aerosol deposition file developed by :ref:`Lamarque et al. (2010) <Lamarqueetal2010>`.

If the user wishes to provide atmospheric forcing data from another source, the data format outlined above will need to be followed with the following exceptions. The data atmosphere model will accept a user-supplied relative humidity :math:`RH` (%) and derive specific humidity :math:`q_{atm}` (kg kg\ :sup:`-1`) from

.. math::
   :label: 31.14

   q_{atm} =\frac{0.622e_{atm} }{P_{atm} -0.378e_{atm} }

where the atmospheric vapor pressure :math:`e_{atm}` (Pa) is derived from the water (:math:`T_{atm} >T_{f}` ) or ice (:math:`T_{atm} \le T_{f}` ) saturation vapor pressure :math:`e_{sat}^{T_{atm} }` as :math:`e_{atm} =\frac{RH}{100} e_{sat}^{T_{atm} }` where :math:`T_{f}` is the freezing temperature of water (K) (:numref:`Table Physical constants`), and :math:`P_{atm}` is the pressure at height :math:`z_{atm}` (Pa). The data atmosphere model will also accept a user-supplied dew point temperature :math:`T_{dew}` (K) and derive specific humidity :math:`q_{atm}` from

.. math::
   :label: 31.15

   q_{atm} = \frac{0.622e_{sat}^{T_{dew} } }{P_{atm} -0.378e_{sat}^{T_{dew} } } .

Here, :math:`e_{sat}^{T}`, the saturation vapor pressure as a function of temperature, is derived from :ref:`Lowe's (1977) <Lowe1977>` polynomials. If not provided by the user, the atmospheric pressure :math:`P_{atm}` (Pa) is set equal to the standard atmospheric pressure :math:`P_{std} =101325` Pa, and surface pressure :math:`P_{srf}` (Pa) is set equal to\ :math:`P_{atm}`.

The user may provide the total direct and diffuse solar radiation, :math:`S_{atm} \, \downarrow ^{\mu }` and :math:`S_{atm} \, \downarrow`. These will be time-interpolated using the procedure described above and then each term equally apportioned into the visible and near-infrared wavebands (e.g., :math:`S_{atm} \, \downarrow _{vis}^{\mu } =0.5S_{atm} \, \downarrow ^{\mu }`, :math:`S_{atm} \, \downarrow _{nir}^{\mu } =0.5S_{atm} \, \downarrow ^{\mu }` ).

.. _Anomaly Forcing:

Anomaly Forcing
-----------------------------

The 'Anomaly Forcing' atmospheric forcing mode provides a means to drive CLM with projections of future climate conditions without the need for large, high-frequency datasets. From an existing climate simulation spanning both the historical and future time periods, a set of anomalies are created by removing a climatological seasonal cycle based on the end of the historical period from each year of the future time period of the simulation. These anomalies can then be applied to a repeating high-frequency forcing dataset of finite duration (e.g. 10 years). State and flux forcing variables are adjusted using additive and multiplicative anomalies, respectively:

.. math::
   :label: 31.16

   \begin{array}{lr}
   S^{'} = S + k_{anomaly} & \quad {\rm state \ variable} \\
   F^{'} = f \times k_{anomaly} & \quad {\rm flux \ variable}
   \end{array}

where :math:`S^{'}` is the adjusted atmospheric state variable, :math:`S` is the state variable from the high-frequency reference atmospheric forcing dataset, and :math:`k_{anomaly}` is an additive anomaly. Similarly, :math:`F^{'}` is the adjusted atmospheric flux variable, :math:`F` is the flux variable from the high-frequency reference atmospheric forcing dataset, and :math:`k_{anomaly}` is a multiplicative anomaly. State variables are temperature :math:`T_{atm}`, pressure :math:`P_{atm}`, humidity :math:`q_{atm}`, and wind :math:`W_{atm}`. Flux variables are precipitation :math:`P`, atmospheric shortwave radiation :math:`S_{atm} \, \downarrow`, and atmospheric longwave radiation :math:`L_{atm} \, \downarrow`.
