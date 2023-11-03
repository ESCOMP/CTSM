.. _rst_Lake Model:

Lake Model
=============

The lake model, denoted the *Lake, Ice, Snow, and Sediment Simulator* (LISSS), is from :ref:`Subin et al. (2012a) <Subinetal2012a>`. It includes extensive modifications to the lake code of :ref:`Zeng et al. (2002) <Zengetal2002>` used in CLM versions 2 through 4, which utilized concepts from the lake models of :ref:`Bonan (1996) <Bonan1996>`, :ref:`Henderson-Sellers (1985) <Henderson-Sellers1985>`, :ref:`Henderson-Sellers (1986) <Henderson-Sellers1986>`, :ref:`Hostetler and Bartlein (1990) <HostetlerBartlein1990>`, and the coupled lake-atmosphere model of :ref:`Hostetler et al. (1993) <Hostetleretal1993>`, :ref:`Hostetler et al. (1993) <Hostetleretal1993>`. Lakes have spatially variable depth prescribed in the surface data (section :ref:`External Data Lake`); the surface data optionally includes lake optical extinction coeffient and horizontal fetch, currently only used for site simulations. Lake physics includes freezing and thawing in the lake body, resolved snow layers, and "soil" and bedrock layers below the lake body. Temperatures and ice fractions are simulated for :math:`N_{levlak} =10` layers (for global simulations) or :math:`N_{levlak} =25` (for site simulations) with discretization described in section :numref:`Vertical Discretization Lake`. Lake albedo is described in section :numref:`Surface Albedo Lake`. Lake surface fluxes (section :numref:`Surface Fluxes and Surface Temperature Lake`) generally follow the formulations for non-vegetated surfaces, including the calculations of aerodynamic resistances (section :numref:`Sensible and Latent Heat Fluxes for Non-Vegetated Surfaces`); however, the lake surface temperature :math:`T_{g}` (representing an infinitesimal interface layer between the top resolved lake layer and the atmosphere) is solved for simultaneously with the surface fluxes. After surface fluxes are evaluated, temperatures are solved simultaneously in the resolved snow layers (if present), the lake body, and the soil and bedrock, using the ground heat flux *G* as a top boundary condition. Snow, soil, and bedrock models generally follow the formulations for non-vegetated surfaces (Chapter :numref:`rst_Soil and Snow Temperatures`), with modifications described below.

.. _Vertical Discretization Lake:

Vertical Discretization
---------------------------

Currently, there is one lake modeled in each grid cell (with prescribed or assumed depth *d*, extinction coefficient :math:`\eta`, and fetch *f*), although this could be modified with changes to the CLM subgrid decomposition algorithm in future model versions. As currently implemented, the lake consists of 0-5 snow layers; water and ice layers (10 for global simulations and 25 for site simulations) comprising the "lake body;" 10 "soil" layers; and 5 bedrock layers. Each lake body layer has a fixed water mass (set by the nominal layer thickness and the liquid density), with frozen mass-fraction *I* a state variable. Resolved snow layers are present if the snow thickness :math:`z_{sno} \ge s_{\min }`, where *s*\ :sub:`min` = 4 cm by default, and is adjusted for model timesteps other than 1800 s in order to maintain numerical stability (section :numref:`Modifications to Snow Layer Logic Lake`). For global simulations with 10 body layers, the default (50 m lake) body layer thicknesses are given by: :math:`\Delta z_{i}` of 0.1, 1, 2, 3, 4, 5, 7, 7, 10.45, and 10.45 m, with node depths :math:`z_{i}` located at the center of each layer (i.e., 0.05, 0.6, 2.1, 4.6, 8.1, 12.6, 18.6, 25.6, 34.325, 44.775 m). For site simulations with 25 layers, the default thicknesses are (m): 0.1 for layer 1; 0.25 for layers 2-5; 0.5 for layers 6-9; 0.75 for layers 10-13; 2 for layers 14-15; 2.5 for layers 16-17; 3.5 for layers 18-21; and 5.225 for layers 22-25. For lakes with depth *d* :math:`\neq` 50 m and *d* :math:`\ge` 1 m, the top layer is kept at 10 cm and the other 9 layer thicknesses are adjusted to maintain fixed proportions. For lakes with *d* :math:`<` 1 m, all layers have equal thickness. Thicknesses of snow, soil, and bedrock layers follow the scheme used over non-vegetated surfaces (Chapter :numref:`rst_Soil and Snow Temperatures`), with modifications to the snow layer thickness rules to keep snow layers at least as thick as *s*\ :sub:`min` (section :numref:`Modifications to Snow Layer Logic Lake`).

.. _External Data Lake:

External Data
-----------------

As discussed in :ref:`Subin et al. (2012a, b) <Subinetal2012a>`, the Global Lake and Wetland Database (:ref:`Lehner and Doll 2004<LehnerDoll2004>`) is currently used to prescribe lake fraction in each land model grid cell, for a total of 2.3 million km\ :sup:`-2`. As in :ref:`Subin et al. (2012a, b) <Subinetal2012a>`, the :ref:`Kourzeneva et al. (2012)<Kourzenevaetal2012>` global gridded dataset is currently used to estimate a mean lake depth in each grid cell, based on interpolated compilations of geographic information.

.. _Surface Albedo Lake:

Surface Albedo
------------------

For direct radiation, the albedo *a* for lakes with ground temperature :math:`{T}_{g}` (K) above freezing is given by (:ref:`Pivovarov, 1972<Pivovarov1972>`)

.. math::
   :label: 12.1

   a=\frac{0.5}{\cos z+0.15}

where *z* is the zenith angle. For diffuse radiation, the expression in eq. is integrated over the full sky to yield *a* = 0.10.

For frozen lakes without resolved snow layers, the albedo at cold temperatures *a*\ :sub:`0` is 0.60 for visible and 0.40 for near infrared radiation. As the temperature at the ice surface, :math:`{T}_{g}`, approaches freezing [ :math:`{T}_{f}` (K) (:numref:`Table Physical Constants`)], the albedo is relaxed towards 0.10 based on :ref:`Mironov et al. (2010)<Mironovetal2010>`:

.. math::
   :label: 12.2

   a=a_{0} \left(1-x\right)+0.10x,x=\exp \left(-95\frac{T_{f} -T_{g} }{T_{f} } \right)

where *a* is restricted to be no less than that given in :eq:`12.1`.

For frozen lakes with resolved snow layers, the reflectance of the ice surface is fixed at *a*\ :sub:`0`, and the snow reflectance is calculated as over non-vegetated surfaces (Chapter :numref:`rst_Surface Albedos`). These two reflectances are combined to obtain the snow-fraction-weighted albedo as in over non-vegetated surfaces (Chapter :numref:`rst_Surface Albedos`).

.. _Surface Fluxes and Surface Temperature Lake:

Surface Fluxes and Surface Temperature
------------------------------------------

.. _Surface Properties Lake:

Surface Properties
^^^^^^^^^^^^^^^^^^^^^^^^

The fraction of shortwave radiation absorbed at the surface, :math:`\beta`, depends on the lake state. If resolved snow layers are present, then :math:`\beta` is set equal to the absorption fraction predicted by the snow-optics submodel (Chapter :numref:`rst_Surface Albedos`) for the top snow layer. Otherwise, :math:`\beta` is set equal to the near infrared fraction of the shortwave radiation reaching the surface simulated by the atmospheric model or atmospheric data model used for offline simulations (Chapter :numref:`rst_Land-only Mode`). The remainder of the shortwave radiation fraction (1 :math:`{-}` :math:`\beta`) is absorbed in the lake body or soil as described in section :numref:`Radiation Penetration`.

The surface roughnesses are functions of the lake state and atmospheric forcing. 

For unfrozen lakes (:math:`T_{g} > T_{f}`), :math:`z_{0m}` is given by (:ref:`Subin et al. (2012a) <Subinetal2012a>`)

.. math::
   :label: 12.3

   z_{0m} =\max \left(\frac{\alpha \nu }{u_{*} } ,C\frac{u_{*} ^{2} }{g} \right)

where :math:`\alpha` = 0.1, :math:`\nu` is the kinematic viscosity of air given below, *C* is the effective Charnock coefficient given below, :math:`u_{*}` is the friction velocity (m/s), and *g* is the acceleration of gravity (:numref:`Table Physical Constants`). The kinematic viscosity is given by

.. math::
   :label: 12.4 

   \nu =\nu _{0} \left(\frac{T_{g} }{T_{0} } \right)^{1.5} \frac{P_{0} }{P_{ref} }

where
:math:`\nu _{0} =1.51\times 10^{-5} {\textstyle\frac{{\rm m}^{{\rm 2}} }{{\rm s}}}`
, :math:`T_{0} ={\rm 293.15\; K}`,
:math:`P_{0} =1.013\times 10^{5} {\rm \; Pa}` , and
:math:`P_{ref}` is the pressure at the atmospheric reference height. The Charnock coefficient *C* is a function of the lake fetch *F* (m), given in the surface data or set to 25 times the lake depth *d* by default:

.. math::
   :label: 12.5 

   \begin{array}{l} {C=C_{\min } +(C_{\max } -C_{\min } )\exp \left\{-\min \left(A,B\right)\right\}} \\ {A={\left(\frac{Fg}{u_{*} ^{2} } \right)^{{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} } \mathord{\left/ {\vphantom {\left(\frac{Fg}{u_{*} ^{2} } \right)^{{1\mathord{\left/ {\vphantom {1 3}} \right.} 3} }  f_{c} }} \right.} f_{c} } } \\ {B=\varepsilon \frac{\sqrt{dg} }{u} } \end{array}

where *A* and *B* define the fetch- and depth-limitation, respectively;
:math:`C_{\min } =0.01` , :math:`C_{\max } =0.01`,
:math:`\varepsilon =1` , :math:`f_{c} =100` , and
*u* (m s\ :sup:`-1`) is the atmospheric forcing wind.

The scalar roughness lengths
(:math:`z_{0q}` for latent heat and :math:`z_{0h}` for sensible heat) are given by
(:ref:`Subin  et al. 2012a<Subinetal2012a>`)

.. math::
   :label: 12.5a

   \begin{array}{l} {R_{0} =(\frac{z_{0m} u_{*} }{\nu })^{0.5} ,} \\ {z_{0h} =z_{0m} \exp \left\{-\frac{k} {Pr} (4 R_{0} ^{0.5} -3.2) \right\},} \\ {z_{0q} =z_{0m} \exp \left\{-\frac{k} {Sc} (4 R_{0} ^{0.5} - 4.2) \right\}}\end{array}

where :math:`R_{0}` is the near-surface atmospheric roughness Reynolds number, :math:`k` is the von Karman constant (:numref:`Table Physical Constants`), :math:`Pr = 0.713` is the molecular Prandt number for air at neutral stability, :math:`Sc = 0.66` is the Schmidt number for water in air at neutral stability.
:math:`z_{0q}` and :math:`z_{0h}` are restricted to be no smaller than :math:`1 \times 10^{-10}`.

For frozen lakes ( :math:`T_{g} \le T_{f}` ) without resolved snow
layers ( :math:`snl = 0` ), :math:`z_{0m} =z_{0m_{ice}} =2.3\times 10^{-3} {\rm m}` (:ref:`Meier et al. (2022) <Meieretal2022>`).

For frozen lakes with resolved
snow layers ( :math:`snl > 0` ), the momentum roughness length is evaluated based on accumulated snow melt :math:`M_{a} {\rm }` (:ref:`Meier et al. (2022) <Meieretal2022>`). 
For :math:`M_{a} >=1\times 10^{-5}`

.. math::
   :label: 12.5b

   z_{0m} =\exp (b_{1} \tan ^{-1} \left[\frac{log_{10} (M_{a}) + 0.23)} {0.08}\right] + b_{4})\times 10^{-3}

where :math:`M_{a}` is accumulated snow melt (meters water equivalent), :math:`b_{1} =1.4` and :math:`b_{4} =-0.31`.
For :math:`M_{a} <1\times 10^{-5}`

.. math::
   :label: 12.5c

   z_{0m} =\exp (-b_{1} 0.5 \pi + b_{4})\times 10^{-3}

Accumulated snow melt :math:`M_{a}` at the current time step :math:`t` is defined as

.. math::
   :label: 12.5d   

   M ^{t}_{a} = M ^{t-1}_{a} - (q ^{t}_{sno} \Delta t + q ^{t}_{snowmelt} \Delta t)\times 10^{-3}

where :math:`M ^{t}_{a}` and :math:`M ^{t-1}_{a}` are the accumulated snowmelt at the current time step and previous time step, respectively (m), :math:`q ^{t}_{sno} \Delta t` is the freshly fallen snow (mm), and :math:`q ^{t}_{snowmelt} \Delta t` is the melted snow (mm).

For frozen lakes without and with resolved snow layers, an initial guess for the scalar roughness lengths is derived by assuming :math:`\theta_{*} = 0 {\rm }` (:ref:`Meier et al. (2022) <Meieretal2022>`)

.. math::
   :label: 12.5e

   z_{0h}=z_{0q}=\frac{70 \nu}{u_{*}}

where :math:`\nu=1.5 \times 10^{-5}` is the kinematic viscosity of air (m\ :sup:`2`  s\ :sup:`-1`), and
:math:`u_{*}` is the friction velocity in the atmospheric surface layer (m s\ :sup:`-1`).
Thereafter, the scalar roughness lengths are updated within the stability iteration described in section :numref:`Surface Flux Solution Lake` as

.. math::
   :label: 12.6

   z_{0h}=z_{0q}=\frac{70 \nu}{u_{*}} \exp (-\beta {u_{*}} ^{0.5} |{\theta_{*}}| ^{0.25} )

where :math:`\beta` = 7.2, and :math:`\theta_{*}` is the potential temperature scale (section :numref:`Surface Flux Solution Lake`).

.. _Surface Flux Solution Lake:

Surface Flux Solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conservation of energy at the lake surface requires

.. math::
   :label: 12.7

   \beta \vec{S}_{g} -\vec{L}_{g} -H_{g} -\lambda E_{g} -G=0

where :math:`\vec{S}_{g}` \ is the absorbed solar radiation in the lake,
:math:`\beta` is the fraction absorbed at the surface,
:math:`\vec{L}_{g}` \ is the net emitted longwave radiation (+ upwards),
:math:`H_{g}` \ is the sensible heat flux (+ upwards),
:math:`E_{g}` \ is the water vapor flux (+ upwards), and
*G* is the ground heat flux (+ downwards). All of these fluxes depend implicitly on the temperature at the lake surface :math:`{T}_{g}`. :math:`\lambda` converts :math:`E_{g}` to an energy flux based on

.. math::
   :label: 12.8

   \lambda =\left\{\begin{array}{l} {\lambda _{sub} \qquad T_{g} \le T_{f} } \\ {\lambda _{vap} \qquad T_{g} >T_{f} } \end{array}\right\}.

The sensible heat flux (W m\ :sup:`-2`) is

.. math::
   :label: 12.9

   H_{g} =-\rho _{atm} C_{p} \frac{\left(\theta _{atm} -T_{g} \right)}{r_{ah} }

where :math:`\rho _{atm}` is the density of moist air (kg m\ :sup:`-3`) (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`),
:math:`C_{p}` is the specific heat capacity of air (J kg\ :sup:`-1` K\ :sup:`-1`) (:numref:`Table Physical Constants`),
:math:`\theta _{atm}` is the atmospheric potential temperature (K) (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`),
:math:`T_{g}` is the lake surface temperature (K) (at an infinitesimal interface just above the top resolved model layer: snow, ice, or water), and
:math:`r_{ah}` is the aerodynamic resistance to sensible heat transfer (s m\ :sup:`-1`) (section :numref:`Monin-Obukhov Similarity Theory`).

The water vapor flux (kg m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 12.10

   E_{g} =-\frac{\rho _{atm} \left(q_{atm} -q_{sat}^{T_{g} } \right)}{r_{aw} }

where :math:`q_{atm}` is the atmospheric specific humidity (kg kg\ :sup:`-1`) (section :numref:`Atmospheric Coupling`),
:math:`q_{sat}^{T_{g} }` \ is the saturated specific humidity (kg kg\ :sup:`-1`) (section :numref:`Saturation Vapor Pressure`) at the lake surface temperature :math:`T_{g}`, and
:math:`r_{aw}` is the aerodynamic resistance to water vapor transfer (s m\ :sup:`-1`) (section :numref:`Monin-Obukhov Similarity Theory`).

The zonal and meridional momentum fluxes are

.. math::
   :label: 12.11

   \tau _{x} =-\rho _{atm} \frac{u_{atm} }{r_{atm} }

.. math::
   :label: 12.12

   \tau _{y} =-\rho _{atm} \frac{v_{atm} }{r_{atm} }

where :math:`u_{atm}` and :math:`v_{atm}` are the zonal and meridional atmospheric winds (m s\ :sup:`-1`) (section :numref:`Atmospheric Coupling`), and
:math:`r_{am}` is the aerodynamic resistance for momentum (s m\ :sup:`-1`) (section :numref:`Monin-Obukhov Similarity Theory`).

The heat flux into the lake surface :math:`G` (W m\ :sup:`-2`) is

.. math::
   :label: 12.13

   G=\frac{2\lambda _{T} }{\Delta z_{T} } \left(T_{g} -T_{T} \right)

where :math:`\lambda _{T}` is the thermal conductivity (W m\ :sup:`-1` K\ :sup:`-1`), :math:`\Delta z_{T}` is the thickness (m), and :math:`T_{T}` is the temperature (K) of the top resolved lake layer (snow, ice, or water). The top thermal conductivity :math:`\lambda _{T}` of unfrozen lakes ( :math:`T_{g} >T_{f}` ) includes conductivities due to molecular ( :math:`\lambda _{liq}` ) and eddy (:math:`\lambda _{K}` ) diffusivities (section :numref:`Eddy Diffusivity and Thermal Conductivities`), as evaluated in the top lake layer at the previous timestep, where :math:`\lambda _{liq}` is the thermal conductivity of water (:numref:`Table Physical Constants`). For frozen lakes without resolved snow layers, :math:`\lambda _{T} =\lambda _{ice}`. When resolved snow layers are present, :math:`\lambda _{T}` \ is calculated based on the water content, ice content, and thickness of the top snow layer, as for non-vegetated surfaces.

The absorbed solar radiation :math:`\vec{S}_{g}`  is

.. math::
   :label: 12.14

   \vec{S}_{g} =\sum _{\Lambda }S_{atm} \, \downarrow _{\Lambda }^{\mu } \left(1-\alpha _{g,\, \Lambda }^{\mu } \right) +S_{atm} \, \downarrow _{\Lambda } \left(1-\alpha _{g,\, \Lambda } \right)

where :math:`S_{atm} \, \downarrow _{\Lambda }^{\mu }` and :math:`S_{atm} \, \downarrow _{\Lambda }` are the incident direct beam and diffuse solar fluxes (W m\ :sup:`-2`) and :math:`\Lambda` denotes the visible (:math:`<` 0.7\ :math:`\mu {\rm m}`) and near-infrared (:math:`\ge` 0.7\ :math:`\mu {\rm m}`) wavebands (section :numref:`Atmospheric Coupling`), and :math:`\alpha _{g,\, \Lambda }^{\mu }` and :math:`\alpha _{g,\, \mu }` are the direct beam and diffuse lake albedos (section :numref:`Surface Albedo Lake`).

The net emitted longwave radiation is

.. math::
   :label: 12.15

   \vec{L}_{g} =L_{g} \, \uparrow -L_{atm} \, \downarrow

where :math:`L_{g} \, \uparrow` is the upward longwave radiation from the surface,
:math:`L_{atm} \, \downarrow` is the downward atmospheric longwave radiation (section :numref:`Atmospheric Coupling`). The upward longwave radiation from the surface is

.. math::
   :label: 12.16

   L\, \uparrow =\left(1-\varepsilon _{g} \right)L_{atm} \, \downarrow +\varepsilon _{g} \sigma \left(T_{g}^{n} \right)^{4} +4\varepsilon _{g} \sigma \left(T_{g}^{n} \right)^{3} \left(T_{g}^{n+1} -T_{g}^{n} \right)

where :math:`\varepsilon _{g} =0.97` is the lake surface emissivity,
:math:`\sigma` is the Stefan-Boltzmann constant (W m\ :sup:`-2` K\ :sup:`-4`) (:numref:`Table Physical Constants`), and
:math:`T_{g}^{n+1} -T_{g}^{n}` is the difference in lake surface temperature between Newton-Raphson iterations (see below).

The sensible heat :math:`H_{g}`, the water vapor flux :math:`E_{g}` through its dependence on the saturated specific humidity, the net longwave radiation :math:`\vec{L}_{g}`, and the ground heat flux :math:`G`, all depend on the lake surface temperature :math:`T_{g}`. Newton-Raphson iteration is applied to solve for :math:`T_{g}` and the surface fluxes as

.. math::
   :label: 12.17

   \Delta T_{g} =\frac{\beta \overrightarrow{S}_{g} -\overrightarrow{L}_{g} -H_{g} -\lambda E_{g} -G}{\frac{\partial \overrightarrow{L}_{g} }{\partial T_{g} } +\frac{\partial H_{g} }{\partial T_{g} } +\frac{\partial \lambda E_{g} }{\partial T_{g} } +\frac{\partial G}{\partial T_{g} } }

where :math:`\Delta T_{g} =T_{g}^{n+1} -T_{g}^{n}` and the subscript "n" indicates the iteration. Therefore, the surface temperature :math:`T_{g}^{n+1}` can be written as

.. math::
   :label: 12.18

   T_{g}^{n+1} =\frac{\beta \overrightarrow{S}_{g} -\overrightarrow{L}_{g} -H_{g} -\lambda E_{g} -G+T_{g}^{n} \left(\frac{\partial \overrightarrow{L}_{g} }{\partial T_{g} } +\frac{\partial H_{g} }{\partial T_{g} } +\frac{\partial \lambda E_{g} }{\partial T_{g} } +\frac{\partial G}{\partial T_{g} } \right)}{\frac{\partial \overrightarrow{L}_{g} }{\partial T_{g} } +\frac{\partial H_{g} }{\partial T_{g} } +\frac{\partial \lambda E_{g} }{\partial T_{g} } +\frac{\partial G}{\partial T_{g} } }

where the partial derivatives are

.. math::
   :label: 12.19

   \frac{\partial \overrightarrow{L}_{g} }{\partial T_{g} } =4\varepsilon _{g} \sigma \left(T_{g}^{n} \right)^{3} ,

.. math::
   :label: 12.20

   \frac{\partial H_{g} }{\partial T_{g} } =\frac{\rho _{atm} C_{p} }{r_{ah} } ,

.. math::
   :label: 12.21

   \frac{\partial \lambda E_{g} }{\partial T_{g} } =\frac{\lambda \rho _{atm} }{r_{aw} } \frac{dq_{sat}^{T_{g} } }{dT_{g} } ,

.. math::
   :label: 12.22

   \frac{\partial G}{\partial T_{g} } =\frac{2\lambda _{T} }{\Delta z_{T} } .

The fluxes of momentum, sensible heat, and water vapor are solved for simultaneously with lake surface temperature as follows. To begin, :math:`z_{0m}` and the scalar roughness lengths are set as described in section :numref:`Surface Properties Lake`.

#. An initial guess for the wind speed :math:`V_{a}` including the convective velocity :math:`U_{c}` is obtained from :eq:`5.24` assuming an initial convective velocity :math:`U_{c} =0` m s\ :sup:`-1` for stable conditions (:math:`\theta _{v,\, atm} -\theta _{v,\, s} \ge 0` as evaluated from :eq:`5.50`) and :math:`U_{c} =0.5` for unstable conditions (:math:`\theta _{v,\, atm} -\theta _{v,\, s} <0`).

#. An initial guess for the Monin-Obukhov length :math:`L` is obtained from the bulk Richardson number using :eq:`5.46` and :eq:`5.48`.

#. The following system of equations is iterated four times:

#. Heat of vaporization / sublimation :math:`\lambda` \ (:eq:`12.8`)

#. Thermal conductivity :math:`\lambda _{T}` \ (above)

#. Friction velocity :math:`u_{*}` (:eq:`5.32`, :eq:`5.33`, :eq:`5.34`, :eq:`5.35`)

#. Potential temperature scale :math:`\theta _{*}` (:eq:`5.37`, :eq:`5.38`, :eq:`5.39`, :eq:`5.40`)

#. Humidity scale :math:`q_{*}` (:eq:`5.41`, :eq:`5.42`, :eq:`5.43`, :eq:`5.44`)

#. Aerodynamic resistances :math:`r_{am}`, :math:`r_{ah}`, and :math:`r_{aw}` (:eq:`5.55`, :eq:`5.56`, :eq:`5.57`)

#. Lake surface temperature :math:`T_{g}^{n+1}` (:eq:`12.18`)

#. Heat of vaporization / sublimation :math:`\lambda` (:eq:`12.8`)

#. Sensible heat flux :math:`H_{g}` is updated for :math:`T_{g}^{n+1}` (:eq:`12.9`)

#. Water vapor flux :math:`E_{g}` is updated for :math:`T_{g}^{n+1}` as

   .. math::
      :label: 12.23

      E_{g} =-\frac{\rho _{atm} }{r_{aw} } \left[q_{atm} -q_{sat}^{T_{g} } -\frac{\partial q_{sat}^{T_{g} } }{\partial T_{g} } \left(T_{g}^{n+1} -T_{g}^{n} \right)\right]

where the last term on the right side of equation :eq:`12.23` is the change in saturated specific humidity due to the change in :math:`T_{g}` between iterations.

#. Saturated specific humidity :math:`q_{sat}^{T_{g} }` and its derivative :math:`\frac{dq_{sat}^{T_{g} } }{dT_{g} }` are updated for :math:`T_{g}^{n+1}` (section :numref:`Monin-Obukhov Similarity Theory`).

#. Virtual potential temperature scale :math:`\theta _{v*}` (:eq:`5.17`)

#. Wind speed including the convective velocity, :math:`V_{a}` (:eq:`5.24`)

#. Monin-Obukhov length :math:`L` (:eq:`5.49`)

#. Roughness lengths (section :numref:`Surface Properties Lake`).

Once the four iterations for lake surface temperature have been yielded a tentative solution :math:`T_{g} ^{{'} }`, several restrictions are imposed in order to maintain consistency with the top lake model layer temperature :math:`T_{T}` \ (:ref:`Subin et al. (2012a) <Subinetal2012a>`).

.. math::
   :label: 12.24

   \begin{array}{l} {{\rm 1)\; }T_{T} \le T_{f} <T_{g} ^{{'} } \Rightarrow T_{g} =T_{f} ,} \\ {{\rm 2)\; }T_{T} >T_{g} ^{{'} } >T_{m} \Rightarrow T_{g} =T_{T} ,} \\ {{\rm 3)\; }T_{m} >T_{g} ^{{'} } >T_{T} >T_{f} \Rightarrow T_{g} =T_{T} } \end{array}

where :math:`T_{m}` \ is the temperature of maximum liquid water density, 3.85°C (:ref:`Hostetler and Bartlein (1990) <HostetlerBartlein1990>`). The first condition requires that, if there is any snow or ice present, the surface temperature is restricted to be less than or equal to freezing. The second and third conditions maintain convective stability in the top lake layer.

If equation :eq:`12.24` is applied, the turbulent fluxes :math:`H_{g}` and :math:`E_{g}` are re-evaluated. The emitted longwave radiation and the momentum fluxes are re-evaluated in any case. The final ground heat flux :math:`G` is calculated from the residual of the energy balance (equation :eq:`12.7`) in order to precisely conserve energy. This ground heat flux is taken as a prescribed flux boundary condition for the lake temperature solution (section :numref:`Boundary Conditions Lake`). A check is included at each timestep to insure that energy balance is obeyed to within 0.1 W m\ :sup:`-2` (see :numref:`Energy Conservation Lake`).

.. _Lake Temperature:

Lake Temperature
--------------------

.. _Introduction Lake:

Introduction
^^^^^^^^^^^^^^^^^^

The (optional-) snow, lake body (water and/or ice), soil, and bedrock system is unified for the lake temperature solution. The governing equation, similar to that for the snow-soil-bedrock system for vegetated land units (Chapter :numref:`rst_Soil and Snow Temperatures`), is

.. math::
   :label: 12.25

   \tilde{c}_{v} \frac{\partial T}{\partial t} =\frac{\partial }{\partial z} \left(\tau \frac{\partial T}{\partial z} \right)-\frac{d\phi }{dz}

where :math:`\tilde{c}_{v}` is the volumetric heat capacity (J m\ :sup:`-3` K\ :sup:`-1`),
:math:`t` is time (s),
*T* is the temperature (K),
:math:`\tau` is the thermal conductivity (W m\ :sup:`-1` K\ :sup:`-1`), and
:math:`\phi` is the solar radiation (W m\ :sup:`-2`) penetrating to depth *z* (m). The system is discretized into *N* layers, where

.. math::
   :label: 12.26

   N=n_{sno} +N_{levlak} +N_{levgrnd}  ,

:math:`n_{sno}` is the number of actively modeled snow layers at the current timestep (Chapter :numref:`rst_Snow Hydrology`), and
:math:`N_{levgrnd}` \ is as for vegetated land units (Chapter :numref:`rst_Soil and Snow Temperatures`). Energy is conserved as

.. math::
   :label: 12.27

   \frac{d}{dt} \sum _{j=1}^{N}\left[\tilde{c}_{v,j} (t)\left(T_{j} -T_{f} \right)+L_{j} (t)\right] \Delta z_{j} =G+\left(1-\beta \right)\vec{S}_{g}

where :math:`\tilde{c}_{v,j} (t)`\ is the volumetric heat capacity of the *j*\ th layer (section :numref:`Radiation Penetration`), :math:`L_{j} (t)`\ is the latent heat of fusion per unit volume of the *j*\ th layer (proportional to the mass of liquid water present), and the right-hand side represents the net influx of energy to the lake system. Note that :math:`\tilde{c}_{v,j} (t)` can only change due to phase change (except for changing snow layer mass, which, apart from energy required to melt snow, represents an untracked energy flux in the land model, along with advected energy associated with water flows in general), and this is restricted to occur at :math:`T_{j} =T_{f}` \ in the snow-lake-soil system, allowing eq. to be precisely enforced and justifying the exclusion of :math:`c_{v,j}` from the time derivative in eq..

.. _Overview of Changes from CLM4 2:

Overview of Changes from CLM4
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Thermal conductivities include additional eddy diffusivity, beyond the :ref:`Hostetler and Bartlein (1990)<HostetlerBartlein1990>` formulation, due to unresolved processes (:ref:`Fang and Stefan 1996<FangStefan1996>`; :ref:`Subin et al. (2012a) <Subinetal2012a>`). Lake water is now allowed to freeze by an arbitrary fraction for each layer, which releases latent heat and changes thermal properties. Convective mixing occurs for all lakes, even if frozen. Soil and bedrock are included beneath the lake. The full snow model is used if the snow thickness exceeds a threshold; if there are resolved snow layers, radiation transfer is predicted by the snow-optics submodel (Chapter :numref:`rst_Surface Albedos`), and the remaining radiation penetrating the bottom snow layer is absorbed in the top layer of lake ice; conversely, if there are no snow layers, the solar radiation penetrating the bottom lake layer is absorbed in the top soil layer. The lakes have variable depth, and all physics is assumed valid for arbitrary depth, except for a depth-dependent enhanced mixing (section :numref:`Eddy Diffusivity and Thermal Conductivities`). Finally, a previous sign error in the calculation of eddy diffusivity (specifically, the Brunt-Väisälä frequency term; eq. ) was corrected.

.. _Boundary Conditions Lake:

Boundary Conditions
^^^^^^^^^^^^^^^^^^^^^^^^^

The top boundary condition, imposed at the top modeled layer :math:`i=j_{top}`, where :math:`j_{top} =-n_{sno} +1`, is the downwards surface flux *G* defined by the energy flux residual during the surface temperature solution (section :numref:`Boundary Conditions Lake`). The bottom boundary condition, imposed at :math:`i=N_{levlak} +N_{levgrnd}`, is zero flux. The 2-m windspeed :math:`u_{2}` \ (m s\ :sup:`-1`) is used in the calculation of eddy diffusivity:

.. math::
   :label: 12.28

   u_{2} =\frac{u_{*} }{k} \ln \left(\frac{2}{z_{0m} } \right)\ge 0.1.

where :math:`u_{*}` \ is the friction velocity calculated in section :numref:`Boundary Conditions Lake` and
*k* is the von Karman constant (:numref:`Table Physical Constants`).

.. _Eddy Diffusivity and Thermal Conductivities:

Eddy Diffusivity and Thermal Conductivities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The total eddy diffusivity :math:`K_{W}`  (m\ :sup:`2` s\ :sup:`-1`) for liquid water in the lake body is given by (:ref:`Subin et al. (2012a) <Subinetal2012a>`)

.. math::
   :label: 12.29

   K_{W} = m_{d} \left(\kappa _{e} +K_{ed} +\kappa _{m} \right)

where :math:`\kappa _{e}` is due to wind-driven eddies (:ref:`Hostetler and Bartlein (1990)<HostetlerBartlein1990>`),
:math:`K_{ed}` is a modest enhanced diffusivity intended to represent unresolved mixing processes (:ref:`Fang and Stefan 1996<FangStefan1996>`),
:math:`\kappa _{m} =\frac{\lambda _{liq} }{c_{liq} \rho _{liq} }` \ is the molecular diffusivity of water (given by the ratio of its thermal conductivity (W m\ :sup:`-1` K\ :sup:`-1`) to the product of its heat capacity (J kg\ :sup:`-1` K\ :sup:`-1`) and density (kg m\ :sup:`-3`), values given in :numref:`Table Physical Constants`), and
:math:`m_{d}` (unitless) is a factor which increases the overall diffusivity for large lakes, intended to represent 3-dimensional mixing processes such as caused by horizontal temperature gradients. As currently implemented,

.. math::
   :label: 12.30

   m_{d} =\left\{\begin{array}{l} {1,\qquad d<25{\rm m}} \\ {10,\qquad d\ge 25{\rm m}} \end{array}\right\}

where *d* is the lake depth.

The wind-driven eddy diffusion coefficient :math:`\kappa _{e,\, i}` (m\ :sup:`2` s\ :sup:`-1`) for layers :math:`1\le i\le N_{levlak}`  is

.. math::
   :label: 12.31

   \kappa _{e,\, i} =\left\{\begin{array}{l} {\frac{kw^{*} z_{i} }{P_{0} \left(1+37Ri^{2} \right)} \exp \left(-k^{*} z_{i} \right)\qquad T_{g} >T_{f} } \\ {0\qquad T_{g} \le T_{f} } \end{array}\right\}

where :math:`P_{0} =1` is the neutral value of the turbulent Prandtl number, :math:`z_{i}` is the node depth (m), the surface friction velocity (m s\ :sup:`-1`) is :math:`w^{*} =0.0012u_{2}`, and :math:`k^{*}` varies with latitude :math:`\phi` as :math:`k^{*} =6.6u_{2}^{-1.84} \sqrt{\left|\sin \phi \right|}`. For the bottom layer, :math:`\kappa _{e,\, N_{levlak} } =\kappa _{e,N_{levlak} -1\, }`. As in :ref:`Hostetler and Bartlein (1990)<HostetlerBartlein1990>`, the 2-m wind speed :math:`u_{2}` (m s\ :sup:`-1`) (eq. ) is used to evaluate :math:`w^{*}` and :math:`k^{*}` rather than the 10-m wind used by :ref:`Henderson-Sellers (1985) <Henderson-Sellers1985>`.

The Richardson number is

.. math::
   :label: 12.32

   R_{i} =\frac{-1+\sqrt{1+\frac{40N^{2} k^{2} z_{i}^{2} }{w^{*^{2} } \exp \left(-2k^{*} z_{i} \right)} } }{20}

where

.. math::
   :label: 12.33

   N^{2} =\frac{g}{\rho _{i} } \frac{\partial \rho }{\partial z}

and :math:`g` is the acceleration due to gravity (m s\ :sup:`-2`) (:numref:`Table Physical Constants`), :math:`\rho _{i}` is the density of water (kg m\ :sup:`-3`), and :math:`\frac{\partial \rho }{\partial z}` is approximated as :math:`\frac{\rho _{i+1} -\rho _{i} }{z_{i+1} -z_{i} }`. Note that because here, *z* is increasing downwards (unlike in :ref:`Hostetler and Bartlein (1990)<HostetlerBartlein1990>`), eq. contains no negative sign; this is a correction from CLM4. The density of water is (:ref:`Hostetler and Bartlein (1990)<HostetlerBartlein1990>`)

.. math::
   :label: 12.34

   \rho _{i} =1000\left(1-1.9549\times 10^{-5} \left|T_{i} -277\right|^{1.68} \right).

The enhanced diffusivity :math:`K_{ed}` is given by (:ref:`Fang and Stefan 1996<FangStefan1996>`)

.. math::
   :label: 12.35

   K_{ed} =1.04\times 10^{-8} \left(N^{2} \right)^{-0.43} ,N^{2} \ge 7.5\times 10^{-5} {\rm s}^{2}

where :math:`N^{2}` \ is calculated as in eq. except for the minimum value imposed in.

The thermal conductivity for the liquid water portion of lake body layer *i*, :math:`\tau _{liq,i}` (W m\ :sup:`-1` K\ :sup:`-1`) is given by

.. math::
   :label: 12.36

   \tau _{liq,i} =K_{W} c_{liq} \rho _{liq}  .

The thermal conductivity of the ice portion of lake body layer *i*, :math:`\tau _{ice,eff}` \ (W m\ :sup:`-1` K\ :sup:`-1`), is constant among layers, and is given by

.. math::
   :label: 12.37

   \tau _{ice,eff} =\tau _{ice} \frac{\rho _{ice} }{\rho _{liq} }

where :math:`\tau _{ice}` \ (:numref:`Table Physical Constants`) is the nominal thermal conductivity of ice: :math:`\tau _{ice,eff}` \ is adjusted for the fact that the nominal model layer thicknesses remain constant even while the physical ice thickness exceeds the water thickness.

The overall thermal conductivity :math:`\tau _{i}` for layer *i* with ice mass-fraction :math:`I_{i}` is the harmonic mean of the liquid and water fractions, assuming that they will be physically vertically stacked, and is given by

.. math::
   :label: 12.38

   \tau _{i} =\frac{\tau _{ice,eff} \tau _{liq,i} }{\tau _{liq,i} I_{i} +\tau _{ice} \left(1-I_{i} \right)}  .

The thermal conductivity of snow, soil, and bedrock layers above and below the lake, respectively, are computed identically to those for vegetated land units (Chapter :numref:`rst_Soil and Snow Temperatures`), except for the adjustment of thermal conductivity for frost heave or excess ice (:ref:`Subin et al., 2012a, Supporting Information<Subinetal2012a>`).

.. _Radiation Penetration:

Radiation Penetration
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If there are no resolved snow layers, the surface absorption fraction :math:`\beta` is set according to the near-infrared fraction simulated by the atmospheric model. This is apportioned to the surface energy budget (section :numref:`Surface Properties Lake`), and thus no additional radiation is absorbed in the top :math:`z_{a}` (currently 0.6 m) of unfrozen lakes, for which the light extinction coefficient :math:`\eta` (m\ :sup:`-1`) varies between lake columns (eq. ). For frozen lakes (:math:`T_{g} \le T_{f}` ), the remaining :math:`\left(1-\beta \right)\vec{S}_{g}` fraction of surface absorbed radiation that is not apportioned to the surface energy budget is absorbed in the top lake body layer. This is a simplification, as lake ice is partially transparent. If there are resolved snow layers, then the snow optics submodel (Chapter :numref:`rst_Surface Albedos`) is used to calculate the snow layer absorption (except for the absorption predicted for the top layer by the snow optics submodel, which is assigned to the surface energy budget), with the remainder penetrating snow layers absorbed in the top lake body ice layer.

For unfrozen lakes, the solar radiation remaining at depth :math:`z>z_{a}` in the lake body is given by

.. math::
   :label: 12.39

   \phi =\left(1-\beta \vec{S}_{g} \right)\exp \left\{-\eta \left(z-z_{a} \right)\right\} .

For all lake body layers, the flux absorbed by the layer *i*,
:math:`\phi _{i}`  , is

.. math::
   :label: 12.40

   \phi _{i} =\left(1-\beta \vec{S}_{g} \right)\left[\exp \left\{-\eta \left(z_{i} -\frac{\Delta z_{i} }{2} -z_{a} \right)\right\}-\exp \left\{-\eta \left(z_{i} +\frac{\Delta z_{i} }{2} -z_{a} \right)\right\}\right] .

The argument of each exponent is constrained to be non-negative (so :math:`\phi _{i}` = 0 for layers contained within :math:`{z}_{a}`). The remaining flux exiting the bottom of layer :math:`i=N_{levlak}` is absorbed in the top soil layer.

The light extinction coefficient :math:`\eta` (m\ :sup:`-1`), if not provided as external data, is a function of depth *d* (m) (:ref:`Subin et al. (2012a) <Subinetal2012a>`):

.. math::
   :label: 12.41

   \eta =1.1925d^{-0.424}  .

.. _Heat Capacities Lake:

Heat Capacities
^^^^^^^^^^^^^^^^^^^^^

The vertically-integrated heat capacity for each lake layer, :math:`\text{c}_{v,i}` (J m\ :sup:`-2`) is determined by the mass-weighted average over the heat capacities for the water and ice fractions:

.. math::
   :label: 12.42

   c_{v,i} =\Delta z_{i} \rho _{liq} \left[c_{liq} \left(1-I_{i} \right)+c_{ice} I_{i} \right] .

Note that the density of water is used for both ice and water fractions, as the thickness of the layer is fixed.

The total heat capacity :math:`c_{v,i}` for each soil, snow, and bedrock layer (J m\ :sup:`-2`) is determined as for vegetated land units (Chapter :numref:`rst_Soil and Snow Temperatures`), as the sum of the heat capacities for the water, ice, and mineral constituents.

.. _Crank-Nicholson Solution Lake:

Crank-Nicholson Solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The solution method for thermal diffusion is similar to that used for soil (Chapter :numref:`rst_Soil and Snow Temperatures`), except that the lake body layers are sandwiched between the snow and soil layers (section :numref:`Introduction Lake`), and radiation flux is absorbed throughout the lake layers. Before solution, layer temperatures :math:`T_{i}` (K), thermal conductivities :math:`\tau _{i}` (W m\ :sup:`-1` K\ :sup:`-1`), heat capacities :math:`c_{v,i}` (J m\ :sup:`-2`), and layer and interface depths from all components are transformed into a uniform set of vectors with length :math:`N=n_{sno} +N_{levlak} +N_{levgrnd}` and consistent units to simplify the solution. Thermal conductivities at layer interfaces are calculated as the harmonic mean of the conductivities of the neighboring layers:

.. math::
   :label: 12.43

   \lambda _{i} =\frac{\tau _{i} \tau _{i+1} \left(z_{i+1} -z_{i} \right)}{\tau _{i} \left(z_{i+1} -\hat{z}_{i} \right)+\tau _{i+1} \left(\hat{z}_{i} -z_{i} \right)}  ,

where :math:`\lambda _{i}` is the conductivity at the interface between layer *i* and layer *i +* 1,
:math:`z_{i}` is the depth of the node of layer *i*, and
:math:`\hat{z}_{i}` is the depth of the interface below layer *i*. Care is taken at the boundaries between snow and lake and between lake and soil. The governing equation is discretized for each layer as

.. math::
   :label: 12.44

   \frac{c_{v,i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)=F_{i-1} -F_{i} +\phi _{i}

where superscripts *n* + 1 and *n* denote values at the end and beginning of the timestep :math:`\Delta t`, respectively,
:math:`F_{i}` (W m\ :sup:`-2`) is the downward heat flux at the bottom of layer *i*, and
:math:`\phi _{i}` is the solar radiation absorbed in layer *i*.

Eq. is solved using the semi-implicit Crank-Nicholson Method, resulting in a tridiagonal system of equations:

.. math::
   :label: 12.45

   \begin{array}{l} {r_{i} =a_{i} T_{i-1}^{n+1} +b_{i} T_{i}^{n+1} +cT_{i+1}^{n+1} ,} \\ {a_{i} =-0.5\frac{\Delta t}{c_{v,i} } \frac{\partial F_{i-1} }{\partial T_{i-1}^{n} } ,} \\ {b_{i} =1+0.5\frac{\Delta t}{c_{v,i} } \left(\frac{\partial F_{i-1} }{\partial T_{i-1}^{n} } +\frac{\partial F_{i} }{\partial T_{i}^{n} } \right),} \\ {c_{i} =-0.5\frac{\Delta t}{c_{v,i} } \frac{\partial F_{i} }{\partial T_{i}^{n} } ,} \\ {r_{i} =T_{i}^{n} +0.5\frac{\Delta t}{c_{v,i} } \left(F_{i-1} -F_{i} \right)+\frac{\Delta t}{c_{v,i} } \phi _{i} .} \end{array}

The fluxes :math:`F_{i}` are defined as follows: for the top layer, :math:`F_{j_{top} -1} =2G;a_{j_{top} } =0`, where *G* is defined as in section :numref:`Boundary Conditions Lake` (the factor of 2 merely cancels out the Crank-Nicholson 0.5 in the equation for :math:`r_{j_{top} }` ). For the bottom layer, :math:`F_{N_{levlak} +N_{levgrnd} } =0`. For all other layers:

.. math::
   :label: 12.46

   F_{i} =\lambda _{i} \frac{T_{i} ^{n} -T_{i+1}^{n} }{z_{n+1} -z_{n} }  .

.. _Phase Change Lake:

Phase Change
^^^^^^^^^^^^^^^^^^

Phase change in the lake, snow, and soil is done similarly to that done for the soil and snow for vegetated land units (Chapter :numref:`rst_Soil and Snow Temperatures`), except without the allowance for freezing point depression in soil underlying lakes. After the heat diffusion is calculated, phase change occurs in a given layer if the temperature is below freezing and liquid water remains, or if the temperature is above freezing and ice remains.

If melting occurs, the available energy for melting, :math:`Q_{avail}` (J m\ :sup:`-2`), is computed as

.. math::
   :label: 12.47

   Q_{avail} =\left(T_{i} -T_{f} \right)c_{v,i}

where :math:`T_{i}` is the temperature of the layer after thermal diffusion (section :numref:`Crank-Nicholson Solution Lake`), and
:math:`c_{v,i}` \ is as calculated in section :numref:`Heat Capacities Lake`. The mass of melt in the layer *M* (kg m\ :sup:`-2`) is given by

.. math::
   :label: 12.48

   M=\min \left\{M_{ice} ,\frac{Q_{avail} }{H_{fus} } \right\}

where :math:`H_{fus}` (J kg\ :sup:`-1`) is the latent heat of fusion of water (:numref:`Table Physical Constants`), and
:math:`M_{ice}` is the mass of ice in the layer: :math:`I_{i} \rho _{liq} \Delta z_{i}` for a lake body layer, or simply the soil / snow ice content state variable (:math:`w_{ice}` ) for a soil / snow layer. The heat remainder, :math:`Q_{rem}` \ is given by

.. math::
   :label: 12.49

   Q_{rem} =Q_{avail} -MH_{fus}  .

Finally, the mass of ice in the layer :math:`M_{ice}` is adjusted downwards by :math:`M`, and the temperature :math:`T_{i}` of the layer is adjusted to

.. math::
   :label: 12.50

   T_{i} =T_{f} +\frac{Q_{rem} }{c'_{v,i} }

where :math:`c'_{v,i} =c_{v,i} +M\left(c_{liq} -c_{ice} \right)`.

If freezing occurs, :math:`Q_{avail}` is again given by but will be negative. The melt :math:`M`, also negative, is given by

.. math::
   :label: 12.51

   M=\max \left\{-M_{liq} ,\frac{Q_{avail} }{H_{fus} } \right\}

where :math:`M_{liq}` is the mass of water in the layer: :math:`\left(1-I_{i} \right)\rho _{liq} \Delta z_{i}` for a lake body layer, or the soil / snow water content state variable (:math:`w_{liq}` ). The heat remainder :math:`Q_{rem}` is given by eq. and will be negative or zero. Finally, :math:`M_{liq}` is adjusted downwards by :math:`-M` and the temperature is reset according to eq..

In the presence of nonzero snow water :math:`W_{sno}` without resolved snow layers over an unfrozen top lake layer, the available energy in the top lake layer :math:`\left(T_{1} -T_{f} \right)c_{v,1}` is used to melt the snow. Similar to above, :math:`W_{sno}` is either completely melted and the remainder of heat returned to the top lake layer, or the available heat is exhausted and the top lake layer is set to freezing. The snow thickness is adjusted downwards in proportion to the amount of melt, maintaining constant density.

.. _Convection Lake:

Convection
^^^^^^^^^^^^^^^^

Convective mixing is based on :ref:`Hostetler et al.'s (1993, 1994)<Hostetleretal1993>` coupled lake-atmosphere model, adjusting the lake temperature after diffusion and phase change to maintain a stable density profile. Unfrozen lakes overturn when :math:`\rho _{i} >\rho _{i+1}`, in which case the layer thickness weighted average temperature for layers 1 to :math:`i+1` is applied to layers 1 to :math:`i+1` and the densities are updated. This scheme is applied iteratively to layers :math:`1\le i<N_{levlak} -1`. Unstable profiles occurring at the bottom of the lake (i.e., between layers :math:`i=N_{levlak} -1` and :math:`i=N_{levlak}` ) are treated separately (:ref:`Subin et al. (2012a) <Subinetal2012a>`), as occasionally these can be induced by heat expelled from the sediments (not present in the original :ref:`Hostetler et al. (1994)<Hostetleretal1994>` model). Mixing proceeds from the bottom upward in this case (i.e., first mixing layers :math:`i=N_{levlak} -1` and :math:`i=N_{levlak}`, then checking :math:`i=N_{levlak} -2` and :math:`i=N_{levlak} -1` and mixing down to :math:`i=N_{levlak}` if needed, and on to the top), so as not to mix in with warmer over-lying layers.\

For frozen lakes, this algorithm is generalized to conserve total enthalpy and ice content, and to maintain ice contiguous at the top of the lake. Thus, an additional mixing criterion is added: the presence of ice in a layer that is below a layer which is not completely frozen. When this occurs, these two lake layers and all those above mix. Total enthalpy *Q* is conserved as

.. math::
   :label: 12.52

   Q=\sum _{j=1}^{i+1}\Delta z_{j} \rho _{liq} \left(T_{j} -T_{f} \right)\left[\left(1-I_{j} \right)c_{liq} +I_{j} c_{ice} \right]  .

Once the average ice fraction :math:`I_{av}`  is calculated from

.. math::
   :label: 12.53

   \begin{array}{l} {I_{av} =\frac{\sum _{j=1}^{i+1}I_{j} \Delta z_{j}  }{Z_{i+1} } ,} \\ {Z_{i+1} =\sum _{j=1}^{i+1}\Delta z_{j}  ,} \end{array}

the temperatures are calculated. A separate temperature is calculated for the frozen (:math:`T_{froz}` ) and unfrozen (:math:`T_{unfr}` ) fractions of the mixed layers. If the total heat content *Q* is positive (e.g. some layers will be above freezing), then the extra heat is all assigned to the unfrozen layers, while the fully frozen layers are kept at freezing. Conversely, if :math:`Q < 0`, the heat deficit will all be assigned to the ice, and the liquid layers will be kept at freezing. For the layer that contains both ice and liquid (if present), a weighted average temperature will have to be calculated.

If :math:`Q > 0`, then :math:`T_{froz} =T_{f}`, and :math:`T_{unfr}` is given by

.. math::
   :label: 12.54

   T_{unfr} =\frac{Q}{\rho _{liq} Z_{i+1} \left[\left(1-I_{av} \right)c_{liq} \right]} +T_{f}  .

If :math:`Q < 0`, then :math:`T_{unfr} =T_{f}`, and :math:`T_{froz}` is given by

.. math::
   :label: 12.55

   T_{froz} =\frac{Q}{\rho _{liq} Z_{i+1} \left[I_{av} c_{ice} \right]} +T_{f}  .

The ice is lumped together at the top. For each lake layer *j* from 1 to *i* + 1, the ice fraction and temperature are set as follows, where :math:`Z_{j} =\sum _{m=1}^{j}\Delta z_{m}` :

#. If :math:`Z_{j} \le Z_{i+1} I_{av}`, then :math:`I_{j} =1` and :math:`T_{j} =T_{froz}`.

#. Otherwise, if :math:`Z_{j-1} <Z_{i+1} I_{av}`, then the layer will contain both ice and water. The ice fraction is given by :math:`I_{j} =\frac{Z_{i+1} I_{av} -Z_{j-1} }{\Delta z_{j} }`. The temperature is set to conserve the desired heat content that would be present if the layer could have two temperatures, and then dividing by the heat capacity of the layer to yield

   .. math::
      :label: 12.56

      T_{j} =\frac{T_{froz} I_{j} c_{ice} +T_{unfr} \left(1-I_{j} \right)c_{liq} }{I_{j} c_{ice} +\left(1-I_{j} \right)c_{liq} }  .

#. Otherwise, :math:`I_{j} =0` and :math:`T_{j} =T_{unfr}` .

.. _Energy Conservation Lake:

Energy Conservation
^^^^^^^^^^^^^^^^^^^^^^^^^^

To check energy conservation, the left-hand side of equation :eq:`12.27` is re-written to yield the total enthalpy of the lake system (J m\ :sup:`-2`) :math:`H_{tot}` :

.. math::
   :label: 12.57

   H_{tot} =\sum _{i=j_{top} }^{N_{levlak} +N_{levgrnd} }\left[c_{v,i} \left(T_{i} -T_{f} \right)+M_{liq,i} H_{fus} \right] -W_{sno,bulk} H_{fus}

where :math:`M_{liq,i}` is the water mass of the *i*\ th layer (similar to section :numref:`Phase Change Lake`), and :math:`W_{sno,bulk}` is the mass of snow-ice not present in resolved snow layers. This expression is evaluated once at the beginning and once at the end of the timestep (re-evaluating each :math:`c_{v,i}` ), and the change is compared with the net surface energy flux to yield the error flux :math:`E_{soi}` (W m\ :sup:`-2`):

.. math::
   :label: 12.58

   E_{soi} =\frac{\Delta H_{tot} }{\Delta t} -G-\sum _{i=j_{top} }^{N_{levlak} +N_{levgrnd} }\phi _{i}

If :math:`\left|E_{soi} \right|<0.1`\ W m\ :sup:`-2`, it is subtracted from the sensible heat flux and added to *G*. Otherwise, the model is aborted.

.. _Lake Hydrology:

Lake Hydrology
------------------

.. _Overview Lake Hydrology:

Overview
^^^^^^^^^^^^^^

Hydrology is done similarly to other impervious non-vegetated columns (e.g., glaciers) where snow layers may be resolved but infiltration into the permanent ground is not allowed. The water mass of lake columns is currently maintained constant, aside from overlying snow. The water budget is balanced with :math:`q_{rgwl}` (eq.; kg m\ :sup:`-2` s\ :sup:`-1`), a generalized runoff term for impervious land units that may be negative.

There are some modifications to the soil and snow parameterizations as compared with the soil in vegetated land units, or the snow overlying other impervious columns. The soil can freeze or thaw, with the allowance for frost heave (or the initialization of excess ice) (sections :numref:`Eddy Diffusivity and Thermal Conductivities` and :numref:`Phase Change Lake`), but no air-filled pore space is allowed in the soil. To preserve numerical stability in the lake model (which uses a slightly different surface flux algorithm than over other non-vegetated land units), two changes are made to the snow model. First, dew or frost is not allowed to be absorbed by a top snow layer which has become completely melted during the timestep. Second, because occasional instabilities occurred during model testing when the Courant–Friedrichs–Lewy (CFL) condition was violated, due to the explicit time-stepping integration of the surface flux solution, resolved snow layers must be a minimum of :math:`s_{\min }` = 4 cm thick rather than 1 cm when the default timestep of 1800 s is used.

.. _Water Balance Lake:

Water Balance
^^^^^^^^^^^^^^^^^^^

The total water balance of the system is given by

.. math::
   :label: 12.59

   \Delta W_{sno} +\sum _{i=1}^{n_{levsoi} }\left(\Delta w_{liq,i} +\Delta w_{ice,i} \right) =\left(q_{rain} +q_{sno} -E_{g} -q_{rgwl} -q_{snwcp,\, ice} \right)\Delta t

where :math:`W_{sno}` (kg m\ :sup:`-2`) is the total mass of snow (both liquid and ice, in resolved snow layers or bulk snow), :math:`w_{liq,i}` and :math:`w_{ice,i}` are the masses of water phases (kg m\ :sup:`-2`) in soil layer *i*, :math:`q_{rain}` and :math:`q_{sno}` are the precipitation forcing from the atmosphere (kg m\ :sup:`-2` s\ :sup:`-1`), :math:`q_{snwcp,\, ice}` is the ice runoff associated with snow-capping (below), :math:`E_{g}` is the ground evaporation (section :numref:`Surface Flux Solution Lake`), and :math:`n_{levsoi}` is the number of hydrologically active soil layers (as opposed to dry bedrock layers).

.. _Precipitation, Evaporation, and Runoff Lake:

Precipitation, Evaporation, and Runoff
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All precipitation reaches the ground, as there is no vegetated fraction. As for other land types, incident snowfall accumulates (with ice mass :math:`W_{sno}` and thickness :math:`z_{sno}` ) until its thickness exceeds a minimum thickness :math:`s_{\min }`, at which point a resolved snow layer is initiated, with water, ice, dissolved aerosol, snow-grain radius, etc., state variables tracked by the Snow Hydrology submodel (Chapter :numref:`rst_Snow Hydrology`). The density of fresh snow is assigned as for other land types (Chapter :numref:`rst_Snow Hydrology`). Solid precipitation is added immediately to the snow, while liquid precipitation is added to snow layers, if they exist, after accounting for dew, frost, and sublimation (below). If :math:`z_{sno}` exceeds :math:`s_{\min }` after solid precipitation is added but no snow layers are present, a new snow layer is initiated immediately, and then dew, frost, and sublimation are accounted for. Snow-capping is invoked if the snow depth :math:`z_{sno} >1000{\rm m}`, in which case additional precipitation and frost deposition is added to :math:`q_{snwcp,\, ice}`.

If there are resolved snow layers, the generalized "evaporation" :math:`E_{g}` (i.e., evaporation, dew, frost, and sublimation) is treated as over other land units, except that the allowed evaporation from the ground is unlimited (though the top snow layer cannot lose more water mass than it contains). If there are no resolved snow layers but :math:`W_{sno} >0` and :math:`E_{g} >0`, sublimation :math:`q_{sub,sno}` \ (kg m\ :sup:`-2` s\ :sup:`-1`) will be given by

.. math::
   :label: 12.60

   q_{sub,sno} =\min \left\{E_{g} ,\frac{W_{sno} }{\Delta t} \right\} .

If :math:`E_{g} <0,T_{g} \le T_{f}`, and there are no resolved snow layers or the top snow layer is not unfrozen, then the rate of frost production :math:`q_{frost} =\left|E_{g} \right|`. If :math:`E_{g} <0` but the top snow layer has completely thawed during the Phase Change step of the Lake Temperature solution (section :numref:`Phase Change Lake`), then frost (or dew) is not allowed to accumulate (:math:`q_{frost} =0`), to insure that the layer is eliminated by the Snow Hydrology (Chapter :numref:`rst_Snow Hydrology`) code. (If :math:`T_{g} >T_{f}`, then no snow is present (section :numref:`Surface Flux Solution Lake`), and evaporation or dew deposition is balanced by :math:`q_{rgwl}`.) The snowpack is updated for frost and sublimation:

.. math::
   :label: 12.61

   W_{sno} =W_{sno} +\Delta t\left(q_{frost} -q_{sub,sno} \right) .

If there are resolved snow layers, then this update occurs using the Snow Hydrology submodel (Chapter :numref:`rst_Snow Hydrology`). Otherwise, the snow ice mass is updated directly, and :math:`z_{sno}` is adjusted by the same proportion as the snow ice (i.e., maintaining the same density), unless there was no snow before adding the frost, in which case the density is assumed to be 250 kg m\ :sup:`-3`.

.. _Soil Hydrology Lake:

Soil Hydrology
^^^^^^^^^^^^^^^^^^^^

The combined water and ice soil volume fraction in a soil layer :math:`\theta _{i}` is given by

.. math::
   :label: 12.62

   \theta _{i} =\frac{1}{\Delta z_{i} } \left(\frac{w_{ice,i} }{\rho _{ice} } +\frac{w_{liq,i} }{\rho _{liq} } \right) .

If :math:`\theta _{i} <\theta _{sat,i}`, the pore volume fraction at saturation (as may occur when ice melts), then the liquid water mass is adjusted to

.. math::
   :label: 12.63

   w_{liq,i} =\left(\theta _{sat,i} \Delta z_{i} -\frac{w_{ice,i} }{\rho _{ice} } \right)\rho _{liq}  .

Otherwise, if excess ice is melting and :math:`w_{liq,i} >\theta _{sat,i} \rho _{liq} \Delta z_{i}`, then the water in the layer is reset to

.. math::
   :label: 12.64

   w_{liq,i} = \theta _{sat,i} \rho _{liq} \Delta z_{i}

This allows excess ice to be initialized (and begin to be lost only after the pore ice is melted, which is realistic if the excess ice is found in heterogeneous chunks) but irreversibly lost when melt occurs.

.. _Modifications to Snow Layer Logic Lake:

Modifications to Snow Layer Logic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A thickness difference :math:`z_{lsa} =s_{\min } -\tilde{s}_{\min }` adjusts the minimum resolved snow layer thickness for lake columns as compared to non-lake columns. The value of :math:`z_{lsa}` is chosen to satisfy the CFL condition for the model timestep. By default, :math:`\tilde{s}_{\min }` \ = 1 cm and :math:`s_{\min }` \ = 4 cm. See :ref:`Subin et al. (2012a; including Supporting Information)<Subinetal2012a>` for further discussion.

The rules for combining and sub-dividing snow layers (section :numref:`Snow Layer Combination and Subdivision`) are adjusted for lakes to maintain minimum thicknesses of :math:`s_{\min }` and to increase all target layer thicknesses by :math:`z_{lsa}`. The rules for combining layers are modified by simply increasing layer thickness thresholds by :math:`z_{lsa}`. The rules for dividing snow layers are contained in a separate subroutine that is modified for lakes, and is a function of the number of layers and the layer thicknesses. There are two types of operations: (a) subdividing layers in half, and (b) shifting some volume from higher layers to lower layers (without increasing the layer number). For subdivisions of type (a), the thickness thresholds triggering subdivision are increased by :math:`2z_{lsa}` for lakes. For shifts of type (b), the thickness thresholds triggering the shifts are increased by :math:`z_{lsa}`. At the end of the modified subroutine, a snow ice and liquid balance check are performed.

In rare instances, resolved snow layers may be present over an unfrozen top lake body layer. In this case, the snow layers may be eliminated if enough heat is present in the top layer to melt the snow: see :ref:`Subin et al. (2012a, Supporting Information) <Subinetal2012a>`.
