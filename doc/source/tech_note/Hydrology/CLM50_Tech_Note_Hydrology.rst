.. _rst_Hydrology:

Hydrology
============

The model parameterizes interception, throughfall, canopy drip, snow accumulation and melt, water transfer between snow layers, infiltration, evaporation, surface runoff, sub-surface drainage, redistribution within the soil column, and groundwater discharge and recharge to simulate changes in canopy water :math:`\Delta W_{can,\,liq}`, canopy snow water :math:`\Delta W_{can,\,sno}` surface water :math:`\Delta W_{sfc}`, snow water :math:`\Delta W_{sno}`, soil water :math:`\Delta w_{liq,\, i}`, and soil ice :math:`\Delta w_{ice,\, i}`, and water in the unconfined aquifer :math:`\Delta W_{a}` (all in kg m\ :sup:`-2` or mm of H\ :sub:`2`\ O) (:numref:`Figure Hydrologic processes`).

The total water balance of the system is

.. math::
   :label: 7.1

   \begin{array}{l} {\Delta W_{can,\,liq} +\Delta W_{can,\,sno} +\Delta W_{sfc} +\Delta W_{sno} +} \\ {\sum _{i=1}^{N_{levsoi} }\left(\Delta w_{liq,\, i} +\Delta w_{ice,\, i} \right)+\Delta W_{a} =\left(\begin{array}{l} {q_{rain} +q_{sno} -E_{v} -E_{g} -q_{over} } \\ {-q_{h2osfc} -q_{drai} -q_{rgwl} -q_{snwcp,\, ice} } \end{array}\right) \Delta t} \end{array}

where :math:`q_{rain}` is the liquid part of precipitation, :math:`q_{sno}` is the solid part of precipitation, :math:`E_{v}` is ET from vegetation (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`), :math:`E_{g}` is ground evaporation (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`), :math:`q_{over}` is surface runoff (section :numref:`Surface Runoff`), :math:`q_{h2osfc}` is runoff from surface water storage (section :numref:`Surface Runoff`), :math:`q_{drai}` is sub-surface drainage (section :numref:`Lateral Sub-surface Runoff`), :math:`q_{rgwl}` and :math:`q_{snwcp,ice}` are liquid and solid runoff from glaciers and lakes, and runoff from other surface types due to snow capping (section :numref:`Runoff from glaciers and snow-capped surfaces`) (all in kg m\ :sup:`-2` s\ :sup:`-1`), :math:`N_{levsoi}` is the number of soil layers (note that hydrology calculations are only done over soil layers 1 to :math:`N_{levsoi}`; ground levels :math:`N_{levsoi} +1` \ to :math:`N_{levgrnd}` are currently hydrologically inactive; :ref:`(Lawrence et al. 2008) <Lawrenceetal2008>` and :math:`\Delta t` is the time step (s).

.. _Figure Hydrologic processes:

.. Figure:: hydrologic.processes.png

 Hydrologic processes represented in CLM.

.. _Canopy Water:

Canopy Water
----------------

Liquid precipitation is either intercepted by the canopy, falls directly to the snow/soil surface (throughfall), or drips off the vegetation (canopy drip). Solid precipitation is treated similarly, with the addition of unloading of previously intercepted snow. Interception by vegetation is divided between liquid and solid phases :math:`q_{intr,\,liq}` and :math:`q_{intr,\,ice}` (kg m\ :sup:`-2` s\ :sup:`-1`)

.. math::
   :label: 7.2

   q_{intr,\,liq} = f_{pi,\,liq} \ q_{rain}

.. math::
   :label: 7.3

   q_{intr,\,ice} = f_{pi,\,ice} \ q_{sno}

where :math:`f_{pi,\,liq}` and :math:`f_{pi,\,ice}` are the fractions of intercepted precipitation of rain and snow, respectively

.. math::
   :label: 7.2b

   f_{pi,\,liq} = \alpha_{liq} \ tanh \left(L+S\right)

.. math::
   :label: 7.3b

   f_{pi,\,ice} =\alpha_{sno} \ \left\{1-\exp \left[-0.5\left(L+S\right)\right]\right\} \ ,

and :math:`L` and :math:`S` are the exposed leaf and stem area index, respectively (section :numref:`Phenology and vegetation burial by snow`), and the :math:`\alpha`\'s scale the fractional area of a leaf that collects water (:ref:`Lawrence et al. 2007 <Lawrenceetal2007>`). Default values of :math:`\alpha_{liq}` and :math:`\alpha_{sno}` are set to 1. Throughfall (kg m\ :sup:`-2` s\ :sup:`-1`) is also divided into liquid and solid phases, reaching the ground (soil or snow surface) as

.. math::
   :label: 7.4

   q_{thru,\, liq} = q_{rain} \left(1 - f_{pi,\,liq}\right)

.. math::
   :label: 7.5

   q_{thru,\, ice} = q_{sno} \left(1 - f_{pi,\,ice}\right)

Similarly, the liquid and solid canopy drip fluxes are

.. math::
   :label: 7.6

   q_{drip,\, liq} =\frac{W_{can,\,liq}^{intr} -W_{can,\,liq}^{max } }{\Delta t} \ge 0

.. math::
   :label: 7.7

   q_{drip,\, ice} =\frac{W_{can,\,sno}^{intr} -W_{can,\,sno}^{max } }{\Delta t} \ge 0

where

.. math::
   :label: 7.8

   W_{can,liq}^{intr} =W_{can,liq}^{n} +q_{intr,\, liq} \Delta t\ge 0

and

.. math::
   :label: 7.9

   W_{can,sno}^{intr} =W_{can,sno}^{n} +q_{intr,\, ice} \Delta t\ge 0

are the the canopy liquid water and snow water equivalent after accounting for interception, :math:`W_{can,\,liq}^{n}` and :math:`W_{can,\,sno}^{n}` are the canopy liquid and snow water from the previous time step, and :math:`W_{can,\,liq}^{max }` and :math:`W_{can,\,snow}^{max }` (kg m\ :sup:`-2` or mm of H\ :sub:`2`\ O) are the maximum amounts of liquid water and snow the canopy can hold. They are defined by

.. math::
   :label: 7.10

   W_{can,\,liq}^{max } =p_{liq}\left(L+S\right)

.. math::
   :label: 7.11

   W_{can,\,sno}^{max } =p_{sno}\left(L+S\right).

The maximum storage of liquid water is :math:`p_{liq}=0.1` kg m\ :sup:`-2` (:ref:`Dickinson et al. 1993 <Dickinsonetal1993>`), and that of snow is :math:`p_{sno}=6` kg m\ :sup:`-2`, consistent with reported field measurements (:ref:`Pomeroy et al. 1998 <Pomeroyetal1998>`).

Canopy snow unloading from wind speed :math:`u` and above-freezing temperatures are modeled from linear fluxes and e-folding times similar to :ref:`Roesch et al. (2001) <Roeschetal2001>`

.. math::
   :label: 7.12

   q_{unl,\, wind} =\frac{u W_{can,sno}}{1.56\times 10^5 \text{ m}}

.. math::
   :label: 7.13

   q_{unl,\, temp} =\frac{W_{can,sno}(T-270 \textrm{ K})}{1.87\times 10^5 \text{ K s}} > 0

.. math::
   :label: 7.14

   q_{unl,\, tot} =\min \left( q_{unl,\, wind} +q_{unl,\, temp} ,W_{can,\, sno} \right)

The canopy liquid water and snow water equivalent are updated as

.. math::
   :label: 7.15

    W_{can,\, liq}^{n+1} =W_{can,liq}^{n} + q_{intr,\, liq} - q_{drip,\, liq} \Delta t - E_{v}^{liq} \Delta t \ge 0

and

.. math::
   :label: 7.16

   W_{can,\, sno}^{n+1} =W_{can,sno}^{n} + q_{intr,\, ice} - \left(q_{drip,\, ice}+q_{unl,\, tot} \right)\Delta t
                         - E_{v}^{ice} \Delta t \ge 0

..   W_{can}^{n+1} =W_{can}^{n} +q_{intr} \Delta t-\left(q_{drip,\, liq} +q_{drip,\, ice} \right)\Delta t-E_{v}^{w} \Delta t\ge 0.

where :math:`E_{v}^{liq}` and :math:`E_{v}^{ice}` are partitioned from the stem and leaf surface evaporation :math:`E_{v}^{w}` (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) based on the vegetation temperature :math:`T_{v}` (K) (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) and its relation to the freezing temperature of water :math:`T_{f}` (K) (:numref:`Table Physical Constants`)

.. math::
   :label: 7.17

   E_{v}^{liq} =
   \left\{\begin{array}{lr}
   E_{v}^{w} &  T_v > T_{f} \\
   0         &  T_v \le T_f
   \end{array}\right\}

.. math::
   :label: 7.18

   E_{v}^{ice} =
   \left\{\begin{array}{lr}
   0         & T_v > T_f \\
   E_{v}^{w} & T_v \le T_f
   \end{array}\right\}.

..    \begin{array}{lr}
..    E_{v}^{liq} = E_{v}^{w} \qquad T > 273 \text{K}  \\
..    E_{v}^{ice} = E_{v}^{w} \qquad T \le 273 \text{K}
..    \end{array}

The total rate of liquid and solid precipitation reaching the ground is then

.. math::
   :label: 7.19

   q_{grnd,liq} =q_{thru,\, liq} +q_{drip,\, liq}

.. math::
   :label: 7.20

   q_{grnd,ice} =q_{thru,\, ice} +q_{drip,\, ice} +q_{unl,\, tot} .

Solid precipitation reaching the soil or snow surface, :math:`q_{grnd,\, ice} \Delta t`, is added immediately to the snow pack (Chapter :numref:`rst_Snow Hydrology`). The liquid part, :math:`q_{grnd,\, liq} \Delta t` is added after surface fluxes (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) and snow/soil temperatures (Chapter :numref:`rst_Soil and Snow Temperatures`) have been determined.

The wetted fraction of the canopy (stems plus leaves), which is required for surface flux (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) calculations, is (:ref:`Dickinson et al.1993 <Dickinsonetal1993>`)

.. math::
   :label: 7.21

   f_{wet} =
   \left\{\begin{array}{lr}
   \left[\frac{W_{can} }{p_{liq}\left(L+S\right)} \right]^{{2\mathord{\left/ {\vphantom {2 3}} \right.} 3} } \le 1 & \qquad L+S > 0 \\
   0 &\qquad L+S = 0
   \end{array}\right\}

while the fraction of the canopy that is dry and transpiring is

.. math::
   :label: 7.22

   f_{dry} =
   \left\{\begin{array}{lr}
   \frac{\left(1-f_{wet} \right)L}{L+S} & \qquad L+S > 0 \\
   0 &\qquad L+S = 0
   \end{array}\right\}.

Similarly, the snow-covered fraction of the canopy is used for surface alebdo when intercepted snow is present (Chapter :numref:`rst_Surface Albedos`)

.. math::
   :label: 7.23

   f_{can,\, sno} =
   \left\{\begin{array}{lr}
   \left[\frac{W_{can,\, sno} }{p_{sno}\left(L+S\right)} \right]^{{3\mathord{\left/ {\vphantom {3 20}} \right.} 20} } \le 1 & \qquad L+S > 0 \\
   0 &\qquad L+S = 0
   \end{array}\right\}.

.. _Surface Runoff, Surface Water Storage, and Infiltration:

Surface Runoff, Surface Water Storage, and Infiltration
-----------------------------------------------------------

The moisture input at the grid cell surface,\ :math:`q_{liq,\, 0}`, is the sum of liquid precipitation reaching the ground and melt water from snow (kg m\ :sup:`-2` s\ :sup:`-1`). The moisture flux is then partitioned between surface runoff, surface water storage, and infiltration into the soil.

.. _Surface Runoff:

Surface Runoff
^^^^^^^^^^^^^^^^^^^^

The simple TOPMODEL-based (:ref:`Beven and Kirkby 1979 <BevenKirkby1979>`) runoff model (SIMTOP) described by :ref:`Niu et al. (2005) <Niuetal2005>` is implemented to parameterize runoff. A key concept underlying this approach is that of fractional saturated area :math:`f_{sat}`, which is determined by the topographic characteristics and soil moisture state of a grid cell. The saturated portion of a grid cell contributes to surface runoff, :math:`q_{over}`, by the saturation excess mechanism (Dunne runoff)

.. math::
   :label: 7.64

   q_{over} =f_{sat} \ q_{liq,\, 0}

The fractional saturated area is a function of soil moisture

.. math::
   :label: 7.65

   f_{sat} =f_{\max } \ \exp \left(-0.5f_{over} z_{\nabla } \right)

where :math:`f_{\max }` is the potential or maximum value of :math:`f_{sat}`, :math:`f_{over}` is a decay factor (m\ :sup:`-1`), and :math:`z_{\nabla}` is the water table depth (m) (section :numref:`Lateral Sub-surface Runoff`). The maximum saturated fraction, :math:`f_{\max }`, is defined as the value of the discrete cumulative distribution function (CDF) of the topographic index when the grid cell mean water table depth is zero. Thus, :math:`f_{\max }` is the percent of pixels in a grid cell whose topographic index is larger than or equal to the grid cell mean topographic index. It should be calculated explicitly from the CDF at each grid cell at the resolution that the model is run. However, because this is a computationally intensive task for global applications, :math:`f_{\max }` is calculated once at 0.125° resolution using the 1-km compound topographic indices (CTIs) based on the HYDRO1K dataset (:ref:`Verdin and Greenlee 1996 <VerdinGreenlee1996>`) from USGS following the algorithm in :ref:`Niu et al. (2005) <Niuetal2005>` and then area-averaged to the desired model resolution (section :numref:`Surface Data`). Pixels with CTIs exceeding the 95 percentile threshold in each 0.125° grid cell are excluded from the calculation to eliminate biased estimation of statistics due to large CTI values at pixels on stream networks. For grid cells over regions without CTIs such as Australia, the global mean :math:`f_{\max }` is used to fill the gaps. See :ref:`Li et al. (2013b) <Lietal2013b>` for additional details. The decay factor :math:`f_{over}` for global simulations was determined through sensitivity analysis and comparison with observed runoff to be 0.5 m\ :sup:`-1`.

.. _Surface Water Storage:

Surface Water Storage
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A surface water store has been added to the model to represent wetlands and small, sub-grid scale water bodies. As a result, the wetland land unit has been removed as of CLM4.5. The state variables for surface water are the mass of water :math:`W_{sfc}` (kg m\ :sup:`-2`) and temperature :math:`T_{h2osfc}` (Chapter :numref:`rst_Soil and Snow Temperatures`). Surface water storage and outflow are functions of fine spatial scale elevation variations called microtopography. The microtopography is assumed to be distributed normally around the grid cell mean elevation. Given the standard deviation of the microtopographic distribution, :math:`\sigma _{micro}` (m), the fractional area of the grid cell that is inundated can be calculated. Surface water storage, :math:`Wsfc`, is related to the height (relative to the grid cell mean elevation) of the surface water, :math:`d`, by

.. math::
   :label: 7.66

   W_{sfc} =\frac{d}{2} \left(1+erf\left(\frac{d}{\sigma _{micro} \sqrt{2} } \right)\right)+\frac{\sigma _{micro} }{\sqrt{2\pi } } e^{\frac{-d^{2} }{2\sigma _{micro} ^{2} } }

where :math:`erf` is the error function. For a given value of :math:`W_{sfc}`, :eq:`7.66` can be solved for :math:`d` using the Newton-Raphson method. Once :math:`d` is known, one can determine the fraction of the area that is inundated as

.. math::
   :label: 7.67

   f_{h2osfc} =\frac{1}{2} \left(1+erf\left(\frac{d}{\sigma _{micro} \sqrt{2} } \right)\right)

No global datasets exist for microtopography, so the default parameterization is a simple function of slope

.. math::
   :label: 7.68

   \sigma _{micro} =\left(\beta +\beta _{0} \right)^{\eta }

where :math:`\beta` is the topographic slope, :math:`\beta_{0} =\left(\sigma_{\max } \right)^{\frac{1}{\eta } }` \ determines the maximum value of :math:`\sigma_{micro}`, and :math:`\eta` is an adjustable parameter. Default values in the model are :math:`\sigma_{\max } =0.4` and :math:`\eta =-3`.

If the spatial scale of the microtopography is small relative to that of the grid cell, one can assume that the inundated areas are distributed randomly within the grid cell. With this assumption, a result from percolation theory can be used to quantify the fraction of the inundated portion of the grid cell that is interconnected

.. math::
   :label: 7.69

   \begin{array}{lr} f_{connected} =\left(f_{h2osfc} -f_{c} \right)^{\mu } & \qquad f_{h2osfc} >f_{c}  \\ f_{connected} =0 &\qquad  f_{h2osfc} \le f_{c}  \end{array}

where :math:`f_{c}` is a threshold below which no single connected inundated area spans the grid cell and :math:`\mu` is a scaling exponent. Default values of :math:`f_{c}` and :math:`\mu` \ are 0.4 and 0.14, respectively. When the inundated fraction of the grid cell surpasses :math:`f_{c}`, the surface water store acts as a linear reservoir

.. math::
   :label: 7.70

   q_{out,h2osfc}=k_{h2osfc} \ f_{connected} \ (Wsfc-Wc)\frac{1}{\Delta t}

where :math:`q_{out,h2osfc}` is the surface water runoff, :math:`k_{h2osfc}` is a constant, :math:`Wc` is the amount of surface water present when :math:`f_{h2osfc} =f_{c}`, and :math:`\Delta t` is the model time step. The linear storage coefficent :math:`k_{h2osfc} = \sin \left(\beta \right)` is a function of grid cell mean topographic slope where :math:`\beta` is the slope in radians.

.. _Infiltration:

Infiltration
^^^^^^^^^^^^^^^^^^

The surface moisture flux remaining after surface runoff has been removed,

.. math::
   :label: 7.71

   q_{in,surface} = (1-f_{sat}) \ q_{liq,\, 0}

is divided into inputs to surface water (:math:`q_{in,\, h2osfc}` ) and the soil :math:`q_{in,soil}`. If :math:`q_{in,soil}` exceeds the maximum soil infiltration capacity (kg m\ :sup:`-2` s\ :sup:`-1`),

.. math::
   :label: 7.72

   q_{infl,\, \max } =(1-f_{sat}) \ \Theta_{ice} k_{sat}

where :math:`\Theta_{ice}` is an ice impedance factor (section :numref:`Hydraulic Properties`), infiltration excess (Hortonian) runoff is generated

.. math::
   :label: 7.73

   q_{infl,\, excess} =\max \left(q_{in,soil} -\left(1-f_{h2osfc} \right)q_{\inf l,\max } ,0\right)

and transferred from :math:`q_{in,soil}` to :math:`q_{in,h2osfc}`. After evaporative losses have been removed, these moisture fluxes are

.. math::
   :label: 7.74

   q_{in,\, h2osfc} = f_{h2osfc} q_{in,surface} + q_{infl,excess} - q_{evap,h2osfc}

and

.. math::
   :label: 7.75

   q_{in,soil} = (1-f_{h2osfc} ) \ q_{in,surface} - q_{\inf l,excess} - (1 - f_{sno} - f_{h2osfc} ) \ q_{evap,soil}.

The balance of surface water is then calculated as

.. math::
   :label: 7.76

   \Delta W_{sfc} =\left(q_{in,h2osfc} - q_{out,h2osfc} - q_{drain,h2osfc} \right) \ \Delta t.

Bottom drainage from the surface water store

.. math::
   :label: 7.77

   q_{drain,h2osfc} = \min \left(f_{h2osfc} q_{\inf l,\max } ,\frac{W_{sfc} }{\Delta t} \right)

is then added to :math:`q_{in,soil}` giving the total infiltration into the surface soil layer

.. math::
   :label: 7.78

   q_{infl} = q_{in,soil} + q_{drain,h2osfc}

Infiltration :math:`q_{infl}` and explicit surface runoff :math:`q_{over}` are not allowed for glaciers.

.. _Soil Water:

Soil Water
--------------

Soil water is predicted from a multi-layer model, in which the vertical soil moisture transport is governed by infiltration, surface and sub-surface runoff, gradient diffusion, gravity, and canopy transpiration through root extraction (:numref:`Figure Hydrologic processes`).

For one-dimensional vertical water flow in soils, the conservation of mass is stated as

.. math::
   :label: 7.79

   \frac{\partial \theta }{\partial t} =-\frac{\partial q}{\partial z} - e

where :math:`\theta` is the volumetric soil water content (mm\ :sup:`3` of water / mm\ :sup:`-3` of soil), :math:`t` is time (s), :math:`z` is height above some datum in the soil column (mm) (positive upwards), :math:`q` is soil water flux (kg m\ :sup:`-2` s\ :sup:`-1` or mm s\ :sup:`-1`) (positive upwards), and :math:`e` is a soil moisture sink term (mm of water mm\ :sup:`-1` of soil s\ :sup:`-1`) (ET loss). This equation is solved numerically by dividing the soil column into multiple layers in the vertical and integrating downward over each layer with an upper boundary condition of the infiltration flux into the top soil layer :math:`q_{infl}` and a zero-flux lower boundary condition at the bottom of the soil column (sub-surface runoff is removed later in the timestep, section :numref:`Lateral Sub-surface Runoff`).

The soil water flux :math:`q` in equation :eq:`7.79` can be described by Darcy's law :ref:`(Dingman 2002) <Dingman2002>`

.. math::
   :label: 7.80

   q = -k \frac{\partial \psi _{h} }{\partial z}

where :math:`k` is the hydraulic conductivity (mm s\ :sup:`-1`), and :math:`\psi _{h}` is the hydraulic potential (mm). The hydraulic potential is

.. math::
   :label: 7.81

   \psi _{h} =\psi _{m} +\psi _{z}

where :math:`\psi _{m}` is the soil matric potential (mm) (which is related to the adsorptive and capillary forces within the soil matrix), and :math:`\psi _{z}` is the gravitational potential (mm) (the vertical distance from an arbitrary reference elevation to a point in the soil). If the reference elevation is the soil surface, then :math:`\psi _{z} =z`. Letting :math:`\psi =\psi _{m}`, Darcy's law becomes

.. math::
   :label: 7.82

   q = -k \left[\frac{\partial \left(\psi +z\right)}{\partial z} \right].

Equation :eq:`7.82` can be further manipulated to yield

.. math::
   :label: 7.83

   q = -k \left[\frac{\partial \left(\psi +z\right)}{\partial z} \right]
   = -k \left(\frac{\partial \psi }{\partial z} + 1 \right) \ .

Substitution of this equation into equation :eq:`7.79`, with :math:`e = 0`, yields the Richards equation :ref:`(Dingman 2002) <Dingman2002>`

.. math::
   :label: 7.84

   \frac{\partial \theta }{\partial t} =
   \frac{\partial }{\partial z} \left[k\left(\frac{\partial \psi }{\partial z} + 1
   \right)\right].

In practice (Section :numref:`Numerical Solution Hydrology`), changes in soil water content are predicted from :eq:`7.79` using finite-difference approximations for :eq:`7.84`.

.. _Hydraulic Properties:

Hydraulic Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

The hydraulic conductivity :math:`k_{i}` (mm s\ :sup:`-1`) and the soil matric potential :math:`\psi _{i}` (mm) for layer :math:`i` vary with volumetric soil water :math:`\theta _{i}` and soil texture. As with the soil thermal properties (section :numref:`Soil And Snow Thermal Properties`) the hydraulic properties of the soil are assumed to be a weighted combination of the mineral properties, which are determined according to sand and clay contents based on work by :ref:`Clapp and Hornberger (1978) <ClappHornberger1978>` and :ref:`Cosby et al. (1984) <Cosbyetal1984>`, and organic properties of the soil (:ref:`Lawrence and Slater 2008 <LawrenceSlater2008>`).

The hydraulic conductivity is defined at the depth of the interface of two adjacent layers :math:`z_{h,\, i}` (:numref:`Figure Water flux schematic`) and is a function of the saturated hydraulic conductivity :math:`k_{sat} \left[z_{h,\, i} \right]`, the liquid volumetric soil moisture of the two layers :math:`\theta _{i}` and :math:`\theta_{i+1}` and an ice impedance factor :math:`\Theta_{ice}`

.. math::
   :label: 7.85

   k\left[z_{h,\, i} \right] =
   \left\{\begin{array}{lr}
   \Theta_{ice} k_{sat} \left[z_{h,\, i} \right]\left[\frac{0.5\left(\theta_{\, i} +\theta_{\, i+1} \right)}{0.5\left(\theta_{sat,\, i} +\theta_{sat,\, i+1} \right)} \right]^{2B_{i} +3} & \qquad 1 \le i \le N_{levsoi} - 1 \\
   \Theta_{ice} k_{sat} \left[z_{h,\, i} \right]\left(\frac{\theta_{\, i} }{\theta_{sat,\, i} } \right)^{2B_{i} +3} & \qquad i = N_{levsoi}
   \end{array}\right\}.

The ice impedance factor is a function of ice content, and is meant to quantify the increased tortuosity of the water flow when part of the pore space is filled with ice. :ref:`Swenson et al. (2012) <Swensonetal2012>` used a power law form

.. math::
   :label: 7.86

   \Theta_{ice} = 10^{-\Omega F_{ice} }

where :math:`\Omega = 6`\ and :math:`F_{ice} = \frac{\theta_{ice} }{\theta_{sat} }` is the ice-filled fraction of the pore space.

Because the hydraulic properties of mineral and organic soil may differ significantly, the bulk hydraulic properties of each soil layer are computed as weighted averages of the properties of the mineral and organic components. The water content at saturation (i.e. porosity) is

.. math::
   :label: 7.90

   \theta_{sat,i} =(1-f_{om,i} )\theta_{sat,\min ,i} +f_{om,i} \theta_{sat,om}

where :math:`f_{om,i}` is the soil organic matter fraction, :math:`\theta_{sat,om}` is the porosity of organic matter, and the porosity of the mineral soil :math:`\theta_{sat,\min,i}` is

.. math::
   :label: 7.91

   \theta_{sat,\min ,i} = 0.489 - 0.00126(\% sand)_{i} .

The exponent :math:`B_{i}` is

.. math::
   :label: 7.92

   B_{i} =(1-f_{om,i} )B_{\min ,i} +f_{om,i} B_{om}

where :math:`B_{om}` is for organic matter and

.. math::
   :label: 7.93

   B_{\min ,i} =2.91+0.159(\% clay)_{i} .

The soil matric potential (mm) is defined at the node depth :math:`z_{i}` of each layer :math:`i` (:numref:`Figure Water flux schematic`)

.. math::
   :label: 7.94

   \psi _{i} =\psi _{sat,\, i} \left(\frac{\theta_{\, i} }{\theta_{sat,\, i} } \right)^{-B_{i} } \ge -1\times 10^{8} \qquad 0.01\le \frac{\theta_{i} }{\theta_{sat,\, i} } \le 1

where the saturated soil matric potential (mm) is

.. math::
   :label: 7.95

   \psi _{sat,i} =(1-f_{om,i} )\psi _{sat,\min ,i} +f_{om,i} \psi _{sat,om}

where :math:`\psi _{sat,om}` \ is the saturated organic matter matric potential and the saturated mineral soil matric potential :math:`\psi _{sat,\min,i}` \ is

.. math::
   :label: 7.96

   \psi _{sat,\, \min ,\, i} =-10.0\times 10^{1.88-0.0131(\% sand)_{i} } .

The saturated hydraulic conductivity, :math:`k_{sat} \left[z_{h,\, i} \right]` (mm s\ :sup:`-1`), for organic soils (:math:`k_{sat,\, om}` ) may be two to three orders of magnitude larger than that of mineral soils (:math:`k_{sat,\, \min }` ). Bulk soil layer values of :math:`k_{sat}` \ calculated as weighted averages based on :math:`f_{om}` may therefore be determined primarily by the organic soil properties even for values of :math:`f_{om}` as low as 1 %. To better represent the influence of organic soil material on the grid cell average saturated hydraulic conductivity, the soil organic matter fraction is further subdivided into "connected" and "unconnected" fractions using a result from percolation theory (:ref:`Stauffer and Aharony 1994 <StaufferAharony1994>`, :ref:`Berkowitz and Balberg 1992 <BerkowitzBalberg1992>`). Assuming that the organic and mineral fractions are randomly distributed throughout a soil layer, percolation theory predicts that above a threshold value :math:`f_{om} = f_{threshold}`, connected flow pathways consisting of organic material only exist and span the soil space. Flow through these pathways interacts only with organic material, and thus can be described by :math:`k_{sat,\, om}`. This fraction of the grid cell is given by

.. math::
   :label: 7.97

   \begin{array}{lr}
   f_{perc} =\; N_{perc} \left(f_{om} {\rm \; }-f_{threshold} \right)^{\beta_{perc} } f_{om} {\rm \; } & \qquad f_{om} \ge f_{threshold}  \\
   f_{perc} = 0 & \qquad f_{om} <f_{threshold}
   \end{array}

where :math:`\beta ^{perc} =0.139`, :math:`f_{threshold} =0.5`, and :math:`N_{perc} =\left(1-f_{threshold} \right)^{-\beta_{perc} }`. In the unconnected portion of the grid cell, :math:`f_{uncon} =\; \left(1-f_{perc} {\rm \; }\right)`, the saturated hydraulic conductivity is assumed to correspond to flow pathways that pass through the mineral and organic components in series

.. math::
   :label: 7.98

   k_{sat,\, uncon} =f_{uncon} \left(\frac{\left(1-f_{om} \right)}{k_{sat,\, \min } } +\frac{\left(f_{om} -f_{perc} \right)}{k_{sat,\, om} } \right)^{-1} .

where saturated hydraulic conductivity for mineral soil depends on soil texture (:ref:`Cosby et al. 1984 <Cosbyetal1984>`) as

.. math::
   :label: 7.99

   k_{sat,\, \min } \left[z_{h,\, i} \right]=0.0070556\times 10^{-0.884+0.0153\left(\% sand\right)_{i} } .

The bulk soil layer saturated hydraulic conductivity is then computed as

.. math::
   :label: 7.100

   k_{sat} \left[z_{h,\, i} \right]=f_{uncon,\, i} k_{sat,\, uncon} \left[z_{h,\, i} \right]+(1-f_{uncon,\, i} )k_{sat,\, om} \left[z_{h,\, i} \right].

The soil organic matter properties implicitly account for the standard observed profile of organic matter properties as

.. math::
   :label: 1.101

   \theta_{sat,om} = max(0.93 - 0.1\times z_{i} / zsapric, 0.83).

.. math::
   :label: 1.102

   B_{om} = min(2.7 + 9.3\times z_{i} / zsapric, 12.0).

.. math::
   :label: 1.103

   \psi_{sat,om} = min(10.3 - 0.2\times z_{i} / zsapric, 10.1).

.. math::
   :label: 1.104

   k_{sat,om} = max(0.28 - 0.2799\times z_{i} / zsapric, k_{sat,\, \min } \left[z_{h,\, i} \right]).

where :math:`zsapric =0.5` \m is the depth that organic matter takes on the characteristics of sapric peat.

.. _Numerical Solution Hydrology:

Numerical Solution
^^^^^^^^^^^^^^^^^^^^^^^^

With reference to :numref:`Figure Water flux schematic`, the equation for conservation of mass (equation :eq:`7.79`) can be integrated over each layer as

.. math::
   :label: 7.101

   \int _{-z_{h,\, i} }^{-z_{h,\, i-1} }\frac{\partial \theta }{\partial t} \,  dz=-\int _{-z_{h,\, i} }^{-z_{h,\, i-1} }\frac{\partial q}{\partial z}  \, dz-\int _{-z_{h,\, i} }^{-z_{h,\, i-1} } e\, dz .

Note that the integration limits are negative since :math:`z` is defined as positive upward from the soil surface. This equation can be written as

.. math::
   :label: 7.102

   \Delta z_{i} \frac{\partial \theta_{liq,\, i} }{\partial t} =-q_{i-1} +q_{i} -e_{i}

where :math:`q_{i}` is the flux of water across interface :math:`z_{h,\, i}`, :math:`q_{i-1}` is the flux of water across interface :math:`z_{h,\, i-1}`, and :math:`e_{i}` is a layer-averaged soil moisture sink term (ET loss) defined as positive for flow out of the layer (mm s\ :sup:`-1`). Taking the finite difference with time and evaluating the fluxes implicitly at time :math:`n+1` yields

.. math::
   :label: 7.103

   \frac{\Delta z_{i} \Delta \theta_{liq,\, i} }{\Delta t} =-q_{i-1}^{n+1} +q_{i}^{n+1} -e_{i}

where :math:`\Delta \theta_{liq,\, i} =\theta_{liq,\, i}^{n+1} -\theta_{liq,\, i}^{n}` is the change in volumetric soil liquid water of layer :math:`i` in time :math:`\Delta t`\ and :math:`\Delta z_{i}` is the thickness of layer :math:`i` (mm).

The water removed by transpiration in each layer :math:`e_{i}` is a function of the total transpiration :math:`E_{v}^{t}` (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`) and the effective root fraction :math:`r_{e,\, i}`

.. math::
   :label: 7.104

   e_{i} =r_{e,\, i} E_{v}^{t} .

.. _Figure Water flux schematic:

.. Figure:: image2.png

 Schematic diagram of numerical scheme used to solve for soil water fluxes.

Shown are three soil layers, :math:`i-1`, :math:`i`, and :math:`i+1`. The soil matric potential :math:`\psi` and volumetric soil water :math:`\theta_{liq}` are defined at the layer node depth :math:`z`. The hydraulic conductivity :math:`k\left[z_{h} \right]` is defined at the interface of two layers :math:`z_{h}`. The layer thickness is :math:`\Delta z`. The soil water fluxes :math:`q_{i-1}` and :math:`q_{i}` are defined as positive upwards. The soil moisture sink term :math:`e` (ET loss) is defined as positive for flow out of the layer.

Note that because more than one plant functional type (PFT) may share a soil column, the transpiration :math:`E_{v}^{t}` is a weighted sum of transpiration from all PFTs whose weighting depends on PFT area as

.. math::
   :label: 7.105

   E_{v}^{t} =\sum _{j=1}^{npft}\left(E_{v}^{t} \right)_{j} \left(wt\right)_{j}

where :math:`npft` is the number of PFTs sharing a soil column, :math:`\left(E_{v}^{t} \right)_{j}` is the transpiration from the :math:`j^{th}` PFT on the column, and :math:`\left(wt\right)_{j}` is the relative area of the :math:`j^{th}` PFT with respect to the column. The effective root fraction :math:`r_{e,\, i}` is also a column-level quantity that is a weighted sum over all PFTs. The weighting depends on the per unit area transpiration of each PFT and its relative area as

.. math::
   :label: 7.106

   r_{e,\, i} =\frac{\sum _{j=1}^{npft}\left(r_{e,\, i} \right)_{j} \left(E_{v}^{t} \right)_{j} \left(wt\right)_{j}  }{\sum _{j=1}^{npft}\left(E_{v}^{t} \right)_{j} \left(wt\right)_{j}  }

where :math:`\left(r_{e,\, i} \right)_{j}` is the effective root fraction for the :math:`j^{th}` PFT

.. math::
   :label: 7.107

   \begin{array}{lr}
   \left(r_{e,\, i} \right)_{j} =\frac{\left(r_{i} \right)_{j} \left(w_{i} \right)_{j} }{\left(\beta _{t} \right)_{j} } & \qquad \left(\beta _{t} \right)_{j} >0 \\
   \left(r_{e,\, i} \right)_{j} =0 & \qquad \left(\beta _{t} \right)_{j} =0
   \end{array}

and :math:`\left(r_{i} \right)_{j}` is the fraction of roots in layer :math:`i` (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`), :math:`\left(w_{i} \right)_{j}` is a soil dryness or plant wilting factor for layer :math:`i` (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`), and :math:`\left(\beta_{t} \right)_{j}` is a wetness factor for the total soil column for the :math:`j^{th}` PFT (Chapter :numref:`rst_Stomatal Resistance and Photosynthesis`).

The soil water fluxes in :eq:`7.103`,, which are a function of :math:`\theta_{liq,\, i}` and :math:`\theta_{liq,\, i+1}` because of their dependence on hydraulic conductivity and soil matric potential, can be linearized about :math:`\theta` using a Taylor series expansion as

.. math::
   :label: 7.108

   q_{i}^{n+1} =q_{i}^{n} +\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } \Delta \theta_{liq,\, i} +\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} } \Delta \theta_{liq,\, i+1}

.. math::
   :label: 7.109

   q_{i-1}^{n+1} =q_{i-1}^{n} +\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} } \Delta \theta_{liq,\, i-1} +\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } \Delta \theta_{liq,\, i} .

Substitution of these expressions for :math:`q_{i}^{n+1}` and :math:`q_{i-1}^{n+1}` into :eq:`7.103` results in a general tridiagonal equation set of the form

.. math::
   :label: 7.110

   r_{i} =a_{i} \Delta \theta_{liq,\, i-1} +b_{i} \Delta \theta_{liq,\, i} +c_{i} \Delta \theta_{liq,\, i+1}

where

.. math::
   :label: 7.111

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} }

.. math::
   :label: 7.112

   b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.113

   c_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} }

.. math::
   :label: 7.114

   r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .

The tridiagonal equation set is solved over :math:`i=1,\ldots,N_{levsoi}`.

The finite-difference forms of the fluxes and partial derivatives in equations :eq:`7.111` - :eq:`7.114` can be obtained from equation :eq:`7.82` as

.. math::
   :label: 7.115

   q_{i-1}^{n} =-k\left[z_{h,\, i-1} \right]\left[\frac{\left(\psi _{i-1} -\psi _{i} \right)+\left(z_{i} - z_{i-1} \right)}{z_{i} -z_{i-1} } \right]

.. math::
   :label: 7.116

   q_{i}^{n} =-k\left[z_{h,\, i} \right]\left[\frac{\left(\psi _{i} -\psi _{i+1} \right)+\left(z_{i+1} - z_{i} \right)}{z_{i+1} -z_{i} } \right]

.. math::
   :label: 7.117

   \frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} } =-\left[\frac{k\left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \frac{\partial \psi _{i-1} }{\partial \theta _{liq,\, i-1} } \right]-\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i-1} } \left[\frac{\left(\psi _{i-1} -\psi _{i} \right)+\left(z_{i} - z_{i-1} \right)}{z_{i} - z_{i-1} } \right]

.. math::
   :label: 7.118

   \frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } =\left[\frac{k\left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \frac{\partial \psi _{i} }{\partial \theta _{liq,\, i} } \right]-\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i} } \left[\frac{\left(\psi _{i-1} -\psi _{i} \right)+\left(z_{i} - z_{i-1} \right)}{z_{i} - z_{i-1} } \right]

.. math::
   :label: 7.119

   \frac{\partial q_{i} }{\partial \theta _{liq,\, i} } =-\left[\frac{k\left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \frac{\partial \psi _{i} }{\partial \theta _{liq,\, i} } \right]-\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i} } \left[\frac{\left(\psi _{i} -\psi _{i+1} \right)+\left(z_{i+1} - z_{i} \right)}{z_{i+1} - z_{i} } \right]

.. math::
   :label: 7.120

   \frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} } =\left[\frac{k\left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \frac{\partial \psi _{i+1} }{\partial \theta _{liq,\, i+1} } \right]-\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i+1} } \left[\frac{\left(\psi _{i} -\psi _{i+1} \right)+\left(z_{i+1} - z_{i} \right)}{z_{i+1} - z_{i} } \right].

The derivatives of the soil matric potential at the node depth are derived from :eq:`7.94`

.. math::
   :label: 7.121

   \frac{\partial \psi _{i-1} }{\partial \theta_{liq,\, \, i-1} } =-B_{i-1} \frac{\psi _{i-1} }{\theta_{\, \, i-1} }

.. math::
   :label: 7.122

   \frac{\partial \psi _{i} }{\partial \theta_{\, liq,\, i} } =-B_{i} \frac{\psi _{i} }{\theta_{i} }

.. math::
   :label: 7.123

   \frac{\partial \psi _{i+1} }{\partial \theta_{liq,\, i+1} } =-B_{i+1} \frac{\psi _{i+1} }{\theta_{\, i+1} }

with the constraint :math:`0.01\, \theta_{sat,\, i} \le \theta_{\, i} \le \theta_{sat,\, i}`.

The derivatives of the hydraulic conductivity at the layer interface are derived from :eq:`7.85`

.. math::
   :label: 7.124

   \begin{array}{l}
   {\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i-1} }
   = \frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i} }
   = \left(2B_{i-1} +3\right) \ \overline{\Theta}_{ice} \ k_{sat} \left[z_{h,\, i-1} \right] \ \left[\frac{\overline{\theta}_{liq}}{\overline{\theta}_{sat}} \right]^{2B_{i-1} +2} \left(\frac{0.5}{\overline{\theta}_{sat}} \right)} \end{array}

where :math:`\overline{\Theta}_{ice} = \Theta(\overline{\theta}_{ice})` :eq:`7.86`,
:math:`\overline{\theta}_{ice} = 0.5\left(\theta_{ice\, i-1} +\theta_{ice\, i} \right)`,
:math:`\overline{\theta}_{liq} = 0.5\left(\theta_{liq\, i-1} +\theta_{liq\, i} \right)`,
and
:math:`\overline{\theta}_{sat} = 0.5\left(\theta_{sat,\, i-1} +\theta_{sat,\, i} \right)`

and

.. math::
   :label: 7.125

   \begin{array}{l}
   {\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i} }
   = \frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i+1} }
   = \left(2B_{i} +3\right) \ \overline{\Theta}_{ice} \ k_{sat} \left[z_{h,\, i} \right] \ \left[\frac{\overline{\theta}_{liq}}{\overline{\theta}_{sat}} \right]^{2B_{i} +2} \left(\frac{0.5}{\overline{\theta}_{sat}} \right)} \end{array}.

where :math:`\overline{\theta}_{liq} = 0.5\left(\theta_{\, i} +\theta_{\, i+1} \right)`,
:math:`\overline{\theta}_{sat} = 0.5\left(\theta_{sat,\, i} +\theta_{sat,\, i+1} \right)`.

Equation set for layer :math:`i=1`
''''''''''''''''''''''''''''''''''''''''''

For the top soil layer (:math:`i=1`), the boundary condition is the infiltration rate (section :numref:`Surface Runoff`), :math:`q_{i-1}^{n+1} =-q_{infl}^{n+1}`, and the water balance equation is

.. math::
   :label: 7.135

   \frac{\Delta z_{i} \Delta \theta_{liq,\, i} }{\Delta t} =q_{infl}^{n+1} +q_{i}^{n+1} -e_{i} .

After grouping like terms, the coefficients of the tridiagonal set of equations for :math:`i=1` are

.. math::
   :label: 7.136

   a_{i} =0

.. math::
   :label: 7.137

   b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.138

   c_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} }

.. math::
   :label: 7.139

   r_{i} =q_{infl}^{n+1} -q_{i}^{n} +e_{i} .

Equation set for layers :math:`i=2,\ldots ,N_{levsoi} -1`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The coefficients of the tridiagonal set of equations for :math:`i=2,\ldots,N_{levsoi} -1` are

.. math::
   :label: 7.140

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} }

.. math::
   :label: 7.141

   b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.142

   c_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i+1} }

.. math::
   :label: 7.143

   r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .

Equation set for layer :math:`i=N_{levsoi}`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

For the lowest soil layer (:math:`i=N_{levsoi}` ), a zero-flux bottom boundary condition is applied (:math:`q_{i}^{n} =0`) and the coefficients of the tridiagonal set of equations for :math:`i=N_{levsoi}` are

.. math::
   :label: 7.148

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i-1} }

.. math::
   :label: 7.149

   b_{i} =\frac{\partial q_{i} }{\partial \theta_{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta_{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.150

   c_{i} =0

.. math::
   :label: 7.151

   r_{i} =q_{i-1}^{n} +e_{i} .

Adaptive Time Stepping
'''''''''''''''''''''''''''''

The length of the time step is adjusted in order to improve the accuracy and stability of the numerical solutions. The difference between two numerical approximations is used to estimate the temporal truncation error, and then the step size :math:`\Delta t_{sub}` is adjusted to meet a user-prescribed error tolerance :ref:`[Kavetski et al., 2002]<Kavetskietal2002>`. The temporal truncation error is estimated by comparing the flux obtained from the first-order Taylor series expansion (:math:`q_{i-1}^{n+1}` and :math:`q_{i}^{n+1}`, equations :eq:`7.108` and :eq:`7.109`) against the flux at the start of the time step (:math:`q_{i-1}^{n}` and :math:`q_{i}^{n}`). Since the tridiagonal solution already provides an estimate of :math:`\Delta \theta_{liq,i}`, it is convenient to compute the error for each of the :math:`i` layers from equation :eq:`7.103` as

.. math::
   :label: 7.152

   \epsilon_{i} = \left[ \frac{\Delta \theta_{liq,\, i} \Delta z_{i}}{\Delta t_{sub}} -
   \left( q_{i-1}^{n} - q_{i}^{n} + e_{i}\right) \right] \ \frac{\Delta t_{sub}}{2}

and the maximum absolute error across all layers as

.. math::
   :label: 7.153

   \begin{array}{lr}
   \epsilon_{crit} = {\rm max} \left( \left| \epsilon_{i} \right| \right) & \qquad 1 \le i \le nlevsoi
   \end{array} \ .

The adaptive step size selection is based on specified upper and lower error tolerances, :math:`\tau_{U}` and :math:`\tau_{L}`. The solution is accepted if :math:`\epsilon_{crit} \le \tau_{U}` and the procedure repeats until the adaptive sub-stepping spans the full model time step (the sub-steps are doubled if :math:`\epsilon_{crit} \le \tau_{L}`, i.e., if the solution is very accurate). Conversely, the solution is rejected if :math:`\epsilon_{crit} > \tau_{U}`. In this case the length of the sub-steps is halved and a new solution is obtained. The halving of substeps continues until either :math:`\epsilon_{crit} \le \tau_{U}` or the specified minimum time step length is reached.

Upon solution of the tridiagonal equation set, the liquid water contents are updated as follows

.. math::
   :label: 7.164

   w_{liq,\, i}^{n+1} =w_{liq,\, i}^{n} +\Delta \theta_{liq,\, i} \Delta z_{i} \qquad i=1,\ldots ,N_{levsoi} .

The volumetric water content is

.. math::
   :label: 7.165

   \theta_{i} =\frac{w_{liq,\, i} }{\Delta z_{i} \rho _{liq} } +\frac{w_{ice,\, i} }{\Delta z_{i} \rho _{ice} } .

.. _Frozen Soils and Perched Water Table:

Frozen Soils and Perched Water Table
----------------------------------------

When soils freeze, the power-law form of the ice impedance factor (section :numref:`Hydraulic Properties`) can greatly decrease the hydraulic conductivity of the soil, leading to nearly impermeable soil layers. When unfrozen soil layers are present above relatively ice-rich frozen layers, the possibility exists for perched saturated zones. Lateral drainage from perched saturated regions is parameterized as a function of the thickness of the saturated zone

.. math::
   :label: 7.166

   q_{drai,perch} =k_{drai,\, perch} \left(z_{frost} -z_{\nabla ,perch} \right)

where :math:`k_{drai,\, perch}` depends on topographic slope and soil hydraulic conductivity,

.. math::
   :label: 7.167

   k_{drai,\, perch} =10^{-5} \sin (\beta )\left(\frac{\sum _{i=N_{perch} }^{i=N_{frost} }\Theta_{ice,i} k_{sat} \left[z_{i} \right]\Delta z_{i}  }{\sum _{i=N_{perch} }^{i=N_{frost} }\Delta z_{i}  } \right)

where :math:`\Theta_{ice}` is an ice impedance factor,
:math:`\beta` is the mean grid cell topographic slope in radians,
:math:`z_{frost}` \ is the depth to the frost table, and
:math:`z_{\nabla,perch}` is the depth to the perched saturated zone. The frost table :math:`z_{frost}` is defined as the shallowest frozen layer having an unfrozen layer above it, while the perched water table :math:`z_{\nabla,perch}` is defined as the depth at which the volumetric water content drops below a specified threshold. The default threshold is set to 0.9. Drainage from the perched saturated zone :math:`q_{drai,perch}` is removed from layers :math:`N_{perch}` through :math:`N_{frost}`, which are the layers containing :math:`z_{\nabla,perch}` and, :math:`z_{frost}` \ respectively.

.. _Lateral Sub-surface Runoff:

Lateral Sub-surface Runoff
---------------------------------------
Lateral sub-surface runoff occurs when saturated soil moisture conditions exist within the soil column. Sub-surface runoff is

.. math::
   :label: 7.168

   q_{drai} = \Theta_{ice} K_{baseflow} tan \left( \beta \right)
   \Delta z_{sat}^{N_{baseflow}} \ ,

where :math:`K_{baseflow}` is a calibration parameter, :math:`\beta` is the topographic slope, the exponent :math:`N_{baseflow}` = 1, and :math:`\Delta z_{sat}` is the thickness of the saturated portion of the soil column.

The saturated thickness is

.. math::
   :label: 7.1681

   \Delta z_{sat} = z_{bedrock} - z_{\nabla},

where the water table :math:`z_{\nabla}` is determined by finding the first soil layer above the bedrock depth (section :numref:`Depth to Bedrock`) in which the volumetric water content drops below a specified threshold. The default threshold is set to 0.9.

The specific yield, :math:`S_{y}`, which depends on the soil properties and the water table location, is derived by taking the difference between two equilibrium soil moisture profiles whose water tables differ by an infinitesimal amount

.. math::
   :label: 7.174

   S_{y} =\theta_{sat} \left(1-\left(1+\frac{z_{\nabla } }{\Psi _{sat} } \right)^{\frac{-1}{B} } \right)

where B is the Clapp-Hornberger exponent. Because :math:`S_{y}` is a function of the soil properties, it results in water table dynamics that are consistent with the soil water fluxes described in section :numref:`Soil Water`.

After the above calculations, two numerical adjustments are implemented to keep the liquid water content of each soil layer (:math:`w_{liq,\, i}` ) within physical constraints of :math:`w_{liq}^{\min } \le w_{liq,\, i} \le \left(\theta_{sat,\, i} -\theta_{ice,\, i} \right)\Delta z_{i}` where :math:`w_{liq}^{\min } =0.01` (mm). First, beginning with the bottom soil layer :math:`i=N_{levsoi}`, any excess liquid water in each soil layer (:math:`w_{liq,\, i}^{excess} =w_{liq,\, i} -\left(\theta_{sat,\, i} -\theta_{ice,\, i} \right)\Delta z_{i} \ge 0`) is successively added to the layer above. Any excess liquid water that remains after saturating the entire soil column is added to drainage :math:`q_{drai}`. Second, to prevent negative :math:`w_{liq,\, i}`, each layer is successively brought up to :math:`w_{liq,\, i} =w_{liq}^{\min }` by taking the required amount of water from the layer below. If this results in :math:`w_{liq,\, N_{levsoi} } <w_{liq}^{\min }`, then the layers above are searched in succession for the required amount of water (:math:`w_{liq}^{\min } -w_{liq,\, N_{levsoi} }` ) and removed from those layers subject to the constraint :math:`w_{liq,\, i} \ge w_{liq}^{\min }`. If sufficient water is not found, then the water is removed from :math:`W_{t}` and :math:`q_{drai}`.

The soil surface layer liquid water and ice contents are then updated for dew :math:`q_{sdew}`, frost :math:`q_{frost}`, or sublimation :math:`q_{subl}` (section :numref:`Update of Ground Sensible and Latent Heat Fluxes`) as

.. math::
   :label: 7.175

   w_{liq,\, 1}^{n+1} =w_{liq,\, 1}^{n} +q_{sdew} \Delta t

.. math::
   :label: 7.176

   w_{ice,\, 1}^{n+1} =w_{ice,\, 1}^{n} +q_{frost} \Delta t

.. math::
   :label: 7.177

   w_{ice,\, 1}^{n+1} =w_{ice,\, 1}^{n} -q_{subl} \Delta t.

Sublimation of ice is limited to the amount of ice available.

.. _Runoff from glaciers and snow-capped surfaces:

Runoff from glaciers and snow-capped surfaces
-------------------------------------------------

All surfaces are constrained to have a snow water equivalent :math:`W_{sno} \le W_{cap} = 10,000` kg m\ :sup:`-2`. For snow-capped columns, any addition of mass at the top (precipitation, dew/riping) is balanced by an equally large mass flux at the bottom of the snow column. This so-called capping flux is separated into solid :math:`q_{snwcp,ice}` \ and liquid :math:`q_{snwcp,liq}` runoff terms. The partitioning of these phases is based on the phase ratio in the bottom snow layer at the time of the capping, such that phase ratio in this layer is unaltered.

The :math:`q_{snwcp,ice}` runoff is sent to the River Transport Model (RTM) (Chapter 11) where it is routed to the ocean as an ice stream and, if applicable, the ice is melted there.

For snow-capped surfaces other than glaciers and lakes the :math:`q_{snwcp,liq}` runoff is assigned to the glaciers and lakes runoff term :math:`q_{rgwl}` (e.g. :math:`q_{rgwl} =q_{snwcp,liq}` ). For glacier surfaces the runoff term :math:`q_{rgwl}` is calculated from the residual of the water balance

.. math::
   :label: 7.180

   q_{rgwl} =q_{grnd,ice} +q_{grnd,liq} -E_{g} -E_{v} -\frac{\left(W_{b}^{n+1} -W_{b}^{n} \right)}{\Delta t} -q_{snwcp,ice}

where :math:`W_{b}^{n}` and :math:`W_{b}^{n+1}` are the water balances at the beginning and ending of the time step defined as

.. math::
   :label: 7.181

   W_{b} =W_{can} +W_{sno} +\sum _{i=1}^{N}\left(w_{ice,i} +w_{liq,i} \right) .

Currently, glaciers are non-vegetated and :math:`E_{v} =W_{can} =0`. The contribution of lake runoff to :math:`q_{rgwl}` is described in section :numref:`Precipitation, Evaporation, and Runoff Lake`. The runoff term :math:`q_{rgwl}` may be negative for glaciers and lakes, which reduces the total amount of runoff available to the river routing model (Chapter :numref:`rst_River Transport Model (RTM)`).
