.. _rst_Hydrology:

Hydrology
============

The model parameterizes interception, throughfall, canopy drip, snow
accumulation and melt, water transfer between snow layers, infiltration,
evaporation, surface runoff, sub-surface drainage, redistribution within
the soil column, and groundwater discharge and recharge to simulate
changes in canopy water :math:`\Delta W_{can}` , surface water
:math:`\Delta W_{sfc}` , snow water :math:`\Delta W_{sno}` , soil water
:math:`\Delta w_{liq,\, i}` , and soil ice :math:`\Delta w_{ice,\, i}` ,
and water in the unconfined aquifer :math:`\Delta W_{a}`  (all in kg
m\ :sup:`-2` or mm of H\ :sub:`2`\ O) (Figure 7.1).

The total water balance of the system is

.. math::
   :label: 7.1) 

   \begin{array}{l} {\Delta W_{can} +\Delta W_{sfc} +\Delta W_{sno} +} \\ {\sum _{i=1}^{N_{levsoi} }\left(\Delta w_{liq,\, i} +\Delta w_{ice,\, i} \right)+\Delta W_{a} =\left(\begin{array}{l} {q_{rain} +q_{sno} -E_{v} -E_{g} -q_{over} } \\ {-q_{h2osfc} -q_{drai} -q_{rgwl} -q_{snwcp,\, ice} } \end{array}\right) \Delta t} \end{array}

where :math:`q_{rain}`  is the liquid part of precipitation,
:math:`q_{sno}`  is the solid part of precipitation, :math:`E_{v}`  is
ET from vegetation (Chapter 5), :math:`E_{g}`  is ground evaporation
(Chapter 5), :math:`q_{over}`  is surface runoff (section 7.3),
:math:`q_{h2osfc}`  is runoff from surface water storage (section 7.3),
:math:`q_{drai}`  is sub-surface drainage (section 7.6),
:math:`q_{rgwl}`  and :math:`q_{snwcp,ice}`  are liquid and solid runoff
from glaciers, wetlands, and lakes, and runoff from other surface types
due to snow capping (section 7.7) (all in kg m\ :sup:`-2`
s\ :sup:`-1`), :math:`N_{levsoi}`  is the number of soil layers
(note that hydrology calculations are only done over soil layers 1 to
:math:`N_{levsoi}` ; ground levels :math:`N_{levsoi} +1`\ to
:math:`N_{levgrnd}`  are currently hydrologically inactive; Lawrence et
al. 2008) and :math:`\Delta t` is the time step (s).

Figure 7.1. Hydrologic processes represented in CLM.

.. image:: image1.png

.. _Canopy Water:

Canopy Water
----------------

Precipitation is either intercepted by the canopy, falls directly to the
snow/soil surface (throughfall), or drips off the vegetation (canopy
drip). Interception by vegetation :math:`q_{intr}`  (kg
m\ :sup:`-2` s\ :sup:`-1`) does not distinguish between
liquid and solid phases

.. math::
   :label: 7.2) 

   q_{intr} =\alpha \left(q_{rain} +q_{sno} \right)\left\{1-\exp \left[-0.5\left(L+S\right)\right]\right\}

where :math:`L` and :math:`S` are the exposed leaf and stem area index,
respectively (section 2.1.4), and :math:`\alpha =0.25` scales
interception from point to grid cell (Lawrence et al. 2007). Throughfall
(kg m\ :sup:`-2` s\ :sup:`-1`), however, is divided into
liquid and solid phases reaching the ground (soil or snow surface) as

.. math::
   :label: 7.3) 

   q_{thru,\, liq} =q_{rain} \left[1-\alpha \left\{1-\exp \left[-0.5\left(L+S\right)\right]\right\}\right]

.. math::
   :label: 7.4) 

   q_{thru,\, ice} =q_{sno} \left[1-\alpha \left\{1-\exp \left[-0.5\left(L+S\right)\right]\right\}\right].

Similarly, the canopy drip is

.. math::
   :label: 7.5) 

   q_{drip,\, liq} =\frac{W_{can}^{intr} -W_{can,\, \max } }{\Delta t} \frac{q_{rain} }{q_{rain} +q_{sno} } \ge 0

.. math::
   :label: 7.6) 

   q_{drip,\, ice} =\frac{W_{can}^{intr} -W_{can,\, \max } }{\Delta t} \frac{q_{sno} }{q_{rain} +q_{sno} } \ge 0

where

.. math::
   :label: 7.7) 

   W_{can}^{intr} =W_{can}^{n} +q_{intr} \Delta t\ge 0

is the canopy water after accounting for interception,
:math:`W_{can}^{n}`  is the canopy water from the previous time step,
and :math:`W_{can,\, \max }`  (kg m\ :sup:`-2`) is the maximum
amount of water the canopy can hold

.. math::
   :label: 7.8) 

   W_{can,\, \max } =p\left(L+S\right).

The maximum storage of solid water is assumed to be the same as that of
liquid water, :math:`p=0.1` kg m\ :sup:`-2` (Dickinson et al.
1993). The canopy water is updated as

.. math::
   :label: 7.9) 

   W_{can}^{n+1} =W_{can}^{n} +q_{intr} \Delta t-\left(q_{drip,\, liq} +q_{drip,\, ice} \right)\Delta t-E_{v}^{w} \Delta t\ge 0.

where :math:`E_{v}^{w}`  is the flux of water vapor from stem and leaf
surfaces (Chapter 5). The total rate of liquid and solid precipitation
reaching the ground is then

.. math::
   :label: ZEqnNum946822 

   q_{grnd,liq} =q_{thru,\, liq} +q_{drip,\, liq}

.. math::
   :label: ZEqnNum339590 

   q_{grnd,ice} =q_{thru,\, ice} +q_{drip,\, ice} .

Solid precipitation reaching the soil or snow surface,
:math:`q_{grnd,\, ice} \Delta t`, is added immediately to the snow pack
(section 7.2). The liquid part, :math:`q_{grnd,\, liq} \Delta t` is
added after surface fluxes (Chapter 5) and snow/soil temperatures
(Chapter 6) have been determined.

The wetted fraction of the canopy (stems plus leaves), which is required
for the surface albedo (section 3.1) and surface flux (Chapter 5)
calculations is (Dickinson et al. 1993)

.. math::
   :label: 7.12) 

   f_{wet} =\left\{\begin{array}{l} {\left[\frac{W_{can} }{p\left(L+S\right)} \right]^{{2\mathord{\left/ {\vphantom {2 3}} \right. \kern-\nulldelimiterspace} 3} } \le 1\qquad L+S>0} \\ {0\qquad L+S=0} \end{array}\right\}

while the fraction of the canopy that is dry and transpiring is

.. math::
   :label: 7.13) 

   f_{dry} =\left\{\begin{array}{l} {\frac{\left(1-f_{wet} \right)L}{L+S} \qquad L+S>0} \\ {0\qquad L+S=0} \end{array}\right\}.

.. _Surface Runoff, Surface Water Storage, and Infiltration:

Surface Runoff, Surface Water Storage, and Infiltration
-----------------------------------------------------------

The moisture input at the grid cell surface ,\ :math:`q_{liq,\, 0}` , is
the sum of liquid precipitation reaching the ground and melt water from
snow (kg m\ :sup:`-2` s\ :sup:`-1`). The moisture flux is
then partitioned between surface runoff, surface water storage, and
infiltration into the soil.

.. _Surface Runoff:

Surface Runoff
^^^^^^^^^^^^^^^^^^^^

The simple TOPMODEL-based (Beven and Kirkby 1979) runoff model (SIMTOP)
described by Niu et al. (2005) is implemented to parameterize runoff. A
key concept underlying this approach is that of fractional saturated
area :math:`f_{sat}` , which is determined by the topographic
characteristics and soil moisture state of a grid cell. The saturated
portion of a grid cell contributes to surface runoff, :math:`q_{over}` ,
by the saturation excess mechanism (Dunne runoff)

.. math::
   :label: ZEqnNum549608 

   q_{over} =f_{sat} q_{liq,\, 0}

The fractional saturated area is a function of soil moisture

.. math::
   :label: 7.65) 

   f_{sat} =f_{\max } \exp \left(-0.5f_{over} z_{\nabla } \right)

where :math:`f_{\max }`  is the potential or maximum value of
:math:`f_{sat}` , :math:`f_{over}`  is a decay factor
(m\ :sup:`-1`), and :math:`z_{\nabla}` is the water table depth
(m) (section 7.6). The maximum saturated fraction, :math:`f_{\max }` ,
is defined as the value of the discrete cumulative distribution function
(CDF) of the topographic index when the grid cell mean water table depth
is zero. Thus, :math:`f_{\max }`  is the percent of pixels in a grid
cell whose topographic index is larger than or equal to the grid cell
mean topographic index. It should be calculated explicitly from the CDF
at each grid cell at the resolution that the model is run. However,
because this is a computationally intensive task for global
applications, :math:`f_{\max }`  is calculated once at
0.125\ :sup:`o` resolution using the 1-km compound topographic
indices (CTIs) based on the HYDRO1K dataset (Verdin and Greenlee 1996)
from USGS following the algorithm in Niu et al. (2005) and then
area-averaged to the desired model resolution (section 2.2.3). Pixels
with CTIs exceeding the 95 percentile threshold in each
0.125\ :sup:`o` grid cell are excluded from the calculation to
eliminate biased estimation of statistics due to large CTI values at
pixels on stream networks. For grid cells over regions without CTIs such
as Australia, the global mean :math:`f_{\max }`  is used to fill the
gaps. See Li et al. (2013b) for additional details. The decay factor
:math:`f_{over}`  for global simulations was determined through
sensitivity analysis and comparison with observed runoff to be 0.5
m\ :sup:`-1`.

.. _Surface Water Storage:

Surface Water Storage
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A surface water store has been added to the model to represent wetlands
and small, sub-grid scale water bodies. As a result, the wetland land
unit has been removed. The state variables for surface water are the
mass of water :math:`W_{sfc}`  (kg m\ :sup:`-2`) and temperature
:math:`T_{h2osfc}`  (Chapter 6). Surface water storage and outflow are
functions of fine spatial scale elevation variations called
microtopography. The microtopography is assumed to be distributed
normally around the grid cell mean elevation. Given the standard
deviation of the microtopographic distribution, :math:`\sigma _{micro}` 
(m), the fractional area of the grid cell that is inundated can be
calculated. Surface water storage, :math:`Wsfc`, is related to the
height (relative to the grid cell mean elevation) of the surface water,
:math:`d`, by

.. math::
   :label: ZEqnNum277892 

   W_{sfc} =\frac{d}{2} \left(1+erf\left(\frac{d}{\sigma _{micro} \sqrt{2} } \right)\right)+\frac{\sigma _{micro} }{\sqrt{2\pi } } e^{\frac{-d^{2} }{2\sigma _{micro} ^{2} } }

where :math:`erf` is the error function. For a given value of
:math:`W_{sfc}` , equation can be solved for :math:`d` using the
Newton-Raphson method. Once :math:`d` is known, one can determine the
fraction of the area that is inundated as

.. math::
   :label: 7.67) 

   f_{h2osfc} =\frac{1}{2} \left(1+erf\left(\frac{d}{\sigma _{micro} \sqrt{2} } \right)\right)

No global datasets exist for microtopography, so the default
parameterization is a simple function of slope

.. math::
   :label: 7.68) 

   \sigma _{micro} =\left(\beta +\beta _{0} \right)^{\eta }

where :math:`\beta`  is the topographic slope,
:math:`\beta _{0} =\left(\sigma _{\max } \right)^{\frac{1}{\eta } }` \ determines
the maximum value of :math:`\sigma _{}` , and :math:`\eta`  is an
adjustable parameter. Default values in the model are
:math:`\sigma _{\max } =0.4` and :math:`\eta =-3`.

If the spatial scale of the microtopography is small relative to that of
the grid cell, one can assume that the inundated areas are distributed
randomly within the grid cell. With this assumption, a result from
percolation theory can be used to quantify the fraction of the inundated
portion of the grid cell that is interconnected

.. math::
   :label: 7.69) 

   \begin{array}{l} {f_{connected} =\left(f_{h2osfc} -f_{c} \right)^{\mu } \qquad for\, f_{h2osfc} >f_{c} } \\ {f_{connected} =0\qquad \qquad \qquad for\, f_{h2osfc} \le f_{c} } \end{array}

where :math:`f_{c}`  is a threshold below which no single connected
inundated area spans the grid cell and :math:`\mu`  is a scaling
exponent. Default values of :math:`f_{c}`  and :math:`\mu` \ are 0.4 and
0.14, respectively. When the inundated fraction of the grid cell
surpasses :math:`f_{c}` , the surface water store acts as a linear
reservoir

.. math::
   :label: 7.70) 

   qout,h2osfc=kh2osfcf_{connected} (Wsfc-Wc)\frac{1}{\Delta t}

where :math:`qout,h2osfc` is the surface water runoff, :math:`kh2osfc`
is a constant, :math:`Wc` is the amount of surface water present when
:math:`f_{h2osfc} =f_{c}` , and :math:`\Delta t` is the model time step.
The linear storage coefficent :math:`kh2osfc=\sin \left(\beta \right)`
is a function of grid cell mean topographic slope where :math:`\beta` 
is the slope in radians.

.. _Infiltration:

Infiltration
^^^^^^^^^^^^^^^^^^

The surface moisture flux remaining after surface runoff has been
removed,

.. math::
   :label: 7.71) 

   qin,surface=(1-f_{sat} )q_{liq,\, 0}

is divided into inputs to surface water (:math:`q_{in,\, h2osfc}` ) and
the soil :math:`q_{in,soil}` . If :math:`q_{in,soil}`  exceeds the
maximum soil infiltration capacity (kg m\ :sup:`-2`
s\ :sup:`-1`),

.. math::
   :label: ZEqnNum569150 

   q_{infl,\, \max } =(1-fsat){\rm \Theta }iceksat

where :math:`{\rm \Theta }ice` is an ice impedance factor (section
7.4), infiltration excess (Hortonian) runoff is generated

.. math::
   :label: 7.73) 

   q_{infl,\, excess} =\max \left(q_{in,soil} -\left(1-f_{h2osfc} \right)q_{\inf l,\max } ,0\right)

and transferred from :math:`q_{in,soil}`  to :math:`q_{in,h2osfc}` .
After evaporative losses have been removed, these moisture fluxes are

.. math::
   :label: 7.74) 

   q_{in,\, h2osfc} =f_{h2osfc} q_{in,surface} +q_{infl,excess} -q_{evap,h2osfc}

and

.. math::
   :label: 7.75) 

   qin,soil=(1-f_{h2osfc} )q_{in,surface} -q_{\inf l,excess} -(1-f_{sno} -f_{h2osfc} )qevap,soil.

The balance of surface water is then calculated as

.. math::
   :label: 7.76) 

   \Delta W_{sfc} =\left(qin,h2osfc-qout,h2osfc-qdrain,h2osfc\right)\Delta t.

Bottom drainage from the surface water store

.. math::
   :label: 7.77) 

   qdrain,h2osfc=\min \left(f_{h2osfc} q_{\inf l,\max } ,\frac{W_{sfc} }{\Delta t} \right)

is then added to :math:`q_{in,soil}`  giving the total infiltration
into the surface soil layer

.. math::
   :label: 7.78) 

   q_{infl} =q_{in,soil} +q_{drain,h2osfc}

Infiltration :math:`q_{infl}`  and explicit surface runoff
:math:`q_{over}`  are not allowed for glaciers.

.. _Soil Water:

Soil Water
--------------

Soil water is predicted from a multi-layer model, in which the vertical
soil moisture transport is governed by infiltration, surface and
sub-surface runoff, gradient diffusion, gravity, canopy transpiration
through root extraction, and interactions with groundwater (Figure 7.1).
The following derivation generally follows that of Z.-L. Yang (1998,
unpublished manuscript) with modifications by Zeng and Decker (2009).

For one-dimensional vertical water flow in soils, the conservation of
mass is stated as

.. math::
   :label: ZEqnNum790844 

   \frac{\partial \theta }{\partial t} =-\frac{\partial q}{\partial z} -Q

where :math:`\theta`  is the volumetric soil water content
(mm:sup:`3` of water mm\ :sup:`-3` of soil), :math:`t` is
time (s), :math:`z` is height above some datum in the soil column (mm)
(positive upwards), :math:`q` is soil water flux (kg m\ :sup:`-2`
s\ :sup:`-1` or mm s\ :sup:`-1`) (positive upwards), and
:math:`Q` is a soil moisture sink term (mm of water mm\ :sup:`-1`
of soil s\ :sup:`-1`) (ET loss). This equation is solved
numerically by dividing the soil column into multiple layers in the
vertical and integrating downward over each layer with an upper boundary
condition of the infiltration flux into the top soil layer
:math:`q_{infl}`  and a lower boundary condition that depends on the
depth of the water table.

The soil water flux :math:`q` in equation can be described by Darcy’s
law

.. math::
   :label: 7.80) 

   q=-k\frac{\partial \psi _{h} }{\partial z}

where :math:`k` is the hydraulic conductivity (mm s\ :sup:`-1`),
and :math:`\psi _{h}`  is the hydraulic potential (mm). The hydraulic
potential is

.. math::
   :label: 7.81) 

   \psi _{h} =\psi _{m} +\psi _{z}

where :math:`\psi _{m}`  is the soil matric potential (mm) (which is
related to the adsorptive and capillary forces within the soil matrix),
and :math:`\psi _{z}`  is the gravitational potential (mm) (the vertical
distance from an arbitrary reference elevation to a point in the soil).
If the reference elevation is the soil surface, then
:math:`\psi _{z} =z`. Letting :math:`\psi =\psi _{m}` , Darcy’s law
becomes

.. math::
   :label: ZEqnNum186573 

   q=-k\left[\frac{\partial \left(\psi +z\right)}{\partial z} \right].

Darcy’s equation can be further manipulated to yield

.. math::
   :label: 7.83) 

   q=-k\left[\frac{\partial \left(\psi +z\right)}{\partial z} \right]=-k\left(\frac{\partial \psi }{\partial z} +1\right)=-k\left(\frac{\partial \theta }{\partial z} \frac{\partial \psi }{\partial \theta } +1\right).

Substitution of this equation into equation , with :math:`Q=0`, yields
the Richards equation

.. math::
   :label: ZEqnNum670361 

   \frac{\partial \theta }{\partial t} =\frac{\partial }{\partial z} \left[k\left(\frac{\partial \theta }{\partial z} \frac{\partial \psi }{\partial \theta } \right)+1\right].

Zeng and Decker (2009) note that this :math:`\theta` -based form of the
Richards equation cannot maintain the hydrostatic equilibrium soil
moisture distribution because of the truncation errors of the
finite-difference numerical scheme. They show that this deficiency can
be overcome by subtracting the equilibrium state from equation as

.. math::
   :label: ZEqnNum936839 

   q=-k\left[\frac{\partial \left(\psi +z-C\right)}{\partial z} \right]

where :math:`C` is a constant hydraulic potential above the water table
:math:`z_{\nabla }` 

.. math::
   :label: ZEqnNum126975 

   C=\psi _{E} +z=\psi _{sat} \left[\frac{\theta _{E} \left(z\right)}{\theta _{sat} } \right]^{-B} +z=\psi _{sat} +z_{\nabla }

so that

.. math::
   :label: ZEqnNum537733 

   q=-k\left[\frac{\partial \left(\psi -\psi _{E} \right)}{\partial z} \right]

where :math:`\psi _{E}` \ is the equilibrium soil matric potential
(mm). Substitution of equations and into equation yields Zeng and
Decker’s (2009) modified Richards equation

.. math::
   :label: 7.88) 

   \frac{\partial \theta }{\partial t} =\frac{\partial }{\partial z} \left[k\left(\frac{\partial \left(\psi -\psi _{E} \right)}{\partial z} \right)\right]-Q

where the soil moisture sink term :math:`Q` is now included.

.. _Hydraulic Properties:

Hydraulic Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^

The hydraulic conductivity :math:`k_{i}`  (mm s\ :sup:`-1`) and
the soil matric potential :math:`\psi _{i}`  (mm) for layer :math:`i`
vary with volumetric soil water :math:`\theta _{i}`  and soil texture.
As with the soil thermal properties (section 6.3) the hydraulic
properties of the soil are assumed to be a weighted combination of the
mineral properties, which are determined according to sand and clay
contents based on work by Clapp and Hornberger (1978) and Cosby et al.
(1984), and organic properties of the soil (Lawrence and Slater 2008).

The hydraulic conductivity is defined at the depth of the interface of
two adjacent layers :math:`z_{h,\, i}`  (Figure 7.3) and is a function
of the saturated hydraulic conductivity
:math:`k_{sat} \left[z_{h,\, i} \right]`, the liquid volumetric soil
moisture of the two layers :math:`\theta _{i}`  and
:math:`\theta _{i+1}`  and an ice impedance factor
:math:`\Theta _{ice}` 

.. math::
   :label: ZEqnNum398074 

   k\left[z_{h,\, i} \right]=\left\{\begin{array}{l} {\Theta _{ice} k_{sat} \left[z_{h,\, i} \right]\left[\frac{0.5\left(\theta _{\, i} +\theta _{\, i+1} \right)}{0.5\left(\theta _{sat,\, i} +\theta _{sat,\, i+1} \right)} \right]^{2B_{i} +3} \qquad 1\le i\le N_{levsoi} -1} \\ {\Theta _{ice} k_{sat} \left[z_{h,\, i} \right]\left(\frac{\theta _{\, i} }{\theta _{sat,\, i} } \right)^{2B_{i} +3} \qquad i=N_{levsoi} } \end{array}\right\}.

The ice impedance factor is a function of ice content, and is meant to
quantify the increased tortuosity of the water flow when part of the
pore space is filled with ice. Swenson et al. (2012) used a power law
form :math:`\Theta _{ice} =10^{-\Omega F_{ice} }`  where
:math:`\Omega =6`\ and
:math:`F_{ice} =\frac{\theta _{ice} }{\theta _{sat} }`  is the
ice-filled fraction of the pore space.

Because the hydraulic properties of mineral and organic soil may differ
significantly, the bulk hydraulic properties of each soil layer are
computed as weighted averages of the properties of the mineral and
organic components. The water content at saturation (i.e. porosity) is

.. math::
   :label: 7.90) 

   \theta _{sat,i} =(1-f_{om,i} )\theta _{sat,\min ,i} +f_{om,i} \theta _{sat,om}

where :math:`f_{om,i}`  is the soil organic matter fraction,
:math:`\theta _{sat,om} =0.9` (Farouki 1981; Letts et al. 2000) is the
porosity of organic matter and the porosity of the mineral soil
:math:`\theta _{sat,\min ,i}`  is

.. math::
   :label: 7.91) 

   \theta _{sai,\min ,i} =0.489-0.00126(\% sand)_{i} .

The exponent “:math:`B`” is

.. math::
   :label: 7.92) 

   B_{i} =(1-f_{om,i} )B_{\min ,i} +f_{om,i} B_{om}

where :math:`B_{om} =2.7`\ (Letts et al. 2000) and

.. math::
   :label: 7.93) 

   B_{\min ,i} =2.91+0.159(\% clay)_{i} .

The soil matric potential (mm) is defined at the node depth
:math:`z_{i}`  of each layer :math:`i` (Figure 7.3)

.. math::
   :label: ZEqnNum316201 

   \psi _{i} =\psi _{sat,\, i} \left(\frac{\theta _{\, i} }{\theta _{sat,\, i} } \right)^{-B_{i} } \ge -1\times 10^{8} \qquad 0.01\le \frac{\theta _{i} }{\theta _{sat,\, i} } \le 1

where the saturated soil matric potential (mm) is

.. math::
   :label: 7.95) 

   \psi _{sat,i} =(1-f_{om,i} )\psi _{sat,\min ,i} +f_{om,i} \psi _{sat,om}

where :math:`\psi _{sat,om} =-10.3`\ mm (Letts et al. 2000) is the
saturated organic matter matric potential and the saturated mineral soil
matric potential :math:`\psi _{sat,\min ,i}` \ is

.. math::
   :label: 7.96) 

   \psi _{sat,\, \min ,\, i} =-10.0\times 10^{1.88-0.0131(\% sand)_{i} } .

The saturated hydraulic conductivity,
:math:`k_{sat} \left[z_{h,\, i} \right]` (mm s\ :sup:`-1`), for
organic soils (:math:`k_{sat,\, om}` ) may be two to three orders of
magnitude larger than that of mineral soils (:math:`k_{sat,\, \min }` ).
Bulk soil layer values of :math:`k_{sat}` \ calculated as weighted
averages based on :math:`f_{om}`  may therefore be determined primarily
by the organic soil properties even for values of :math:`f_{om}`  as low
as 1 %. To better represent the influence of organic soil material on
the grid cell average saturated hydraulic conductivity, the soil organic
matter fraction is further subdivided into “connected” and “unconnected”
fractions using a result from percolation theory (Stauffer and Aharony
1994, Berkowitz and Balberg 1992). Assuming that the organic and mineral
fractions are randomly distributed throughout a soil layer, percolation
theory predicts that above a threshold value
:math:`f_{om} =f_{threshold}` , connected flow pathways consisting of
organic material only exist and span the soil space. Flow through these
pathways interacts only with organic material, and thus can be described
by :math:`k_{sat,\, om}` . This fraction of the grid cell is given by

.. math::
   :label: 7.97) 

   \begin{array}{l} {f_{perc} =\; N_{perc} \left(f_{om} {\rm \; }-f_{threshold} \right)^{\beta _{perc} } f_{om} {\rm \; }\qquad f_{om} \ge f_{threshold} } \\ {f_{perc} =0\qquad f_{om} <f_{threshold} } \end{array}

where :math:`\beta ^{perc} =0.139`, :math:`f_{threshold} =0.5`, and
:math:`N_{perc} =\left(1-f_{threshold} \right)^{-\beta _{perc} }` . In
the unconnected portion of the grid cell,
:math:`f_{uncon} =\; \left(1-f_{perc} {\rm \; }\right)`, the saturated
hydraulic conductivity is assumed to correspond to flow pathways that
pass through the mineral and organic components in series

.. math::
   :label: 7.98) 

   k_{sat,\, uncon} =f_{uncon} \left(\frac{\left(1-f_{om} \right)}{k_{sat,\, \min } } +\frac{\left(f_{om} -f_{perc} \right)}{k_{sat,\, om} } \right)^{-1} .

where saturated hydraulic conductivity for mineral soil depends on soil
texture (Cosby et al. 1984) as

.. math::
   :label: 7.99) 

   k_{sat,\, \min } \left[z_{h,\, i} \right]=0.0070556\times 10^{-0.884+0.0153\left(\% sand\right)_{i} } .

The bulk soil layer saturated hydraulic conductivity is then computed
as

.. math::
   :label: 7.100) 

   k_{sat} \left[z_{h,\, i} \right]=f_{uncon,\, i} k_{sat,\, uncon} \left[z_{h,\, i} \right]+(1-f_{uncon,\, i} )k_{sat,\, om} \left[z_{h,\, i} \right].

.. _Numerical Solution Hydrology:

Numerical Solution
^^^^^^^^^^^^^^^^^^^^^^^^

With reference to Figure 7.3, the equation for conservation of mass
(equation ) can be integrated over each layer as

.. math::
   :label: 7.101) 

   \int _{-z_{h,\, i} }^{-z_{h,\, i-1} }\frac{\partial \theta }{\partial t} \,  dz=-\int _{-z_{h,\, i} }^{-z_{h,\, i-1} }\frac{\partial q}{\partial z}  \, dz-\int _{-z_{h,\, i} }^{-z_{h,\, i-1} }Q\, dz .

Note that the integration limits are negative since :math:`z` is defined
as positive upward from the soil surface. This equation can be written
as

.. math::
   :label: ZEqnNum225309 

   \Delta z_{i} \frac{\partial \theta _{liq,\, i} }{\partial t} =-q_{i-1} +q_{i} -e_{i}

where :math:`q_{i}`  is the flux of water across interface
:math:`z_{h,\, i}` , :math:`q_{i-1}`  is the flux of water across
interface :math:`z_{h,\, i-1}` , and :math:`e_{i}`  is a layer-averaged
soil moisture sink term (ET loss) defined as positive for flow out of
the layer (mm s\ :sup:`-1`). Taking the finite difference with
time and evaluating the fluxes implicitly at time :math:`n+1` yields

.. math::
   :label: ZEqnNum181361 

   \frac{\Delta z_{i} \Delta \theta _{liq,\, i} }{\Delta t} =-q_{i-1}^{n+1} +q_{i}^{n+1} -e_{i}

where
:math:`\Delta \theta _{liq,\, i} =\theta _{liq,\, i}^{n+1} -\theta _{liq,\, i}^{n}` 
is the change in volumetric soil liquid water of layer :math:`i` in time
:math:`\Delta t`\ and :math:`\Delta z_{i}`  is the thickness of layer
:math:`i` (mm).

The water removed by transpiration in each layer :math:`e_{i}`  is a
function of the total transpiration :math:`E_{v}^{t}`  (Chapter 5) and
the effective root fraction :math:`r_{e,\, i}` 

.. math::
   :label: ZEqnNum357392 

   e_{i} =r_{e,\, i} E_{v}^{t} .

Figure 7.3. Schematic diagram of numerical scheme used to solve for soil
water fluxes.

Shown are three soil layers, :math:`i-1`, :math:`i`, and :math:`i+1`.
The soil matric potential :math:`\psi`  and volumetric soil water
:math:`\theta _{liq}`  are defined at the layer node depth :math:`z`.
The hydraulic conductivity :math:`k\left[z_{h} \right]` is defined at
the interface of two layers :math:`z_{h}` . The layer thickness is
:math:`\Delta z`. The soil water fluxes :math:`q_{i-1}`  and
:math:`q_{i}`  are defined as positive upwards. The soil moisture sink
term :math:`e` (ET loss) is defined as positive for flow out of the
layer.

.. image:: image2.png

Note that because more than one plant functional type (PFT) may share a
soil column, the transpiration :math:`E_{v}^{t}`  is a weighted sum of
transpiration from all PFTs whose weighting depends on PFT area as

.. math::
   :label: 7.105) 

   E_{v}^{t} =\sum _{j=1}^{npft}\left(E_{v}^{t} \right)_{j} \left(wt\right)_{j}

where :math:`npft` is the number of PFTs sharing a soil column,
:math:`\left(E_{v}^{t} \right)_{j}`  is the transpiration from the
:math:`j^{th}`  PFT on the column, and :math:`\left(wt\right)_{j}`  is
the relative area of the :math:`j^{th}`  PFT with respect to the column.
The effective root fraction :math:`r_{e,\, i}`  is also a column-level
quantity that is a weighted sum over all PFTs. The weighting depends on
the per unit area transpiration of each PFT and its relative area as

.. math::
   :label: 7.106) 

   r_{e,\, i} =\frac{\sum _{j=1}^{npft}\left(r_{e,\, i} \right)_{j} \left(E_{v}^{t} \right)_{j} \left(wt\right)_{j}  }{\sum _{j=1}^{npft}\left(E_{v}^{t} \right)_{j} \left(wt\right)_{j}  }

where :math:`\left(r_{e,\, i} \right)_{j}`  is the effective root
fraction for the :math:`j^{th}`  PFT

.. math::
   :label: 7.107) 

   \begin{array}{l} {\left(r_{e,\, i} \right)_{j} =\frac{\left(r_{i} \right)_{j} \left(w_{i} \right)_{j} }{\left(\beta _{t} \right)_{j} } \qquad \left(\beta _{t} \right)_{j} >0} \\ {\left(r_{e,\, i} \right)_{j} =0\qquad \left(\beta _{t} \right)_{j} =0} \end{array}

and :math:`\left(r_{i} \right)_{j}`  is the fraction of roots in layer
:math:`i` (Chapter 8), :math:`\left(w_{i} \right)_{j}`  is a soil
dryness or plant wilting factor for layer :math:`i` (Chapter 8), and
:math:`\left(\beta _{t} \right)_{j}`  is a wetness factor for the total
soil column for the :math:`j^{th}`  PFT (Chapter 8).

The soil water fluxes in equation , which are a function of
:math:`\theta _{liq,\, i}`  and :math:`\theta _{liq,\, i+1}`  because of
their dependence on hydraulic conductivity and soil matric potential,
can be linearized about :math:`\theta`  using a Taylor series expansion
as

.. math::
   :label: 7.108) 

   q_{i}^{n+1} =q_{i}^{n} +\frac{\partial q_{i} }{\partial \theta _{liq,\, i} } \Delta \theta _{liq,\, i} +\frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} } \Delta \theta _{liq,\, i+1}

.. math::
   :label: 7.109) 

   q_{i-1}^{n+1} =q_{i-1}^{n} +\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} } \Delta \theta _{liq,\, i-1} +\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } \Delta \theta _{liq,\, i} .

Substitution of these expressions for :math:`q_{i}^{n+1}`  and
:math:`q_{i-1}^{n+1}`  into equation results in a general tridiagonal
equation set of the form

.. math::
   :label: 7.110) 

   r_{i} =a_{i} \Delta \theta _{liq,\, i-1} +b_{i} \Delta \theta _{liq,\, i} +c_{i} \Delta \theta _{liq,\, i+1}

where

.. math::
   :label: ZEqnNum557934 

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} }

.. math::
   :label: 7.112) 

   b_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.113) 

   c_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} }

.. math::
   :label: ZEqnNum981892 

   r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .

The tridiagonal equation set is solved over
:math:`i=1,\ldots ,N_{levsoi} +1` where the layer
:math:`i=N_{levsoi} +1` is a virtual layer representing the aquifer.

The finite-difference forms of the fluxes and partial derivatives in
equations - can be obtained from equation as

.. math::
   :label: 7.115) 

   q_{i-1}^{n} =-k\left[z_{h,\, i-1} \right]\left[\frac{\left(\psi _{i-1} -\psi _{i} \right)+\left(\psi _{E,\, i} -\psi _{E,\, i-1} \right)}{z_{i} -z_{i-1} } \right]

.. math::
   :label: 7.116) 

   q_{i}^{n} =-k\left[z_{h,\, i} \right]\left[\frac{\left(\psi _{i} -\psi _{i+1} \right)+\left(\psi _{E,\, i+1} -\psi _{E,\, i} \right)}{z_{i+1} -z_{i} } \right]

.. math::
   :label: 7.117) 

   \frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} } =-\left[\frac{k\left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \frac{\partial \psi _{i-1} }{\partial \theta _{liq,\, i-1} } \right]-\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i-1} } \left[\frac{\left(\psi _{i-1} -\psi _{i} \right)+\left(\psi _{E,\, i} -\psi _{E,\, i-1} \right)}{z_{i} -z_{i-1} } \right]

.. math::
   :label: 7.118) 

   \frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } =\left[\frac{k\left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \frac{\partial \psi _{i} }{\partial \theta _{liq,\, i} } \right]-\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i} } \left[\frac{\left(\psi _{i-1} -\psi _{i} \right)+\left(\psi _{E,\, i} -\psi _{E,\, i-1} \right)}{z_{i} -z_{i-1} } \right]

.. math::
   :label: 7.119) 

   \frac{\partial q_{i} }{\partial \theta _{liq,\, i} } =-\left[\frac{k\left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \frac{\partial \psi _{i} }{\partial \theta _{liq,\, i} } \right]-\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i} } \left[\frac{\left(\psi _{i} -\psi _{i+1} \right)+\left(\psi _{E,\, i+1} -\psi _{E,\, i} \right)}{z_{i+1} -z_{i} } \right]

.. math::
   :label: 7.120) 

   \frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} } =\left[\frac{k\left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \frac{\partial \psi _{i+1} }{\partial \theta _{liq,\, i+1} } \right]-\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i+1} } \left[\frac{\left(\psi _{i} -\psi _{i+1} \right)+\left(\psi _{E,\, i+1} -\psi _{E,\, i} \right)}{z_{i+1} -z_{i} } \right].

The derivatives of the soil matric potential at the node depth are
derived from equation

.. math::
   :label: 7.121) 

   \frac{\partial \psi _{i-1} }{\partial \theta _{liq,\, \, i-1} } =-B_{i-1} \frac{\psi _{i-1} }{\theta _{\, \, i-1} }

.. math::
   :label: 7.122) 

   \frac{\partial \psi _{i} }{\partial \theta _{\, liq,\, i} } =-B_{i} \frac{\psi _{i} }{\theta _{i} }

.. math::
   :label: 7.123) 

   \frac{\partial \psi _{i+1} }{\partial \theta _{liq,\, i+1} } =-B_{i+1} \frac{\psi _{i+1} }{\theta _{\, i+1} }

with the constraint
:math:`0.01\, \theta _{sat,\, i} \le \theta _{\, i} \le \theta _{sat,\, i}` .

The derivatives of the hydraulic conductivity at the layer interface are
derived from equation

.. math::
   :label: 7.124) 

   \begin{array}{l} {\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i-1} } =\frac{\partial k\left[z_{h,\, i-1} \right]}{\partial \theta _{liq,\, i} } =\left(1-\frac{f_{frz,\, i-1} +f_{frz,\, i} }{2} \right)\left(2B_{i-1} +3\right)k_{sat} \left[z_{h,\, i-1} \right]\times } \\ {\qquad \left[\frac{0.5\left(\theta _{\, i-1} +\theta _{\, i} \right)}{0.5\left(\theta _{sat,\, i-1} +\theta _{sat,\, i} \right)} \right]^{2B_{i-1} +2} \left(\frac{0.5}{0.5\left(\theta _{sat,\, i-1} +\theta _{sat,\, i} \right)} \right)} \end{array}

.. math::
   :label: 7.125) 

   \begin{array}{l} {\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i} } =\frac{\partial k\left[z_{h,\, i} \right]}{\partial \theta _{liq,\, i+1} } =\left(1-\frac{f_{frz,\, i} +f_{frz,\, i+1} }{2} \right)\left(2B_{i} +3\right)k_{sat} \left[z_{h,\, i} \right]\times } \\ {\qquad \left[\frac{0.5\left(\theta _{\, i} +\theta _{\, i+1} \right)}{0.5\left(\theta _{sat,\, i} +\theta _{sat,\, i+1} \right)} \right]^{2B_{i} +2} \left(\frac{0.5}{0.5\left(\theta _{sat,\, i} +\theta _{sat,\, i+1} \right)} \right)} \end{array}.

Equilibrium soil matric potential and volumetric moisture
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The equilibrium soil matric potential :math:`\psi _{E}`  can be derived
from equation as

.. math::
   :label: ZEqnNum899028 

   \psi _{E} =\psi _{sat} \left(\frac{\theta _{E} \left(z\right)}{\theta _{sat} } \right)^{-B}

and the equilibrium volumetric water content
:math:`\theta _{E} \left(z\right)` at depth :math:`z` can also be
derived as

.. math::
   :label: 7.127) 

   \theta _{E} \left(z\right)=\theta _{sat} \left(\frac{\psi _{sat} +z_{\nabla } -z}{\psi _{sat} } \right)^{-\frac{1}{B} } .

Here, the soil matric potentials, the water table depth
:math:`z_{\nabla }`  and the soil depths have units of mm. For the
finite-difference scheme, a layer-average equilibrium volumetric water
content is used in equation and can be obtained from

.. math::
   :label: 7.128) 

   \overline{\theta _{E,\, i} }=\int _{z_{h,\, i-1} }^{z_{h,\, i} }\frac{\theta _{E} \left(z\right)}{z_{h,\, i} -z_{h,\, i-1} }  \, dz

which when integrated yields

.. math::
   :label: ZEqnNum445442 

   \overline{\theta _{E,\, i} }=\frac{\theta _{sat,\, i} \psi _{sat,\, i} }{\left(z_{h,\, i} -z_{h,\, i-1} \right)\left(1-\frac{1}{B_{i} } \right)} \left[\left(\frac{\psi _{sat,\, i} -z_{\nabla } +z_{h,\, i} }{\psi _{sat,\, i} } \right)^{1-\frac{1}{B_{i} } } -\left(\frac{\psi _{sat,\, i} -z_{\nabla } +z_{h,\, i-1} }{\psi _{sat,\, i} } \right)^{1-\frac{1}{B_{i} } } \right].

Equation is valid when the water table :math:`z_{\nabla }`  is deeper
than both interface depths :math:`z_{h,\, i-1}`  and
:math:`z_{h,\, i}` . Since the water table can be within the soil
column, the equation is modified if the water table is within soil layer
:math:`i` (:math:`z_{h,\, i-1} <z_{\nabla } <z_{h,\, i}` ) as a weighted
average of the saturated part and the unsaturated part

.. math::
   :label: 7.130) 

   \overline{\theta _{E,\, i} }=\overline{\theta _{E,\, sat,\, i} }\left(\frac{z_{h,\, i} -z_{\nabla } }{z_{h,\, i} -z_{h,\, i-1} } \right)+\overline{\theta _{E,\, unsat,\, i} }\left(\frac{z_{\nabla } -z_{h,\, i-1} }{z_{h,\, i} -z_{h,\, i-1} } \right)

where :math:`\overline{\theta _{E,\, sat,\, i} }=\theta _{sat,\, i}` 
and the unsaturated part :math:`\overline{\theta _{E,\, unsat,\, i} }`
is

.. math::
   :label: 7.131) 

   \overline{\theta _{E,\, unsat,\, i} }=\frac{\theta _{sat,\, i} \psi _{sat,\, i} }{\left(z_{\nabla } -z_{h,\, i-1} \right)\left(1-\frac{1}{B_{i} } \right)} \left[1-\left(\frac{\psi _{sat,\, i} -z_{\nabla } +z_{h,\, i-1} }{\psi _{sat,\, i} } \right)^{1-\frac{1}{B_{i} } } \right].

If :math:`z_{\nabla } <z_{h,\, i-1}` , then
:math:`\overline{\theta _{E,\, i} }=\overline{\theta _{E,\, sat,\, i} }=\theta _{sat,\, i}` .
If the water table is below the soil column
(:math:`z_{\nabla } >z_{h,\, N_{levsoi} }` ), an equilibrium volumetric
soil moisture is calculated for a virtual layer :math:`i=N_{levsoi} +1`
as

.. math::
   :label: ZEqnNum235293 

   \overline{\theta _{E,\, i=N_{levsoi+1} } }=\frac{\theta _{sat,i-1} \psi _{sat,\, i-1} }{\left(z_{\nabla } -z_{h,\, i-1} \right)\left(1-\frac{1}{B_{i-1} } \right)} \left[1-\left(\frac{\psi _{sat,\, i-1} -z_{\nabla } +z_{h,\, i-1} }{\psi _{sat,\, i-1} } \right)^{1-\frac{1}{B_{i-1} } } \right]

The equilibrium volumetric soil moisture is constrained by

.. math::
   :label: 7.133) 

   0\le \overline{\theta _{E,\, i} }\le \theta _{sat,\, i}

The equilibrium soil matric potential is then

.. math::
   :label: ZEqnNum533842 

   \psi _{E,\, i} =\psi _{sat,\, i} \left(\frac{\overline{\theta _{E,\, i} }}{\theta _{sat,\, i} } \right)^{-B_{i} } \ge -1\times 10^{8} \qquad \frac{\overline{\theta _{E,\, i} }}{\theta _{sat,\, i} } \ge 0.01

Equation set for layer :math:`i=1`
''''''''''''''''''''''''''''''''''''''''''

For the top soil layer (:math:`i=1`), the boundary condition is the
infiltration rate (section 7.3),
:math:`q_{i-1}^{n+1} =-q_{infl}^{n+1}` , and the water balance equation
is

.. math::
   :label: 7.135) 

   \frac{\Delta z_{i} \Delta \theta _{liq,\, i} }{\Delta t} =q_{infl}^{n+1} +q_{i}^{n+1} -e_{i} .

After grouping like terms, the coefficients of the tridiagonal set of
equations for :math:`i=1` are

.. math::
   :label: 7.136) 

   a_{i} =0

.. math::
   :label: 7.137) 

   b_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.138) 

   c_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} }

.. math::
   :label: 7.139) 

   r_{i} =q_{infl}^{n+1} -q_{i}^{n} +e_{i} .

Equation set for layers :math:`i=2,\ldots ,N_{levsoi} -1`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The coefficients of the tridiagonal set of equations for
:math:`i=2,\ldots ,N_{levsoi} -1` are

.. math::
   :label: 7.140) 

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} }

.. math::
   :label: 7.141) 

   b_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.142) 

   c_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} }

.. math::
   :label: 7.143) 

   r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .

Equation set for layers :math:`i=N_{levsoi} ,\ldots N_{levsoi} +1`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

For the lowest soil layer (:math:`i=N_{levsoi}` ), the bottom boundary
condition depends on the depth of the water table. If the water table is
within the soil column (:math:`z_{\nabla } \le z_{h,\, N_{levsoi} }` ),
a zero-flux bottom boundary condition is applied (:math:`q_{i}^{n} =0`)
and the coefficients of the tridiagonal set of equations for
:math:`i=N_{levsoi}`  are

.. math::
   :label: 7.144) 

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} }

.. math::
   :label: 7.145) 

   b_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.146) 

   c_{i} =0

.. math::
   :label: 7.147) 

   r_{i} =q_{i-1}^{n} +e_{i} .

The coefficients for the aquifer layer :math:`i=N_{levsoi} +1` are then

.. math::
   :label: 7.148) 

   a_{i} =0

.. math::
   :label: 7.149) 

   b_{i} =-\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.150) 

   c_{i} =0

.. math::
   :label: 7.151) 

   r_{i} =0.

If the water table is below the soil column
(:math:`z_{\nabla } >z_{h,\, N_{levsoi} }` ), the coefficients for
:math:`i=N_{levsoi}`  are

.. math::
   :label: 7.152) 

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} }

.. math::
   :label: 7.153) 

   b_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i} } -\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.154) 

   c_{i} =\frac{\partial q_{i} }{\partial \theta _{liq,\, i+1} }

.. math::
   :label: 7.155) 

   r_{i} =q_{i-1}^{n} -q_{i}^{n} +e_{i} .

The :math:`i=N_{levsoi} +1` terms are evaluated using

.. math::
   :label: 7.156) 

   \psi _{N_{levsoi} +1} =\psi _{sat,\, N_{levsoi} } \left[s_{N_{levsoi} +1} \right]^{-B_{N_{levsoi} } } \ge -1\times 10^{8}

.. math::
   :label: 7.157) 

   z_{N_{levsoi} +1} =0.5\left(z_{\nabla } +z_{N_{levsoi} } \right)

where

.. math::
   :label: 7.158) 

   s_{N_{levsoi} +1} =0.5\left(\frac{\theta _{sat,\, N_{levsoi} } +\theta _{N_{levsoi} } }{\theta _{sat,\, N_{levsoi} } } \right)\qquad 0.01\le s_{N_{levsoi} +1} \le 1,

:math:`\psi _{E,\, N_{levsoi} +1}`  is evaluated from equations and ,

and

.. math::
   :label: 7.159) 

   \frac{\partial \psi _{N_{levsoi} +1} }{\partial \theta _{liq,\, N_{levsoi} +1} } =-B_{N_{levsoi} } \frac{\psi _{N_{levsoi} +1} }{s_{\, N_{levsoi} } \theta _{sat,\, N_{levsoi} } } .

The coefficients for the aquifer layer :math:`i=N_{levsoi} +1` are then

.. math::
   :label: 7.160) 

   a_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i-1} }

.. math::
   :label: 7.161) 

   b_{i} =-\frac{\partial q_{i-1} }{\partial \theta _{liq,\, i} } -\frac{\Delta z_{i} }{\Delta t}

.. math::
   :label: 7.162) 

   c_{i} =0

.. math::
   :label: 7.163) 

   r_{i} =q_{i-1}^{n} .

Upon solution of the tridiagonal equation set (Press et al. 1992), the
liquid water contents are updated as follows

.. math::
   :label: 7.164) 

   w_{liq,\, i}^{n+1} =w_{liq,\, i}^{n} +\Delta \theta _{liq,\, i} \Delta z_{i} \qquad i=1,\ldots ,N_{levsoi} .

The volumetric water content is

.. math::
   :label: 7.165) 

   \theta _{i} =\frac{w_{liq,\, i} }{\Delta z_{i} \rho _{liq} } +\frac{w_{ice,\, i} }{\Delta z_{i} \rho _{ice} } .

.. _Frozen Soils and Perched Water Table:

Frozen Soils and Perched Water Table
----------------------------------------

When soils freeze, the power-law form of the ice impedance factor
(section 7.4.1) can greatly decrease the hydraulic conductivity of the
soil, leading to nearly impermeable soil layers. When unfrozen soil
layers are present above relatively ice-rich frozen layers, the
possibility exists for perched saturated zones. Lateral drainage from
perched saturated regions is parameterized as a function of the
thickness of the saturated zone

.. math::
   :label: 7.166) 

   q_{drai,perch} =k_{drai,\, perch} \left(z_{frost} -z_{\nabla ,perch} \right)

where :math:`k_{drai,\, perch}`  depends on topographic slope and soil
hydraulic conductivity,

.. math::
   :label: 7.167) 

   k_{drai,\, perch} =10^{-5} \sin (\beta )\left(\frac{\sum _{i=N_{perch} }^{i=N_{frost} }\Theta _{ice,i} k_{sat} \left[z_{i} \right]\Delta z_{i}  }{\sum _{i=N_{perch} }^{i=N_{frost} }\Delta z_{i}  } \right)

where :math:`\Theta _{ice}`  is an ice impedance factor determined from
the ice content of the soil layers interacting with the water table
(section 7.6), :math:`\beta`  is the mean grid cell topographic slope in
radians, :math:`z_{frost}` \ is the depth to the frost table, and
:math:`z_{\nabla ,perch}`  is the depth to the perched saturated zone.
The frost table :math:`z_{frost}`  is defined as the shallowest frozen
layer having an unfrozen layer above it, while the perched water table
:math:`z_{\nabla ,perch}`  is defined as the depth at which the
volumetric water content drops below a specified threshold. The default
threshold is set to 0.9. Drainage from the perched saturated zone
:math:`q_{drai,perch}`  is removed from layers :math:`N_{perch}` 
through :math:`N_{frost}` , which are the layers containing
:math:`z_{\nabla ,perch}`  and, :math:`z_{frost}` \ respectively.

.. _Groundwater-Soil Water Interactions:

Groundwater-Soil Water Interactions
---------------------------------------

Drainage or sub-surface runoff is based on the SIMTOP scheme (Niu et al.
2005) with a modification to account for reduced drainage in frozen
soils. In the work of Niu et al. (2005), the drainage :math:`q_{drai}` 
(kg m\ :sup:`-2` s\ :sup:`-1`) was formulated as

.. math::
   :label: ZEqnNum924767 

   q_{drai} =q_{drai,\, \max } \exp \left(-f_{drai} z_{\nabla } \right).

Here, the water table depth :math:`z_{\nabla }`  has units of meters. To
restrict drainage in frozen soils, Niu et al. (2005) added the following
condition

.. math::
   :label: 7.169) 

   q_{drai} =0\qquad {\rm for\; \; }w_{ice,\, N_{levsoi} } >w_{liq,\, N_{levsoi} } .

In preliminary testing it was found that a more gradual restriction of
drainage was required so that the water table depth remained dynamic
under partially frozen conditions. The following modification is made to
equation

.. math::
   :label: ZEqnNum150955 

   q_{drai} =\Theta _{ice} q_{drai,\, \max } \exp \left(-f_{drai} z_{\nabla } \right)

where :math:`\Theta _{ice}`  is an ice impedance factor determined from
the ice content of the soil layers interacting with the water table

.. math::
   :label: 7.171) 

   \Theta _{ice} =10^{-\Omega \left(\frac{\sum _{i=jwt}^{i=N_{levsoi} }F_{ice,i} \Delta z_{i}  }{\sum _{i=jwt}^{i=N_{levsoi} }\Delta z_{i}  } \right)}

where :math:`\Omega =6`\ is an adjustable parameter, :math:`jwt` is the
index of the layer directly above the water table,
:math:`F_{ice} =\frac{\theta _{ice} }{\theta _{sat} }`  is the
ice-filled fraction of the pore space of soil layer :math:`i` (kg
m\ :sup:`-2`), and :math:`\Delta z_{i}`  is the layer thickness
(mm). This expression is functionally the same as that used to determine
the ice impedance factor in section 7.4. In equation , the decay factor
:math:`f_{drai} =2.5` m\ :sup:`-1` and the maximum drainage when
the water table depth is at the surface
:math:`q_{drai,\, \max } =10\sin (\beta )` kg m\ :sup:`-2`
s\ :sup:`-1` , where :math:`\beta`  is the mean grid cell
topographic slope in radians, were determined for global simulations
through sensitivity analysis and comparison with observed runoff.

Determination of water table depth :math:`z_{\nabla }`  is based on work
by Niu et al. (2007). In this approach, a groundwater component is added
in the form of an unconfined aquifer lying below the soil column (Figure
7.1). The groundwater solution is dependent on whether the water table
is within or below the soil column. The water stored in the unconfined
aquifer :math:`W_{a}`  has a prescribed maximum value (5000 mm). When
the water table is within the soil column, :math:`W_{a}`  is constant
because there is no water exchange between the soil column and the
underlying aquifer. In this case, recharge to the water table is
diagnosed by applying Darcy’s law across the water table

.. math::
   :label: 7.172) 

   q_{rech\arg e} =-k_{aq} \frac{\left(\Psi _{\nabla } -\Psi _{jwt} \right)}{\left(z_{\nabla } -z_{jwt} \right)}

where :math:`\Psi _{\nabla } =0` is the matric potential at the water
table and\ :math:`k_{aq} =\Theta _{ice,jwt+1} k\left[z_{jwt+1} \right]`
is the hydraulic conductivity of the layer containing the water table.
Change in the water table is then calculated as the difference between
recharge and drainage, scaled by the specific yield of the layer
containing the water table

.. math::
   :label: ZEqnNum287831 

   \Delta z_{\nabla } =\frac{\left(q_{rech\arg e} -q_{drai} \right)}{S_{y} } \Delta t.

The specific yield, :math:`S_{y}` , which depends on the soil
properties and the water table location, is derived by taking the
difference between two equilibrium soil moisture profiles whose water
tables differ by an infinitesimal amount

.. math::
   :label: 7.174) 

   S_{y} =\theta _{sat} \left(1-\left(1+\frac{z_{\nabla } }{\Psi _{sat} } \right)^{\frac{-1}{B} } \right)

where B is the Clapp-Hornberger exponent. Because :math:`S_{y}`  is a
function of the soil properties, it results in water table dynamics that
are consistent with the soil water fluxes described in section 7.4.

For the case when the water table is below the soil column, the change
in water stored in the unconfined aquifer :math:`W_{a}`  (mm) is updated
as

.. math::
   :label: 7.174a) 

   \Delta W_{a}^{} =\left(q_{recharge} -q_{drai} \right)\Delta t

and the water table is updated using equation with the specific yield of
layer :math:`N_{levsoi}` .

The recharge rate is defined as positive when water enters the aquifer

.. math::
   :label: 7.174b) 

   q_{recharge} =\frac{\Delta \theta _{liq,\, N_{levsoi} +1} \Delta z_{N_{levsoi} +1} }{\Delta t}

where
:math:`\Delta \theta _{liq,\, N_{levsoi} +1} =\theta _{liq,\, N_{levsoi} +1}^{n+1} -\theta _{liq,\, N_{levsoi} +1}^{n}` 
is the change in liquid water content for layer :math:`i=N_{levsoi} +1`
calculated from the solution of the soil water equations (section 7.4),
and :math:`\Delta z_{N_{levsoi} +1}`  (mm) is

.. math::
   :label: 7.174c) 

   \Delta z_{N_{levsoi} +1} =z_{\nabla }^{n} -z_{h,\, N_{levsoi} } .

After the above calculations, two numerical adjustments are implemented
to keep the liquid water content of each soil layer
(:math:`w_{liq,\, i}` ) within physical constraints of
:math:`w_{liq}^{\min } \le w_{liq,\, i} \le \left(\theta _{sat,\, i} -\theta _{ice,\, i} \right)\Delta z_{i}` 
where :math:`w_{liq}^{\min } =0.01` (mm). First, beginning with the
bottom soil layer :math:`i=N_{levsoi}` , any excess liquid water in each
soil layer
(:math:`w_{liq,\, i}^{excess} =w_{liq,\, i} -\left(\theta _{sat,\, i} -\theta _{ice,\, i} \right)\Delta z_{i} \ge 0`)
is successively added to the layer above. Any excess liquid water that
remains after saturating the entire soil column (plus a maximum surface
ponding depth :math:`w_{liq}^{pond} =10` kg m\ :sup:`-2`), is
added to drainage :math:`q_{drai}` . Second, to prevent negative
:math:`w_{liq,\, i}` , each layer is successively brought up to
:math:`w_{liq,\, i} =w_{liq}^{\min }`  by taking the required amount of
water from the layer below. If this results in
:math:`w_{liq,\, N_{levsoi} } <w_{liq}^{\min }` , then the layers above
are searched in succession for the required amount of water
(:math:`w_{liq}^{\min } -w_{liq,\, N_{levsoi} }` ) and removed from
those layers subject to the constraint
:math:`w_{liq,\, i} \ge w_{liq}^{\min }` . If sufficient water is not
found, then the water is removed from :math:`W_{t}`  and
:math:`q_{drai}` .

The soil surface layer liquid water and ice contents are then updated
for dew :math:`q_{sdew}` , frost :math:`q_{frost}` , or sublimation
:math:`q_{subl}`  (section 5.4) as

.. math::
   :label: 7.175) 

   w_{liq,\, 1}^{n+1} =w_{liq,\, 1}^{n} +q_{sdew} \Delta t

.. math::
   :label: 7.176) 

   w_{ice,\, 1}^{n+1} =w_{ice,\, 1}^{n} +q_{frost} \Delta t

.. math::
   :label: 7.177) 

   w_{ice,\, 1}^{n+1} =w_{ice,\, 1}^{n} -q_{subl} \Delta t.

Sublimation of ice is limited to the amount of ice available.

.. _Runoff from glaciers and snow-capped surfaces:

Runoff from glaciers and snow-capped surfaces
-------------------------------------------------

All surfaces are constrained to have a snow water equivalent
:math:`W_{sno} \le 1000` kg m\ :sup:`-2`. For snow-capped
surfaces, the solid and liquid precipitation reaching the snow surface
and dew in solid or liquid form, is separated into solid
:math:`q_{snwcp,ice}` \ and liquid :math:`q_{snwcp,liq}`  runoff terms

.. math::
   :label: 7.178) 

   q_{snwcp,ice} =q_{grnd,ice} +q_{frost}

.. math::
   :label: 7.179) 

   q_{snwcp,liq} =q_{grnd,liq} +q_{dew}

and snow pack properties are unchanged. The :math:`q_{snwcp,ice}` 
runoff is sent to the River Transport Model (RTM) (Chapter 11) where it
is routed to the ocean as an ice stream and, if applicable, the ice is
melted there.

For snow-capped surfaces other than glaciers and lakes the
:math:`q_{snwcp,liq}`  runoff is assigned to the glaciers and lakes
runoff term :math:`q_{rgwl}`  (e.g. :math:`q_{rgwl} =q_{snwcp,liq}` ).
For glacier surfaces the runoff term :math:`q_{rgwl}`  is calculated
from the residual of the water balance

.. math::
   :label: 7.180) 

   q_{rgwl} =q_{grnd,ice} +q_{grnd,liq} -E_{g} -E_{v} -\frac{\left(W_{b}^{n+1} -W_{b}^{n} \right)}{\Delta t} -q_{snwcp,ice}

where :math:`W_{b}^{n}`  and :math:`W_{b}^{n+1}`  are the water balances
at the beginning and ending of the time step defined as

.. math::
   :label: 7.181) 

   W_{b} =W_{can} +W_{sno} +\sum _{i=1}^{N}\left(w_{ice,i} +w_{liq,i} \right) .

Currently, glaciers are non-vegetated and :math:`E_{v} =W_{can} =0`.
The contribution of lake runoff to :math:`q_{rgwl}`  is described in
section 9.6.3. The runoff term :math:`q_{rgwl}`  may be negative for
glaciers and lakes, which reduces the total amount of runoff available
to the RTM.

.. _The Variable Infiltration Capacity parameterizations as a hydrologic option:

The Variable Infiltration Capacity parameterizations as a hydrologic option
-------------------------------------------------------------------------------

The hydrologic parameterizations from the Variable Infiltration Capacity
(VIC) land surface model (Liang et al. 1994) have been implemented as a
hydrologic option. VIC includes two different time scales of runoff
generation. To capture such dynamics, the soil column in the original
VIC model
(http://www.hydro.washington.edu/Lettenmaier/Models/VIC/Overview/ModelOverview.shtml)
is typically divided into three layers with variable soil depths. The
upper two layers are designed to represent the dynamic responses of the
soil to rainfall events for surface runoff generation, and the lower
layer is used to characterize the seasonal soil moisture behavior and
subsurface runoff generation. The implementation of the VIC
parameterizations are as described in Li et al. (2011) except where
updated for consistency with modifications to CLM hydrology in CLM4.5.
Note that unless explicitly mentioned in this section, any descriptions
from sections 7.1-7.7 are intact and remain valid when the VIC hydrology
option is turned on.

Three VIC soil layers are defined by aggregating multiple layers in the
CLM soil column with thicknesses of :math:`\sum^3_{i=1}{\Delta z_i}`,
:math:`\sum^6_{i=4}{\Delta z_i}`,
:math:`\sum^{N_{levsoi}}_{i=7}{\Delta z_i}`, respectively. At each time
step, the soil moisture profile is determined following the algorithms
detailed in section 7.4, and aggregated to the three VIC layers for
runoff generation calculations. The surface runoff generated by the
saturation excess runoff mechanism, q\ :sub:`over`, is
calculated using equation , but with the fractional saturated area
defined as

.. math::
   :label: ZEqnNum627546 

   f_{sat} =1-\left(1-{w_{top} \mathord{\left/ {\vphantom {w_{top}  w_{m,\, top} }} \right. \kern-\nulldelimiterspace} w_{m,\, top} } \right)^{{1\mathord{\left/ {\vphantom {1 \left(1+b_{inf} \right)}} \right. \kern-\nulldelimiterspace} \left(1+b_{inf} \right)} }

where :math:`w_{top}` and :math:`w_{m,top}` are calculated as
:math:`\sum^6_{i=1}{{\theta }_i\Delta z_i}` and
:math:`\sum^6_{i=1}{{\theta }_{s,i}\Delta z_i}`, respectively, and
represent the soil moisture (kg m\ :sup:`-2`) and maximum soil
moisture (kg m\ :sup:`-2`) in the top two VIC layers combined.

In equation , it is hypothesized that the spatial heterogeneity of soil
moisture holding capacity in the top VIC layers can be represented by a
soil moisture holding capacity curve as defined in equation , in which
:math:`b_{inf}` is a parameter that controls the shape of the curve.
That is, if one assumes that a grid cell (or soil column in this case)
is composed of many pixels (or points) with varying soil moisture
capacity, this variation across the grid cell can be represented
conceptually as

.. math::
   :label: ZEqnNum929519 

   i=i_m\left(1-{\left(1-A\right)}^{1/b_{inf}}\right)

where :math:`i` and :math:`i_{m}`  are the point and maximum point soil
moisture holding capacities (kg m\ :sup:`-2`), respectively;
:math:`A` is the fraction of a grid cell for which the soil moisture
holding capacity is less than or equal to :math:`i`; and
:math:`i_m=w_{m,top}\left(1+b_{inf}\right)`. When :math:`A` is equal to
:math:`f_{sat}`, the corresponding point soil moisture holding capacity
is denoted as :math:`i_0`. The maximum soil infiltration capacity (kg
m\ :sup:`-2` s\ :sup:`-1`) in equation becomes

.. math::
   :label: ZEqnNum202398 

   q_{inf,\, max} =\left\{\begin{array}{l} {\left(1-f_{sat} \right)\theta _{ice,\, top} \left(\frac{w_{m,\, top} -w_{top} }{\Delta t} \right)\qquad i_{o} +q_{in,\, soil} \Delta t\ge i_{m} } \\ {\left(1-f_{sat} \right)\frac{\theta _{ice,\, top} }{\Delta t} \left\{\begin{array}{l} {\left(w_{m,\, top} -w_{top} \right)-w_{m,\, top} \times } \\ {\left[1-\max \left(1,\, \frac{\left(i_{o} +q_{in,\, soil} \Delta t\right)}{i_{m} } \right)^{1+b_{inf} } \right]} \end{array}\right\}\qquad i_{o} +q_{in,\, soil} \Delta t<i_{m} } \end{array}\right\}

where :math:`\theta _{ice,\, top}`  is an ice impedance factor
determined from the ice content of the top two VIC layers combined,
similar to the one used in equation . Interested readers are referred to
Wood et al. (1992) for a schematic representation of equation and
derivations associated with equations , , and .

The subsurface runoff in equation is parameterized as

.. math::
   :label: ZEqnNum411417 

   q_{drai}={\mathrm{\Theta }}_{ice,bot}\left[ \begin{array}{c} \frac{D_sD_{smax}}{W_sw_{m,bot}}w_{bot}+ \\ 
   max\left(0,\frac{w_{bot}-W_sw_{m,bot}}{{w_{m,bot}-W}_sw_{m,bot}}\right)\left(D_{smax}-\frac{D_sD_{smax}}{W_s}\right) \end{array}
   \right]/\mathrm{\Delta }t

where :math:`w_{bot}` and :math:`w_{m,bot}` are the soil moisture (kg
m\ :sup:`-2`) and maximum soil moisture (kg m\ :sup:`-2`) in
the bottom VIC layer, respectively, :math:`D_{smax}` is the maximum
subsurface flow rate (kg m\ :sup:`-2` s\ :sup:`-1`),
:math:`D_s` is a fraction of :math:`D_{smax}`, :math:`W_s` is a fraction
of :math:`w_{m,bot}`, and :math:`{\mathrm{\Theta }}_{ice,bot}` is an ice
impedance factor determined from the ice content of the bottom VIC
layer, similar to the ones in equations and .

As the VIC parameterizations are based on conceptual models, Huang and
Liang (2006) recommended calibrating the VIC parameters, including
:math:`b_{inf}` , :math:`D_{smax}` , :math:`D_s` , :math:`W_s` , and the
second and third layer soil thicknesses using observations. In this
implementation, the thicknesses of the VIC soil layers are fixed to
maintain consistency with the soil water algorithms in section 7.4. The
other four parameters, :math:`b_{inf}` , :math:`D_{smax}` , :math:`D_s`
, and :math:`W_s` are prescribed and are included in the CLM surface
dataset. Users can provide calibrated parameter values determined
manually or automatically by modifying the surface dataset. Note that
the units of :math:`D_{smax}` on the surface dataset are mm
d\ :sup:`-1` (the traditional units for other standard VIC
applications) which are then converted to kg m\ :sup:`-2`
s\ :sup:`-1` for use in CLM. A preliminary calibration was
performed by perturbing the three parameters :math:`b_{inf}` ,
:math:`D_{smax}` , and :math:`W_s`, and fixing :math:`D_s=0.1` globally.
The parameter space for :math:`b_{inf}` , :math:`D_{smax}` , and
:math:`W_s` was sampled using the global sensitivity analysis framework
described by Hou et al. (2012) to produce 64 combinations of parameter
values based on *a priori* information about the parameters. For each
set of parameter values, a global simulation was performed using the
compset I\_2000 (i.e., driven by satellite phenology) at a resolution of
0.9\ :sup:`o`\ x1.25\ :sup:`o` on the basis of the
development tag betr\_m\_sci10\_clm45sci13\_clm4\_0\_54. At each model
grid cell, the set of :math:`b_{inf}` , :math:`D_{smax}` , and
:math:`W_s` values corresponding to the simulation that produced the
lowest absolute bias compared to the climatological mean annual total
runoff from the Global Runoff Data Center (GRDC) was selected as the
calibrated values. These values are provided only as a reference due to
the preliminary nature of the calibration. Interested users of the VIC
hydrology option are encouraged to calibrate the parameters for their
applications for improved performance.

