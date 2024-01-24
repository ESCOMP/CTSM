.. _rst_Soil and Snow Temperatures:

Soil and Snow Temperatures
=============================

The first law of heat conduction is

.. math::
   :label: 6.1

   F=-\lambda \nabla T

where :math:`F` is the amount of heat conducted across a unit cross-sectional area in unit time (W m\ :sup:`-2`), :math:`\lambda` is thermal conductivity (W m\ :sup:`-1` K\ :sup:`-1`), and :math:`\nabla T` is the spatial gradient of temperature (K m\ :sup:`-1`). In one-dimensional form

.. math::
   :label: 6.2

   F_{z} =-\lambda \frac{\partial T}{\partial z}

where :math:`z` is in the vertical direction (m) and is positive downward and :math:`F_{z}` is positive upward. To account for non-steady or transient conditions, the principle of energy conservation in the form of the continuity equation is invoked as

.. math::
   :label: 6.3

   c\frac{\partial T}{\partial t} =-\frac{\partial F_{z} }{\partial z}

where :math:`c` is the volumetric snow/soil heat capacity (J m\ :sup:`-3` K\ :sup:`-1`) and :math:`t` is time (s). Combining equations and yields the second law of heat conduction in one-dimensional form

.. math::
   :label: 6.4

   c\frac{\partial T}{\partial t} =\frac{\partial }{\partial z} \left[\lambda \frac{\partial T}{\partial z} \right].

This equation is solved numerically to calculate the soil, snow, and surface water temperatures for a 25-layer soil column with up to twelve overlying layers of snow and a single surface water layer with the boundary conditions of :math:`h` as the heat flux into the top soil, snow, and surface water layers from the overlying atmosphere (section :numref:`Numerical Solution Temperature`) and zero heat flux at the bottom of the soil column. The temperature profile is calculated first without phase change and then readjusted for phase change (section :numref:`Phase Change`).

.. _Numerical Solution Temperature:

Numerical Solution
----------------------

The soil column is discretized into 25 layers (section :numref:`Vertical Discretization`) where :math:`N_{levgrnd} = 25` is the number of soil layers (:numref:`Table Soil layer structure`).

The overlying snow pack is modeled with up to twelve layers depending on the total snow depth. The layers from top to bottom are indexed in the Fortran code as :math:`i=-4,-3,-2,-1,0`, which permits the accumulation or ablation of snow at the top of the snow pack without renumbering the layers. Layer :math:`i=0` is the snow layer next to the soil surface and layer :math:`i=snl+1` is the top layer, where the variable :math:`snl` is the negative of the number of snow layers. The number of snow layers and the thickness of each layer is a function of snow depth :math:`z_{sno}` (m) as follows.

.. math::

   \left\{ \begin{array}{l}
   snl=-1 \\
   \Delta z_{0} = z_{sno}
   \end{array} \right\} & \qquad {\rm for\; 0.01}\le {\rm z}_{{\rm sno}} \le 0.03 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-2  \\
   \Delta z_{-1} ={z_{sno} \mathord{\left/ {\vphantom {z_{sno}  2}} \right.} 2} \\
   \Delta z_{0} = \Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.03}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.04 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-2 \\
   \Delta z_{-1} = 0.02 \\
   \Delta z_{0} = z_{sno} -\Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.04}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.07 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-3 \\
   \Delta z_{-2} = 0.02 \\
   \Delta z_{-1} = {\left(z_{sno} -0.02\right)\mathord{\left/ {\vphantom {\left(z_{sno} -0.02\right) 2}} \right.} 2} \\
   \Delta z_{0} = \Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.07}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.12 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-3 \\
   \Delta z_{-2} = 0.02 \\
   \Delta z_{-1} = 0.05 \\
   \Delta z_{0} = z_{sno} -\Delta z_{-2} -\Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.12}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.18 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-4  \\
   \Delta z_{-3} = 0.02 \\
   \Delta z_{-2} = 0.05 \\
   \Delta z_{-1} = {\left(z_{sno} -\Delta z_{-3} -\Delta z_{-2} \right)\mathord{\left/ {\vphantom {\left(z_{sno} -\Delta z_{-3} -\Delta z_{-2} \right) 2}} \right.} 2}  \\
   \Delta z_{0} =\Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.18}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.29 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-4 \\
   \Delta z_{-3} = 0.02  \\
   \Delta z_{-2} = 0.05  \\
   \Delta z_{-1} = 0.11  \\
   \Delta z_{0} = z_{sno} -\Delta z_{-3} -\Delta z_{-2} -\Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.29}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.41 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-5  \\
   \Delta z_{-4} = 0.02  \\
   \Delta z_{-3} = 0.05  \\
   \Delta z_{-2} = 0.11  \\
   \Delta z_{-1} = {\left(z_{sno} -\Delta z_{-4} -\Delta z_{-3} -\Delta z_{-2} \right)\mathord{\left/ {\vphantom {\left(z_{sno} -\Delta z_{-4} -\Delta z_{-3} -\Delta z_{-2} \right) 2}} \right.} 2}  \\
   \Delta z_{0} = \Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.41}\, {\rm <}\, {\rm z}_{{\rm sno}} \le 0.64 \\

.. math::

   \left\{ \begin{array}{l}
   snl=-5 \\
   \Delta z_{-4} = 0.02  \\
   \Delta z_{-3} = 0.05  \\
   \Delta z_{-2} = 0.11  \\
   \Delta z_{-1} = 0.23  \\
   \Delta z_{0} = z_{sno} -\Delta z_{-4} -\Delta z_{-3} -\Delta z_{-2} -\Delta z_{-1}
   \end{array} \right\} & \qquad {\rm for\; 0.64}\, {\rm <}\, {\rm z}_{{\rm sno}}

The node depths, which are located at the midpoint of the snow layers, and the layer interfaces are both referenced from the soil surface and are defined as negative values

.. math::
   :label: 6.8

   z_{i} =z_{h,\, i} -0.5\Delta z_{i} \qquad i=snl+1,\ldots ,0

.. math::
   :label: 6.9

   z_{h,\, i} =z_{h,\, i+1} -\Delta z_{i+1} \qquad i=snl,\ldots ,-1.

Note that :math:`z_{h,\, 0}`, the interface between the bottom snow layer and the top soil layer, is zero. Thermal properties (i.e., temperature :math:`T_{i}` [K]; thermal conductivity :math:`\lambda _{i}` [W m\ :sup:`-1` K\ :sup:`-1`]; volumetric heat capacity :math:`c_{i}` [J m\ :sup:`-3` K\ :sup:`-1`]) are defined for soil layers at the node depths (:numref:`Figure Soil Temperature Schematic`) and for snow layers at the layer midpoints. When present, snow occupies a fraction of a grid cell's area, therefore snow depth represents the thickness of the snowpack averaged over only the snow covered area. The grid cell average snow depth is related to the depth of the snow covered area as :math:`\bar{z}_{sno} =f_{sno} z_{sno}`. By default, the grid cell average snow depth is written to the history file.

The heat flux :math:`F_{i}`  (W m\ :sup:`-2`) from layer :math:`i`
to layer :math:`i+1` is

.. math::
   :label: 6.10

   F_{i} =-\lambda \left[z_{h,\, i} \right]\left(\frac{T_{i} -T_{i+1} }{z_{i+1} -z_{i} } \right)

where the thermal conductivity at the interface :math:`\lambda \left[z_{h,\, i} \right]` is

.. math::
   :label: 6.11

   \lambda \left[z_{h,\, i} \right]=\left\{\begin{array}{l} {\frac{\lambda _{i} \lambda _{i+1} \left(z_{i+1} -z_{i} \right)}{\lambda _{i} \left(z_{i+1} -z_{h,\, i} \right)+\lambda _{i+1} \left(z_{h,\, i} -z_{i} \right)} \qquad i=snl+1,\ldots ,N_{levgrnd} -1} \\ {0\qquad i=N_{levgrnd} } \end{array}\right\}.

These equations are derived, with reference to :numref:`Figure Soil Temperature Schematic`, assuming that the heat flux from :math:`i` (depth :math:`z_{i}` ) to the interface between :math:`i` and :math:`i+1` (depth :math:`z_{h,\, i}` ) equals the heat flux from the interface to :math:`i+1` (depth :math:`z_{i+1}` ), i.e.,

.. math::
   :label: 6.12

   -\lambda _{i} \frac{T_{i} -T_{m} }{z_{h,\, i} -z_{i} } =-\lambda _{i+1} \frac{T_{m} -T_{i+1} }{z_{i+1} -z_{h,\, i} }

where :math:`T_{m}` is the temperature at the interface of layers :math:`i` and :math:`i+1`.

Shown are three soil layers, :math:`i-1`, :math:`i`, and :math:`i+1`. The thermal conductivity :math:`\lambda`, specific heat capacity :math:`c`, and temperature :math:`T` are defined at the layer node depth :math:`z`. :math:`T_{m}` is the interface temperature. The thermal conductivity :math:`\lambda \left[z_{h} \right]` is defined at the interface of two layers :math:`z_{h}`. The layer thickness is :math:`\Delta z`. The heat fluxes :math:`F_{i-1}` and :math:`F_{i}` are defined as positive upwards.

.. _Figure Soil Temperature Schematic:

.. figure:: image1.png

 Schematic diagram of numerical scheme used to solve for soil temperature.

The energy balance for the :math:`i^{th}`  layer is

.. math::
   :label: 6.13

   \frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)=-F_{i-1} +F_{i}

where the superscripts :math:`n` and :math:`n+1` indicate values at the beginning and end of the time step, respectively, and :math:`\Delta t` is the time step (s). This equation is solved using the Crank-Nicholson method, which combines the explicit method with fluxes evaluated at :math:`n` (:math:`F_{i-1}^{n},F_{i}^{n}` ) and the implicit method with fluxes evaluated at :math:`n+1` (:math:`F_{i-1}^{n+1},F_{i}^{n+1}` )

.. math::
   :label: 6.14

   \frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)=\alpha \left(-F_{i-1}^{n} +F_{i}^{n} \right)+\left(1-\alpha \right)\left(-F_{i-1}^{n+1} +F_{i}^{n+1} \right)

where :math:`\alpha =0.5`, resulting in a tridiagonal system of equations

.. math::
   :label: 6.15

   r_{i} =a_{i} T_{i-1}^{n+1} +b_{i} T_{i}^{n+1} +c_{i} T_{i+1}^{n+1}

where :math:`a_{i}`, :math:`b_{i}`, and :math:`c_{i}` are the subdiagonal, diagonal, and superdiagonal elements in the tridiagonal matrix and :math:`r_{i}` is a column vector of constants. When surface water is present, the equation for the top soil layer has an additional term representing the surface water temperature; this results in a four element band-diagonal system of equations.

For the top soil layer :math:`i=1`, top snow layer :math:`i=snl+1`, or surface water layer, the heat flux from the overlying atmosphere :math:`h` (W m\ :sup:`-2`, defined as positive into the surface) is

.. math::
   :label: 6.16

   h^{n+1} =-\alpha F_{i-1}^{n} -\left(1-\alpha \right)F_{i-1}^{n+1} .

The energy balance for these layers is then

.. math::
   :label: 6.17

   \frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)=h^{n+1} +\alpha F_{i}^{n} +\left(1-\alpha \right)F_{i}^{n+1} .

The heat flux :math:`h` at :math:`n+1` may be approximated as follows

.. math::
   :label: 6.18

   h^{n+1} =h^{n} +\frac{\partial h}{\partial T_{i} } \left(T_{i}^{n+1} -T_{i}^{n} \right).

The resulting equations are then

.. math::
   :label: 6.19

   \begin{array}{rcl} {\frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)} & {=} & {h^{n} +\frac{\partial h}{\partial T_{i} } \left(T_{i}^{n+1} -T_{i} \right)} \\ {} & {} & {-\alpha \frac{\lambda \left[z_{h,\, i} \right]\left(T_{i}^{n} -T_{i+1}^{n} \right)}{z_{i+1} -z_{i} } -\left(1-\alpha \right)\frac{\lambda \left[z_{h,\, i} \right]\left(T_{i}^{n+1} -T_{i+1}^{n+1} \right)}{z_{i+1} -z_{i} } } \end{array}

For the top snow layer, :math:`i=snl+1`, the coefficients are

.. math::
   :label: 6.20

   a_{i} =0

.. math::
   :label: 6.21

   b_{i} =1+\frac{\Delta t}{c_{i} \Delta z_{i} } \left[\left(1-\alpha \right)\frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } -\frac{\partial h}{\partial T_{i} } \right]

.. math::
   :label: 6.22

   c_{i} =-\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} }

.. math::
   :label: 6.23

   r_{i} =T_{i}^{n} +\frac{\Delta t}{c_{i} \Delta z_{i} } \left[h_{sno} ^{n} -\frac{\partial h}{\partial T_{i} } T_{i}^{n} +\alpha F_{i} \right]

where

.. math::
   :label: 6.24

   F_{i} =-\lambda \left[z_{h,\, i} \right]\left(\frac{T_{i}^{n} -T_{i+1}^{n} }{z_{i+1} -z_{i} } \right).

The heat flux into the snow surface from the overlying atmosphere :math:`h` is

.. math::
   :label: 6.25

   h=\overrightarrow{S}_{sno} -\overrightarrow{L}_{sno} -H_{sno} -\lambda E_{sno}

where :math:`\overrightarrow{S}_{sno}` is the solar radiation absorbed by the top snow layer (section :numref:`Snow Albedo`), :math:`\overrightarrow{L}_{sno}` is the longwave radiation absorbed by the snow (positive toward the atmosphere) (section :numref:`Longwave Fluxes`), :math:`H_{sno}` is the sensible heat flux from the snow (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`), and :math:`\lambda E_{sno}` is the latent heat flux from the snow (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`). The partial derivative of the heat flux :math:`h` with respect to temperature is

.. math::
   :label: 6.26

   \frac{\partial h}{\partial T_{} } =-\frac{\partial \overrightarrow{L}_{} }{\partial T_{} } -\frac{\partial H_{} }{\partial T_{} } -\frac{\partial \lambda E_{} }{\partial T_{} }

where the partial derivative of the net longwave radiation is

.. math::
   :label: 6.27

   \frac{\partial \overrightarrow{L}_{} }{\partial T_{} } =4\varepsilon _{g} \sigma \left(T_{}^{n} \right)^{3}

and the partial derivatives of the sensible and latent heat fluxes are given by equations and for non-vegetated surfaces, and by equations and for vegetated surfaces. :math:`\sigma` is the Stefan-Boltzmann constant (W m\ :sup:`-2` K\ :sup:`-4`) (:numref:`Table Physical Constants`) and :math:`\varepsilon _{g}` is the ground emissivity (section :numref:`Longwave Fluxes`). For purposes of computing :math:`h` and :math:`\frac{\partial h}{\partial T_{g} }`, the term :math:`\lambda` is arbitrarily assumed to be

.. math::
   :label: 6.28

   \lambda =\left\{\begin{array}{l} {\lambda _{sub} \qquad {\rm if\; }w_{liq,\, snl+1} =0{\rm \; and\; }w_{ice,\, snl+1} >0} \\ {\lambda _{vap} \qquad {\rm otherwise}} \end{array}\right\}

where :math:`\lambda _{sub}` and :math:`\lambda _{vap}` are the latent heat of sublimation and vaporization, respectively (J kg\ :sup:`-1`) (:numref:`Table Physical Constants`), and :math:`w_{liq,\, snl+1}` and :math:`w_{ice,\, snl+1}` are the liquid water and ice contents of the top snow/soil layer, respectively (kg m\ :sup:`-2`) (Chapter :numref:`rst_Hydrology`).

For the top soil layer, :math:`i=1`, the coefficients are

.. math::
   :label: 6.29

   a_{i} =-f_{sno} \left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} }

.. math::
   :label: 6.30

   b_{i} =1+\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \left[f_{sno} \frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } +\frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \right]-\left(1-f_{sno} \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\partial h}{\partial T}

.. math::
   :label: 6.31

   c_{i} =-\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} }

.. math::
   :label: 6.32

   r_{i} =T_{i}^{n} +\frac{\Delta t}{c_{i} \Delta z_{i} } \left[\left(1-f_{sno} \right)\left(h_{soil} ^{n} -\frac{\partial h}{\partial T_{} } T_{i}^{n} \right)+\alpha \left(F_{i} -f_{sno} F_{i-1} \right)\right]

The heat flux into the soil surface from the overlying atmosphere
:math:`h` is

.. math::
   :label: 6.33

   h=\overrightarrow{S}_{soil} -\overrightarrow{L}_{soil} -H_{soil} -\lambda E_{soil}

It can be seen that when no snow is present (:math:`f_{sno} =0`), the expressions for the coefficients of the top soil layer have the same form as those for the top snow layer.

The surface snow/soil layer temperature computed in this way is the layer-averaged temperature and hence has somewhat reduced diurnal amplitude compared with surface temperature. An accurate surface temperature is provided that compensates for this effect and numerical error by tuning the heat capacity of the top layer (through adjustment of the layer thickness) to give an exact match to the analytic solution for diurnal heating. The top layer thickness for :math:`i=snl+1` is given by

.. math::
   :label: 6.34

   \Delta z_{i*} =0.5\left[z_{i} -z_{h,\, i-1} +c_{a} \left(z_{i+1} -z_{h,\, i-1} \right)\right]

where :math:`c_{a}` is a tunable parameter, varying from 0 to 1, and is taken as 0.34 by comparing the numerical solution with the analytic solution (:ref:`Z.-L. Yang 1998, unpublished manuscript<Yang1998>`). :math:`\Delta z_{i*}` is used in place of :math:`\Delta z_{i}` for :math:`i=snl+1` in equations -. The top snow/soil layer temperature computed in this way is the ground surface temperature :math:`T_{g}^{n+1}`.

The boundary condition at the bottom of the snow/soil column is zero heat flux, :math:`F_{i} =0`, resulting in, for :math:`i=N_{levgrnd}`,

.. math::
   :label: 6.35

   \frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)=\alpha \frac{\lambda \left[z_{h,\, i-1} \right]\left(T_{i-1}^{n} -T_{i}^{n} \right)}{z_{i} -z_{i-1} } +\left(1-\alpha \right)\frac{\lambda \left[z_{h,\, i-1} \right]\left(T_{i-1}^{n+1} -T_{i}^{n+1} \right)}{z_{i} -z_{i-1} }

.. math::
   :label: 6.36

   a_{i} =-\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} }

.. math::
   :label: 6.37

   b_{i} =1+\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} }

.. math::
   :label: 6.38

   c_{i} =0

.. math::
   :label: 6.39

   r_{i} =T_{i}^{n} -\alpha \frac{\Delta t}{c_{i} \Delta z_{i} } F_{i-1}

where

.. math::
   :label: 6.40

   F_{i-1} =-\frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } \left(T_{i-1}^{n} -T_{i}^{n} \right).

For the interior snow/soil layers, :math:`snl+1<i<N_{levgrnd}`, excluding the top soil layer,

.. math::
   :label: 6.41

   \begin{array}{rcl} {\frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{i}^{n+1} -T_{i}^{n} \right)} & {=} & {-\alpha \frac{\lambda \left[z_{h,\, i} \right]\left(T_{i}^{n} -T_{i+1}^{n} \right)}{z_{i+1} -z_{i} } +\alpha \frac{\lambda \left[z_{h,\, i-1} \right]\left(T_{i-1}^{n} -T_{i}^{n} \right)}{z_{i} -z_{i-1} } } \\ {} \end{array}

.. math::
   :label: 6.42

   a_{i} =-\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} }

.. math::
   :label: 6.43

   b_{i} =1+\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \left[\frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } +\frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \right]

.. math::
   :label: 6.44

   c_{i} =-\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} }

.. math::
   :label: 6.45

   r_{i} =T_{i}^{n} +\alpha \frac{\Delta t}{c_{i} \Delta z_{i} } \left(F_{i} -F_{i-1} \right)+\frac{\Delta t}{c_{i} \Delta z_{i} } \vec{S}_{g,i} .

where :math:`\vec{S}_{g,i}` is the absorbed solar flux in layer :math:`i` (section :numref:`Snow Albedo`).

When surface water exists, the following top soil layer coefficients are modified

.. math::
   :label: 6.46

   \begin{array}{l} {b_{i} =1+\left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \left[f_{h2osfc} \frac{\lambda _{h2osfc} }{z_{i} -z_{h2osfc} } +f_{sno} \frac{\lambda \left[z_{h,\, i-1} \right]}{z_{i} -z_{i-1} } +\frac{\lambda \left[z_{h,\, i} \right]}{z_{i+1} -z_{i} } \right]} \\ {\quad \quad -\left(1-f_{sno} -f_{h2osfc} \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\partial h}{\partial T} } \end{array}

.. math::
   :label: 6.47

   r_{i} =T_{i}^{n} +\frac{\Delta t}{c_{i} \Delta z_{i} } \left[\begin{array}{l} {\left(1-f_{sno} -f_{h2osfc} \right)\left(h_{soil} ^{n} -\frac{\partial h}{\partial T_{} } T_{i}^{n} \right)} \\ {+\alpha \left(F_{i} -f_{sno} F_{i-1} +f_{h2osfc} \frac{\lambda _{h2osfc} }{z_{1} -z_{h2osfc} } \left(T_{1} -T_{h2osfc} \right)\right)} \end{array}\right]

.. math::
   :label: 6.48

   d_{i} =-f_{h2osfc} \left(1-\alpha \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \left[\frac{\lambda _{h2osfc} }{z_{i} -z_{h2osfc} } \right]

where :math:`d_{i}` is an additional coefficient representing the heat flux from the surface water layer. The surface water layer coefficients are

.. math::
   :label: 6.49

   a_{h2osfc} =0

.. math::
   :label: 6.50

   b_{h2osfc} =1+\frac{\Delta t}{c_{h2osfc} \Delta z_{h2osfc} } \left[\left(1-\alpha \right)\frac{\lambda _{h2osfc} }{z_{1} -z_{h2osfc} } -\frac{\partial h}{\partial T} \right]

.. math::
   :label: 6.51

   c_{h2osfc} =-\left(1-\alpha \right)\frac{\Delta t}{c_{h2osfc} \Delta z_{h2osfc} } \frac{\lambda _{h2osfc} }{z_{1} -z_{h2osfc} }

.. math::
   :label: 6.52

   r_{h2osfc} =T_{h2osfc}^{n} +\frac{\Delta t}{c_{i} \Delta z_{i} } \left[h_{h2osfc} ^{n} -\frac{\partial h}{\partial T_{} } T_{h2osfc}^{n} +\alpha \frac{\lambda _{h2osfc} }{z_{1} -z_{h2osfc} } \left(T_{1} -T_{h2osfc} \right)\right]_{}

.. _Phase Change:

Phase Change
----------------

.. _Soil and Snow Layers:

Soil and Snow Layers
^^^^^^^^^^^^^^^^^^^^^^^^^^

Upon update, the snow/soil temperatures are evaluated to determine if phase change will take place as

.. math::
   :label: 6.53a

   \begin{array}{lr}
   T_{i}^{n+1} >T_{f} {\rm \; and\; }w_{ice,\, i} >0 & \qquad i=snl+1,\ldots ,N_{levgrnd} \qquad {\rm melting}   \end{array}

.. math::
   :label: 6.53b

   \begin{array}{lr}
   \begin{array}{lr}
   T_{i}^{n+1} <T_{f} {\rm \; and\; }w_{liq,\, i} >0 & \qquad i=snl+1,\ldots ,0 \\
   T_{i}^{n+1} <T_{f} {\rm \; and\; }w_{liq,\, i} >w_{liq,\, \max ,\, i} & \quad i=1,\ldots ,N_{levgrnd}
   \end{array} & \quad {\rm freezing}
   \end{array}

where :math:`T_{i}^{n+1}` is the soil layer temperature after solution of the tridiagonal equation set, :math:`w_{ice,\, i}` and :math:`w_{liq,\, i}` are the mass of ice and liquid water (kg m\ :sup:`-2`) in each snow/soil layer, respectively, and :math:`T_{f}` is the freezing temperature of water (K) (:numref:`Table Physical Constants`). For the freezing process in soil layers, the concept of supercooled soil water from :ref:`Niu and Yang (2006)<NiuYang2006>` is adopted. The supercooled soil water is the liquid water that coexists with ice over a wide range of temperatures below freezing and is implemented through a freezing point depression equation

.. math::
   :label: 6.54

   w_{liq,\, \max ,\, i} =\Delta z_{i} \theta _{sat,\, i} \left[\frac{10^{3} L_{f} \left(T_{f} -T_{i} \right)}{gT_{i} \psi _{sat,\, i} } \right]^{{-1\mathord{\left/ {\vphantom {-1 B_{i} }} \right.} B_{i} } } \qquad T_{i} <T_{f}

where :math:`w_{liq,\, \max,\, i}` is the maximum liquid water in layer :math:`i` (kg m\ :sup:`-2`) when the soil temperature :math:`T_{i}` is below the freezing temperature :math:`T_{f}`, :math:`L_{f}` is the latent heat of fusion (J kg\ :sup:`-1`) (:numref:`Table Physical Constants`), :math:`g` is the gravitational acceleration (m s\ :sup:`-2`) (:numref:`Table Physical Constants`), and :math:`\psi _{sat,\, i}` and :math:`B_{i}` are the soil texture-dependent saturated matric potential (mm) and :ref:`Clapp and Hornberger (1978)<ClappHornberger1978>` exponent (section :numref:`Soil Water`).

For the special case when snow is present (snow mass :math:`W_{sno} >0`) but there are no explicit snow layers (:math:`snl=0`) (i.e., there is not enough snow present to meet the minimum snow depth requirement of 0.01 m), snow melt will take place for soil layer :math:`i=1` if the soil layer temperature is greater than the freezing temperature (:math:`T_{1}^{n+1} >T_{f}` ).

The rate of phase change is assessed from the energy excess (or deficit) needed to change :math:`T_{i}` to freezing temperature, :math:`T_{f}`. The excess or deficit of energy :math:`H_{i}` (W m\ :sup:`-2`) is determined as follows

.. math::
   :label: 6.55

   H_{i} =\left\{\begin{array}{lr}
   \frac{\partial h}{\partial T} \left(T_{f} -T_{i}^{n} \right)-\frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{f} -T_{i}^{n} \right) & \quad \quad i=snl+1 \\
   \left(1-f_{sno} -f_{h2osfc} \right)\frac{\partial h}{\partial T} \left(T_{f} -T_{i}^{n} \right)-\frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{f} -T_{i}^{n} \right)\quad {\kern 1pt} {\kern 1pt} {\kern 1pt} {\kern 1pt} & i=1 \\
   -\frac{c_{i} \Delta z_{i} }{\Delta t} \left(T_{f} -T_{i}^{n} \right) & \quad \quad i\ne \left\{1,snl+1\right\}
   \end{array}\right\}.

If the melting criteria is met :eq:`6.53a` and :math:`H_{m} =\frac{H_{i} \Delta t}{L_{f} } >0`, then the ice mass is readjusted as

.. math::
   :label: 6.56

   w_{ice,\, i}^{n+1} =w_{ice,\, i}^{n} -H_{m} \ge 0\qquad i=snl+1,\ldots ,N_{levgrnd} .

If the freezing criteria is met :eq:`6.53b` and :math:`H_{m} <0`, then the ice mass is readjusted for :math:`i=snl+1,\ldots,0` as

.. math::
   :label: 6.57

   w_{ice,\, i}^{n+1} =\min \left(w_{liq,\, i}^{n} +w_{ice,\, i}^{n} ,w_{ice,\, i}^{n} -H_{m} \right)

and for :math:`i=1,\ldots,N_{levgrnd}` as

.. math::
   :label: 6.58

   w_{ice,\, i}^{n+1} =
   \left\{\begin{array}{lr}
   \min \left(w_{liq,\, i}^{n} +w_{ice,\, i}^{n} -w_{liq,\, \max ,\, i}^{n} ,\, w_{ice,\, i}^{n} -H_{m} \right) & \qquad w_{liq,\, i}^{n} +w_{ice,\, i}^{n} \ge w_{liq,\, \max ,\, i}^{n} {\rm \; } \\
   {\rm 0} & \qquad w_{liq,\, i}^{n} +w_{ice,\, i}^{n} <w_{liq,\, \max ,\, i}^{n} {\rm \; \; }\,
   \end{array}\right\}.

Liquid water mass is readjusted as

.. math::
   :label: 6.59

   w_{liq,\, i}^{n+1} =w_{liq,\, i}^{n} +w_{ice,\, i}^{n} -w_{ice,\, i}^{n+1} \ge 0.

Because part of the energy :math:`H_{i}` may not be consumed in melting or released in freezing, the energy is recalculated as

.. math::
   :label: 6.60

   H_{i*} =H_{i} -\frac{L_{f} \left(w_{ice,\, i}^{n} -w_{ice,\, i}^{n+1} \right)}{\Delta t}

and this energy is used to cool or warm the snow/soil layer (if :math:`\left|H_{i*} \right|>0`) as

.. math::
   :label: 6.61

   T_{i}^{n+1} =
   \left\{\begin{array}{lr}
   T_{f} +{\frac{\Delta t}{c_{i} \Delta z_{i} } H_{i*} \mathord{\left/ {\vphantom {\frac{\Delta t}{c_{i} \Delta z_{i} } H_{i*}  \left(1-\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\partial h}{\partial T} \right)}} \right.} \left(1-\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\partial h}{\partial T} \right)} & \quad \quad \quad \quad \, i=snl+1 \\
   T_{f} +{\frac{\Delta t}{c_{i} \Delta z_{i} } H_{i*} \mathord{\left/ {\vphantom {\frac{\Delta t}{c_{i} \Delta z_{i} } H_{i*}  \left(1-\left(1-f_{sno} -f_{h2osfc} \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\partial h}{\partial T} \right)}} \right.} \left(1-\left(1-f_{sno} -f_{h2osfc} \right)\frac{\Delta t}{c_{i} \Delta z_{i} } \frac{\partial h}{\partial T} \right)} & \qquad i=1 \\
   T_{f} +\frac{\Delta t}{c_{i} \Delta z_{i} } H_{i*} & \quad \quad \quad \quad \, i\ne \left\{1,snl+1\right\}
   \end{array}\right\}.

For the special case when snow is present (:math:`W_{sno} >0`), there are no explicit snow layers (:math:`snl=0`), and :math:`\frac{H_{1} \Delta t}{L_{f} } >0` (melting), the snow mass :math:`W_{sno}` (kg m\ :sup:`-2`) is reduced according to

.. math::
   :label: 6.62

   W_{sno}^{n+1} =W_{sno}^{n} -\frac{H_{1} \Delta t}{L_{f} } \ge 0.

The snow depth is reduced proportionally

.. math::
   :label: 6.63

   z_{sno}^{n+1} =\frac{W_{sno}^{n+1} }{W_{sno}^{n} } z_{sno}^{n} .

Again, because part of the energy may not be consumed in melting, the energy for the surface soil layer :math:`i=1` is recalculated as

.. math::
   :label: 6.64

   H_{1*} =H_{1} -\frac{L_{f} \left(W_{sno}^{n} -W_{sno}^{n+1} \right)}{\Delta t} .

If there is excess energy (:math:`H_{1*} >0`), this energy becomes available to the top soil layer as

.. math::
   :label: 6.65

   H_{1} =H_{1*} .

The ice mass, liquid water content, and temperature of the top soil layer are then determined from :eq:`6.56`, :eq:`6.59`, and :eq:`6.61` using the recalculated energy from :eq:`6.65`. Snow melt :math:`M_{1S}` (kg m\ :sup:`-2` s\ :sup:`-1`) and phase change energy :math:`E_{p,\, 1S}` (W m\ :sup:`-2`) for this special case are

.. math::
   :label: 6.66

   M_{1S} =\frac{W_{sno}^{n} -W_{sno}^{n+1} }{\Delta t} \ge 0

.. math::
   :label: 6.67

   E_{p,\, 1S} =L_{f} M_{1S} .

The total energy of phase change :math:`E_{p}` (W m\ :sup:`-2`) for the snow/soil column is

.. math::
   :label: 6.68

   E_{p} =E_{p,\, 1S} +\sum _{i=snl+1}^{N_{levgrnd} }E_{p,i}

where

.. math::
   :label: 6.69

   E_{p,\, i} =L_{f} \frac{\left(w_{ice,\, i}^{n} -w_{ice,\, i}^{n+1} \right)}{\Delta t} .

The total snow melt :math:`M` (kg m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 6.70

   M=M_{1S} +\sum _{i=snl+1}^{i=0}M_{i}

where

.. math::
   :label: 6.71

   M_{i} =\frac{\left(w_{ice,\, i}^{n} -w_{ice,\, i}^{n+1} \right)}{\Delta t} \ge 0.

The solution for snow/soil temperatures conserves energy as

.. math::
   :label: 6.72

   G-E_{p} -\sum _{i=snl+1}^{i=N_{levgrnd} }\frac{c_{i} \Delta z_{i} }{\Delta t}  \left(T_{i}^{n+1} -T_{i}^{n} \right)=0

where :math:`G` is the ground heat flux (section :numref:`Update of Ground Sensible and Latent Heat Fluxes`).

.. _Surface Water:

Surface Water
^^^^^^^^^^^^^^^^^^^

Phase change of surface water takes place when the surface water temperature, :math:`T_{h2osfc}`, becomes less than :math:`T_{f}`. The energy available for freezing is

.. math::
   :label: 6.73

   H_{h2osfc} =\frac{\partial h}{\partial T} \left(T_{f} -T_{h2osfc}^{n} \right)-\frac{c_{h2osfc} \Delta z_{h2osfc} }{\Delta t} \left(T_{f} -T_{h2osfc}^{n} \right)

where :math:`c_{h2osfc}` is the volumetric heat capacity of water, and :math:`\Delta z_{h2osfc}` is the depth of the surface water layer. If :math:`H_{m} =\frac{H_{h2osfc} \Delta t}{L_{f} } >0` then :math:`H_{m}` is removed from surface water and added to the snow column as ice

.. math::
   :label: 6.74

   H^{n+1} _{h2osfc} =H^{n} _{h2osfc} -H_{m}

.. math::
   :label: 6.75

   w_{ice,\, 0}^{n+1} =w_{ice,\, 0}^{n} +H_{m}

The snow depth is adjusted to account for the additional ice mass

.. math::
   :label: 6.76

   \Delta z_{sno} =\frac{H_{m} }{\rho _{ice} }

If :math:`H_{m}` \ is greater than :math:`W_{sfc}`, the excess heat :math:`\frac{L_{f} \left(H_{m} -W_{sfc} \right)}{\Delta t}` is used to cool the snow layer.

.. _Soil and Snow Thermal Properties:

Soil and Snow Thermal Properties
------------------------------------

The thermal properties of the soil are assumed to be a weighted combination of the mineral and organic properties of the soil (:ref:`Lawrence and Slater 2008 <LawrenceSlater2008>`). The soil layer organic matter fraction :math:`f_{om,i}` is

.. math::
   :label: 6.77

   f_{om,i} =\rho _{om,i} /\rho _{om,\max } .

Soil thermal conductivity :math:`\lambda _{i}` (W m\ :sup:`-1` K\ :sup:`-1`) is from :ref:`Farouki (1981) <Farouki1981>`

.. math::
   :label: 6.78

   \begin{array}{lr}
   \lambda _{i} = \left\{
   \begin{array}{lr}
   K_{e,\, i} \lambda _{sat,\, i} +\left(1-K_{e,\, i} \right)\lambda _{dry,\, i} &\qquad S_{r,\, i} > 1\times 10^{-7}  \\
   \lambda _{dry,\, i} &\qquad S_{r,\, i} \le 1\times 10^{-7}
   \end{array}\right\} &\qquad i=1,\ldots ,N_{levsoi}  \\
   \lambda _{i} =\lambda _{bedrock} &\qquad i=N_{levsoi} +1,\ldots N_{levgrnd}
   \end{array}

where :math:`\lambda _{sat,\, i}` is the saturated thermal conductivity, :math:`\lambda _{dry,\, i}` is the dry thermal conductivity, :math:`K_{e,\, i}` is the Kersten number, :math:`S_{r,\, i}` is the wetness of the soil with respect to saturation, and :math:`\lambda _{bedrock} =3` W m\ :sup:`-1` K\ :sup:`-1` is the thermal conductivity assumed for the deep ground layers (typical of saturated granitic rock; :ref:`Clauser and Huenges 1995 <ClauserHuenges1995>`). For glaciers,

.. math::
   :label: 6.79

   \lambda _{i} =\left\{\begin{array}{l} {\lambda _{liq,\, i} \qquad T_{i} \ge T_{f} } \\ {\lambda _{ice,\, i} \qquad T_{i} <T_{f} } \end{array}\right\}

where :math:`\lambda _{liq}` and :math:`\lambda _{ice}` are the thermal conductivities of liquid water and ice, respectively (:numref:`Table Physical Constants`). The saturated thermal conductivity :math:`\lambda _{sat,\, i}` (W m\ :sup:`-1` K\ :sup:`-1`) depends on the thermal conductivities of the soil solid, liquid water, and ice constituents

.. math::
   :label: 6.80

   \lambda _{sat} =\lambda _{s}^{1-\theta _{sat} } \lambda _{liq}^{\frac{\theta _{liq} }{\theta _{liq} +\theta _{ice} } \theta _{sat} } \lambda _{ice}^{\theta _{sat} \left(1-\frac{\theta _{liq} }{\theta _{liq} +\theta _{ice} } \right)}

where the thermal conductivity of soil solids :math:`\lambda _{s,\, i}` varies with the sand, clay, and organic matter content

.. math::
   :label: 6.81

   \lambda _{s,i} =(1-f_{om,i} )\lambda _{s,\min ,i} +f_{om,i} \lambda _{s,om}

where the mineral soil solid thermal conductivity :math:`\lambda _{s,\min,i}` \ is

.. math::
   :label: 6.82

   \lambda _{s,\, \min ,i} =\frac{8.80{\rm \; }\left(\% sand\right)_{i} +{\rm 2.92\; }\left(\% clay\right)_{i} }{\left(\% sand\right)_{i} +\left(\% clay\right)_{i} } ,

and :math:`\lambda _{s,om} =0.25`\ W m\ :sup:`-1` K\ :sup:`-1` (:ref:`Farouki 1981 <Farouki1981>`). :math:`\theta _{sat,\, i}` is the volumetric water content at saturation (porosity) (section :numref:`Hydraulic Properties`).

The thermal conductivity of dry soil is

.. math::
   :label: 6.83

   \lambda _{dry,i} =(1-f_{om,i} )\lambda _{dry,\min ,i} +f_{om,i} \lambda _{dry,om}

where the thermal conductivity of dry mineral soil :math:`\lambda _{dry,\min,i}` \ (W m\ :sup:`-1` K\ :sup:`-1`) depends on the bulk density :math:`\rho _{d,\, i} =2700\left(1-\theta _{sat,\, i} \right)` (kg m\ :sup:`-3`) as

.. math::
   :label: 6.84

   \lambda _{dry,\, \min ,i} =\frac{0.135\rho _{d,\, i} +64.7}{2700-0.947\rho _{d,\, i} }

and :math:`\lambda _{dry,om} =0.05` W m\ :sup:`-1` K\ :sup:`-1` (:ref:`Farouki 1981 <Farouki1981>`) is the dry thermal conductivity of organic matter. The Kersten number :math:`K_{e,\, i}` is a function of the degree of saturation :math:`S_{r}` and phase of water

.. math::
   :label: 6.85

   K_{e,\, i} = \left\{
   \begin{array}{lr}
   \log \left(S_{r,\, i} \right)+1\ge 0 &\qquad T_{i} \ge T_{f}  \\
   S_{r,\, i} &\qquad T_{i} <T_{f}
   \end{array}\right\}

where

.. math::
   :label: 6.86

   S_{r,\, i} =\left(\frac{w_{liq,\, i} }{\rho _{liq} \Delta z_{i} } +\frac{w_{ice,\, i} }{\rho _{ice} \Delta z_{i} } \right)\frac{1}{\theta _{sat,\, i} } =\frac{\theta _{liq,\, i} +\theta _{ice,\, i} }{\theta _{sat,\, i} } \le 1.

Thermal conductivity :math:`\lambda _{i}` (W m\ :sup:`-1` K\ :sup:`-1`) for snow is from :ref:`Jordan (1991) <Jordan1991>`

.. math::
   :label: 6.87

   \lambda _{i} =\lambda _{air} +\left(7.75\times 10^{-5} \rho _{sno,\, i} +1.105\times 10^{-6} \rho _{sno,\, i}^{2} \right)\left(\lambda _{ice} -\lambda _{air} \right)

where :math:`\lambda _{air}` is the thermal conductivity of air (:numref:`Table Physical Constants`) and :math:`\rho _{sno,\, i}` is the bulk density of snow (kg m\ :sup:`-3`)

.. math::
   :label: 6.88

   \rho _{sno,\, i} =\frac{w_{ice,\, i} +w_{liq,\, i} }{\Delta z_{i} } .

The volumetric heat capacity :math:`c_{i}` (J m\ :sup:`-3` K\ :sup:`-1`) for soil is from :ref:`de Vries (1963) <deVries1963>` and depends on the heat capacities of the soil solid, liquid water, and ice constituents

.. math::
   :label: 6.89

   c_{i} =c_{s,\, i} \left(1-\theta _{sat,\, i} \right)+\frac{w_{ice,\, i} }{\Delta z_{i} } C_{ice} +\frac{w_{liq,\, i} }{\Delta z_{i} } C_{liq}

where :math:`C_{liq}` and :math:`C_{ice}` are the specific heat capacities (J kg\ :sup:`-1` K\ :sup:`-1`) of liquid water and ice, respectively (:numref:`Table Physical Constants`). The heat capacity of soil solids :math:`c_{s,i}` \ (J m\ :sup:`-3` K\ :sup:`-1`) is

.. math::
   :label: 6.90

   c_{s,i} =(1-f_{om,i} )c_{s,\min ,i} +f_{om,i} c_{s,om}

where the heat capacity of mineral soil solids :math:`c_{s,\min,\, i}` (J m\ :sup:`-3` K\ :sup:`-1`) is

.. math::
   :label: 6.91

   \begin{array}{lr}
   c_{s,\min ,\, i} =\left(\frac{2.128{\rm \; }\left(\% sand\right)_{i} +{\rm 2.385\; }\left(\% clay\right)_{i} }{\left(\% sand\right)_{i} +\left(\% clay\right)_{i} } \right)\times 10^{6} &\qquad i=1,\ldots ,N_{levsoi}  \\
   c_{s,\, \min ,i} =c_{s,\, bedrock} &\qquad i=N_{levsoi} +1,\ldots ,N_{levgrnd}
   \end{array}

where :math:`c_{s,bedrock} =2\times 10^{6}` J m\ :sup:`-3` K\ :sup:`-1` is the heat capacity of bedrock and :math:`c_{s,om} =2.5\times 10^{6}` \ J m\ :sup:`-3` K\ :sup:`-1` (:ref:`Farouki 1981 <Farouki1981>`) is the heat capacity of organic matter. For glaciers and snow

.. math::
   :label: 6.92

   c_{i} =\frac{w_{ice,\, i} }{\Delta z_{i} } C_{ice} +\frac{w_{liq,\, i} }{\Delta z_{i} } C_{liq} .

For the special case when snow is present (:math:`W_{sno} >0`) but there are no explicit snow layers (:math:`snl=0`), the heat capacity of the top layer is a blend of ice and soil heat capacity

.. math::
   :label: 6.93

   c_{1} =c_{1}^{*} +\frac{C_{ice} W_{sno} }{\Delta z_{1} }

where :math:`c_{1}^{*}`  is calculated from :eq:`6.89` or :eq:`6.92`.

.. _Excess Ground Ice:

Excess Ground Ice
------------------------------------

An optional parameterization of excess ground ice melt and respective subsidence based on (:ref:`Lee et al., (2014) <Leeetal2014>`). Initial excess ground ice concentrations for soil columns are derived from (:ref:`Brown et al., (1997) <Brownetal1997>`). When the excess ground ice is present in the soil column, soil depth for a given layer (:math:`z_{i}`) is adjusted by the amount of excess ice in the column:

.. math::
   :label: 6.94

   z_{i}^{'}=\Sigma_{j=1}^{i} \ z_{j}^{'}+\frac{w_{exice,\, j}}{\rho_{ice} }

where :math:`w_{exice,\,j}` is excess ground ice amount (kg m :sup:`-2`) in layer :math:`j` and :math:`\rho_{ice}` is the density of ice (kg m :sup:`-3`). After adjustment of layer depths have been made, all of the soil temperature equations (from :eq:`6.80` to :eq:`6.89`) are calculted based on the adjusted depths. Thermal properties are additionally adjusted (:eq:`6.8` and :eq:`6.8`) in the following way:

.. math::
   :label: 6.95

   \begin{array}{lr}
    \theta_{sat}^{'} =\frac{\theta _{liq} }{\theta _{liq} +\theta _{ice} +\theta_{exice}}{\theta_{sat}} \\
    \lambda _{sat}^{'} =\lambda _{s}^{1-\theta _{sat}^{'} } \lambda _{liq}^{\frac{\theta _{liq} }{\theta _{liq} +\theta _{ice} +\theta_{exice}} \theta _{sat}^{'} } \lambda _{ice}^{\theta _{sat}^{'} \left(1-\frac{\theta _{liq} }{\theta _{liq} +\theta _{ice} +\theta_{exice}} \right)} \\
    c_{i}^{'} =c_{s,\, i} \left(1-\theta _{sat,\, i}^{'} \right)+\frac{w_{ice,\, i} +w_{exice,\,j}}{\Delta z_{i}^{'} } C_{ice} +\frac{w_{liq,\, i} }{\Delta z_{i}^{'} } C_{liq}
   \end{array}

Soil subsidence at the timestep :math:`n+1` (:math:`z_{exice}^{n+1}`, m) is then calculated as:

.. math::
   :label: 6.96

   z_{exice}^{n+1}=\Sigma_{i=1}^{N_{levgrnd}} \ z_{j}^{',\ ,n+1}-z_{j}^{',\ ,n }

With regards to hydraulic counductivity, excess ground ice is treated the same way normal soil ice is treated in :numref:`Frozen Soils and Perched Water Table`. When a soil layer thaws, excess ground ice is only allowed to melt when no normals soil ice is present in the layer. When a soil layer refreezes, liquid soil water can only turn into normal soil ice, thus, no new of excess ice can be created but only melted. The excess liquid soil moisture from excess ice melt is distributed within the soil column according to :numref:`Lateral Sub-surface Runoff`.

