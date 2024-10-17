.. _rst_River Transport Model (RTM):

Model for Scale Adaptive River Transport (MOSART)
=================================================

.. _Overview MOSART:

Overview
---------

MOSART is a river transport model designed for applications across local, regional and global scales :ref:`(Li et al., 2013b) <Lietal2013b>`. A major purpose of MOSART is to provide freshwater input for the ocean model in coupled Earth System Models. MOSART also provides an effective way of evaluating and diagnosing the soil hydrology simulated by land surface models through direct comparison of the simulated river flow with observations of natural streamflow at gauging stations :ref:`(Li et al., 2015a)<Lietal2015a>`. Moreover, MOSART provides a modeling framework for representing riverine transport and transformation of energy and biogeochemical fluxes under both natural and human-influenced conditions ( :ref:`(Li et al., 2015b) <Lietal2015b>`.

.. _Routing Processes:

Routing Processes
------------------

MOSART divides each spatial unit such as a lat/lon grid or watershed into three categories of hydrologic units (as shown in :numref:`Figure MOSART conceptual diagram`): hillslopes that convert both surface and subsurface runoff into tributaries, tributaries that discharge into a single main channel, and the main channel that connects the local spatial unit with upstream/downstream units through the river network. MOSART assumes that all the tributaries within a spatial unit can be treated as a single hypothetical sub-network channel with a transport capacity equivalent to all the tributaries combined. Correspondingly, three routing processes are represented in MOSART: 1) hillslope routing: in each spatial unit, surface runoff is routed as overland flow into the sub-network channel, while subsurface runoff generated in the spatial unit directly enters the sub-network channel; 2) sub-network channel routing: the sub-network channel receives water from the hillslopes, routes water through the channel and discharges it into the main channel; 3) main channel routing: the main channel receives water from the sub-network channel and/or inflow, if any, from the upstream spatial units, and discharges the water to its downstream spatial unit or the ocean.

.. Figure 14.1. MOSART conceptual diagram

.. _Figure MOSART conceptual diagram:

.. figure:: mosart_diagram.png
    :width: 800px
    :height: 400px

MOSART only routes positive runoff, although negative runoff can be generated occasionally by the land model (e.g., :math:`q_{gwl}`). Negative runoff in any runoff component including :math:`q_{sur}`, :math:`q_{sub}`, :math:`q_{gwl}` is not routed through MOSART, but instead is mapped directly from the spatial unit where it is generated at any time step to the coupler.

In MOSART, the travel velocities of water across hillslopes, sub-network and main channel are all estimated using Manning's equation with different levels of simplifications. Generally the Manning's equation is in the form of

.. math::
   :label: 14.1

   V = \frac{R^{\frac{2}{3}} S_{f}}{n}

where :math:`V` is the travel velocity (m s :sup:`-1` ), :math:`R` is the hydraulic radius (m). :math:`S_{f}` is the friction slope that accounts for the effects of gravity, friction, inertia and other forces on the water. If the channel slope is steep enough, the gravity force dominates over the others so one can approximate :math:`S_{f}` by the channel bed slope :math:`S`, which is the key assumption underpinning the kinematic wave method. :math:`n` is the Manning's roughness coefficient, which is mainly controlled by surface roughness and sinuosity of the flow path.

If the water surface is sufficiently large or the water depth :math:`h` is sufficiently shallow, the hydraulic radius can be approximated by the water depth. This is the case for both hillslope and sub-network channel routing.

.. math::
   :label: 14.2

   R_{h} = h_{h}
   R_{t} = h_{t}

Here :math:`R_{h}` (m) and :math:`R_{t}` (m) are hydraulic radius for hillslope and sub-network channel routing respectively, and :math:`h_{h}` (m) and :math:`h_{t}` (m) are water depth during hillslope and sub-network channel routing respectively.

For the main channel, the hydraulic radius is given by

.. math::
   :label: 14.3

   R_{r} = \frac{A_{r}}{P_{r}}

where :math:`A_{r}` (m :sup:`2` ) is the wetted area defined as the part of the channel cross-section area below the water surface, :math:`P_{r}` (m) is the wetted perimeter, the perimeter confined in the wetted area.

For hillslopes, sub-network and main channels, a common continuity equation can be written as

.. math::
   :label: 14.4

   \frac{dS}{dt} = Q_{in} - Q_{out} + R

where :math:`Q_{in}` (m :sup:`3` s :sup:`-1` ) is the main channel flow from the upstream grid(s) into the main channel of the current grid, which is zero for hillslope and sub-network routing. :math:`Q_{out}` (m :sup:`3` s :sup:`-1` ) is the outflow rate from hillslope into the sub-network, from the sub-network into the main channel, or from the current main channel to the main channel of its downstream grid (if not the outlet grid) or ocean (if the current grid is the basin outlet). :math:`R` (m :sup:`3` s :sup:`-1` ) is a source term, which could be the surface runoff generation rate for hillslopes, or lateral inflow (from hillslopes) into sub-network channel or water-atmosphere exchange fluxes such as precipitation and evaporation. It is assumed that surface runoff is generated uniformly across all the hillslopes. Currently, MOSART does not exchange water with the atmosphere or return water to the land model so its function is strictly to transport water from runoff generation through the hillslope, tributaries, and main channels to the basin outlets.

.. _Numerical Solution MOSART:

Numerical Solution
----------------------------

The numerical implementation of MOSART is mainly based on a subcycling scheme and a local time-stepping algorithm. There are two levels of subcycling. For convenience, we denote :math:`T_{inputs}` (s), :math:`T_{mosart}` (s), :math:`T_{hillslope}` (s) and :math:`T_{channel}` (s) as the time steps of runoff inputs (from CLM to MOSART via the flux coupler), MOSART routing, hillslope routing, and channel routing, respectively. The first level of subcycling is between the runoff inputs and MOSART routing. If :math:`T_{inputs}` is 10800s and :math:`T_{mosart}` is 3600s, three MOSART time steps will be invoked each time the runoff inputs are updated. The second level of subcycling is between the hillslope routing and channel routing. This is to account for the fact that the travel velocity of water across hillslopes is usually much slower than that in the channels. :math:`T_{hillslope}` is usually set as the same as :math:`T_{mosart}`, but within each time step of hillslope routing there are a few time steps for channel routing, i.e., :math:`T_{hillslope} = D_{levelH2R} \cdot T_{channel}`. The local time-stepping algorithm is to account for the fact that the travel velocity of water is much faster in some river channels (e.g., with steeper bed slope, narrower channel width) than others. That is, for each channel (either a sub-network or main channel), the final time step of local channel routing is given as :math:`T_{local}=T_{channel}/D_{local}`. :math:`D_{local}` is currently estimated empirically as a function of local channel slope, width, length and upstream drainage area. If MOSART crashes due to a numerical issue, we recommend increasing :math:`D_{levelH2R}` and, if the issue remains, reducing :math:`T_{mosart}`.

.. _Parameters and Input Data:

Parameters and Input Data
---------------------------------

MOSART is supported by a comprehensive, global hydrography dataset at 0.5 ° resolution. As such, the fundamental spatial unit of MOSART is a 0.5 ° lat/lon grid. The topographic parameters (such as flow direction, channel length, topographic and channel slopes, etc.) were derived using the Dominant River Tracing (DRT) algorithm (:ref:`Wu et al., 2011<Wuetal2011>`; :ref:`Wu et al. 2012 <Wuetal2012>`). The DRT algorithm produces the topographic parameters in a scale-consistent way to preserve/upscale the key features of a baseline high-resolution hydrography dataset at multiple coarser spatial resolutions. Here the baseline high-resolution hydrography dataset is the 1km resolution Hydrological data and maps based on SHuttle Elevation Derivatives at multiple Scales (HydroSHEDS) (:ref:`Lehner and Döll, 2004 <LehnerDoll2004>`; :ref:`Lehner et al., 2008 <Lehneretal2008>`). The channel geometry parameters, e.g., bankfull width and depth, were estimated from empirical hydraulic geometry relationships as functions of the mean annual discharge. The Manning roughness coefficients for overland and channel flow were calculated as functions of landcover and water depth. For more details on the methodology to derive channel geometry and the Manning's roughness coefficients, please refer to :ref:`Getirana et al. (2012) <Getiranaetal2012>`. The full list of parameters included in this global hydrography dataset is provided in :numref:`Table MOSART Parameters`. Evaluation of global simulations by MOSART using the aforementioned parameters is described in :ref:`Li et al. (2015b) <Lietal2015b>`.

.. _Table MOSART Parameters:

.. table:: List of parameters in the global hydrography dataset

 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | Name                    | Unit          | Description                                                                                                                        |
 +=========================+===============+====================================================================================================================================+
 | :math:`F_{dir}`         | \-            | The D8 single flow direction for each coarse grid cell coded using 1 (E), 2 (SE), 4 (S), 8 (SW), 16 (W), 32 (NW), 64 (N), 128 (NE) |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`A_{total}`       | km :sup:`2`   | The upstream drainage area of each coarse grid cell                                                                                |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`F_{dis}`         | m             | The dominant river length for each coarse grid cell                                                                                |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`S_{channel}`     | \-            | The average channel slope for each coarse grid cell                                                                                |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`S_{topographic}` | \-            | The average topographic slope (for overland flow routing) for each coarse grid cell                                                |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`A_{local}`       | km :sup:`2`   | The surface area for each coarse grid cell                                                                                         |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`D_{p}`           | m :sup:`-1`   | Drainage density, calculated  as the total channel length within each coarse grid cell divided by the local cell area              |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`D_{r}`           | m             | The bankfull depth of main channel                                                                                                 |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`W_{r}`           | m             | The bankfull width of main channel                                                                                                 |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`D_{t}`           | m             | The average bankfull depth of tributary channels                                                                                   |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`W_{t}`           | m             | The average bankfull width of tributary channels                                                                                   |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`n_{r}`           | \-            | Manning's roughness coefficient for channel flow routing                                                                           |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+
 | :math:`n_{h}`           | \-            | Manning's roughness coefficient for overland flow routing                                                                          |
 +-------------------------+---------------+------------------------------------------------------------------------------------------------------------------------------------+

Difference between CLM5.0 and CLM4.5
-------------------------------------

1. Routing methods: RTM, a linear reservoir method, is used in CLM4.5 for river routing, whilst in CLM5.0, MOSART is an added option for river routing based on the more physically-based kinematic wave method.

2. Runoff treatment: In RTM runoff is routed regardless of its sign so negative streamflow can be simulated at times. MOSART routes only non-negative runoff and always produces positive streamflow, which is important for future extensions to model riverine heat and biogeochemical fluxes.

3. Input parameters: RTM in CLM4.5 only requires one layer of a spatially varying variable of channel velocity, whilst MOSART in CLM5.0 requires 13 parameters that are all available globally at 0.5 ° resolution.

4. Outputs: RTM only produces streamflow simulation, whilst MOSART additionally simulates the time-varying channel velocities, channel water depth, and channel surface water variations.

