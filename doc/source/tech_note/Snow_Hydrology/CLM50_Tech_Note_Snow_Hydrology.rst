.. _rst_Snow Hydrology:

Snow Hydrology
===============

The parameterizations for snow are based primarily on :ref:`Anderson (1976) <Anderson1976>`, :ref:`Jordan (1991) <Jordan1991>`, and :ref:`Dai and Zeng (1997) <DaiZeng1997>`. The snowpack can have up to twelve layers. These layers are indexed in the Fortran code as :math:`i=-11,-10,...,-1,0` where layer :math:`i=0` is the snow layer next to the top soil layer and layer :math:`i=-11` is the top layer of a twelve-layer snow pack. Since the number of snow layers varies according to the snow depth, we use the notation :math:`snl+1` to describe the top layer of snow for the variable layer snow pack, where :math:`snl` is the negative of the number of snow layers. Refer to :numref:`Figure three layer snow pack` for an example of the snow layer structure for a three layer snow pack.

.. _Figure three layer snow pack:

.. Figure:: image1.png

 Example of three layer snow pack (:math:`snl=-3`).

Shown are three snow layers, :math:`i=-2`, :math:`i=-1`, and :math:`i=0`. The layer node depth is :math:`z`, the layer interface is :math:`z_{h}`, and the layer thickness is :math:`\Delta z`.

The state variables for snow are the mass of water :math:`w_{liq,i}` (kg m\ :sup:`-2`), mass of ice :math:`w_{ice,i}` (kg m\ :sup:`-2`), layer thickness :math:`\Delta z_{i}` (m), and temperature :math:`T_{i}` (Chapter :numref:`rst_Soil and Snow Temperatures`). The water vapor phase is neglected. Snow can also exist in the model without being represented by explicit snow layers. This occurs when the snowpack is less than a specified minimum snow depth (:math:`z_{sno} < 0.01` m). In this case, the state variable is the mass of snow :math:`W_{sno}` (kg m\ :sup:`-2`).

Section :numref:`Snow Covered Area Fraction` describes the calculation of fractional snow covered area, which is used in the surface albedo calculation (Chapter :numref:`rst_Surface Albedos`) and the surface flux calculations (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`). The following two sections (:numref:`Ice Content` and :numref:`Water Content`) describe the ice and water content of the snow pack assuming that at least one snow layer exists. Section :numref:`Black and organic carbon and mineral dust within snow` describes how black and organic carbon and mineral dust particles are represented within snow, including meltwater flushing. See Section :numref:`Initialization of snow layer` for a description of how a snow layer is initialized.

.. _Snow Covered Area Fraction:

Snow Covered Area Fraction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The fraction of the ground covered by snow, :math:`f_{sno}`, is based on the method of :ref:`Swenson and Lawrence (2012) <SwensonLawrence2012>`. Because the processes governing snowfall and snowmelt differ, changes in :math:`f_{sno}` are calculated separately for accumulation and depletion. When snowfall occurs, :math:`f_{sno}` is updated as

.. math::
   :label: 8.14

   f^{n+1} _{sno} =1-\left(\left(1-\tanh (k_{accum} q_{sno} \Delta t)\right)\left(1-f^{n} _{sno} \right)\right)

where :math:`k_{accum}` is a constant whose default value is 0.1, :math:`q_{sno} \Delta t` is the amount of new snow, :math:`f^{n+1} _{sno}` is the updated snow covered fraction (SCF), and :math:`f^{n} _{sno}` is the SCF from the previous time step.

When snow melt occurs, :math:`f_{sno}` is calculated from the depletion curve

.. math::
   :label: 8.15

   f_{sno} =1-\left(\frac{\cos ^{-1} \left(2R_{sno} -1\right)}{\pi } \right)^{N_{melt} }

where :math:`R_{sno}` is the ratio of :math:`W_{sno}` to the maximum accumulated snow :math:`W_{\max }`, and :math:`N_{melt}` is a parameter that depends on the topographic variability within the grid cell. Whenever :math:`W_{sno}` reaches zero, :math:`W_{\max }` is reset to zero. The depletion curve shape parameter is defined as

.. math::
   :label: 8.16

   N_{melt} =\frac{200}{\min \left(10,\sigma _{topo} \right)}

The standard deviation of the elevation within a grid cell, :math:`\sigma _{topo}`, is calculated from a high resolution DEM (a 1km DEM is used for CLM). Note that *glacier\_mec* columns (section :numref:`Multiple elevation class scheme`) are treated differently in this respect, as they already account for the subgrid topography in a grid cell in their own way. Therefore, in each *glacier\_mec* column very flat terrain is assumed, implemented as :math:`N_{melt}=10`.

.. _Ice Content:

Ice Content
^^^^^^^^^^^^^^^^^

The conservation equation for mass of ice in snow layers is

.. math::
   :label: 8.17

   \frac{\partial w_{ice,\, i} }{\partial t} =
   \left\{\begin{array}{lr}
   f_{sno} \ q_{ice,\, i-1} -\frac{\left(\Delta w_{ice,\, i} \right)_{p} }{\Delta t} & \qquad i=snl+1 \\
   -\frac{\left(\Delta w_{ice,\, i} \right)_{p} }{\Delta t} & \qquad i=snl+2,\ldots ,0
   \end{array}\right\}

where :math:`q_{ice,\, i-1}` is the rate of ice accumulation from precipitation or frost or the rate of ice loss from sublimation (kg m\ :sup:`-2` s\ :sup:`-1`) in the top layer and :math:`{\left(\Delta w_{ice,\, i} \right)_{p} \mathord{\left/ {\vphantom {\left(\Delta w_{ice,\, i} \right)_{p} \Delta t}} \right.} \Delta t}` is the change in ice due to phase change (melting rate) (section :numref:`Phase Change`). The term :math:`q_{ice,\, i-1}` is computed in two steps as

.. math::
   :label: 8.18

   q_{ice,\, i-1} =q_{grnd,\, ice} +\left(q_{frost} -q_{subl} \right)

where :math:`q_{grnd,\, ice}` is the rate of solid precipitation reaching the ground (section :numref:`Canopy Water`) and :math:`q_{frost}` and :math:`q_{subl}` are gains due to frost and losses due to sublimation, respectively (sectio :numref:`Update of Ground Sensible and Latent Heat Fluxes`). In the first step, immediately after :math:`q_{grnd,\, ice}` has been determined after accounting for interception (section :numref:`Canopy Water`), a new snow depth :math:`z_{sno}` (m) is calculated from

.. math::
   :label: 8.19

   z_{sno}^{n+1} =z_{sno}^{n} +\Delta z_{sno}

where

.. math::
   :label: 8.20

   \Delta z_{sno} =\frac{q_{grnd,\, ice} \Delta t}{f_{sno} \rho _{sno} }

and :math:`\rho _{sno}` is the bulk density of newly fallen snow (kg m\ :sup:`-3`), which parameterized by a temperature-dependent and a wind-dependent term:

.. math::
   :label: 8.21a

   \rho_{sno} = \rho_{T} + \rho_{w}.

The temperature dependent term is given by (:ref:`van Kampenhout et al. (2017) <vanKampenhoutetal2017>`)

.. math::
   :label: 8.21b

   \rho_{T} =
   \left\{\begin{array}{lr}
   50 + 1.7 \left(17\right)^{1.5} & \qquad T_{atm} >T_{f} +2 \ \\
   50+1.7 \left(T_{atm} -T_{f} + 15\right)^{1.5} & \qquad T_{f} - 15 < T_{atm} \le T_{f} + 2 \ \\
   -3.833 \ \left( T_{atm} -T_{f} \right) - 0.0333 \ \left( T_{atm} -T_{f} \right)^{2}
   &\qquad T_{atm} \le T_{f} - 15
   \end{array}\right\}

.. bifall(c) = -(50._r8/15._r8 + 0.0333_r8*15_r8)*(forc_t(c)-tfrz) - 0.0333_r8*(forc_t(c)-tfrz)**2

where :math:`T_{atm}` is the atmospheric temperature (K), and :math:`T_{f}` is the freezing temperature of water (K) (:numref:`Table Physical Constants`). When 10 m wind speed :math:`W_{atm}` is greater than 0.1 m\ :sup:`-1`, snow density increases due to wind-driven compaction according to :ref:`van Kampenhout et al. 2017 <vanKampenhoutetal2017>`

.. math::
   :label: 8.21c

   \rho_{w} = 266.861 \left(\frac{1 + \tanh(\frac{W_{atm}}{5})}{2}\right)^{8.8}

.. bifall(c) = bifall(c) + (266.861_r8 * ((1._r8 + TANH(forc_wind(g)/5.0_r8))/2._r8)**8.8_r8)

which is added to the temperature-dependent term (cf. equation :eq:`8.21a`).

The mass of snow :math:`W_{sno}`  is

.. math::
   :label: 8.22

   W_{sno}^{n+1} =W_{sno}^{n} +q_{grnd,\, ice} \Delta t.

The ice content of the top layer and the layer thickness are updated as

.. math::
   :label: 8.23

   w_{ice,\, snl+1}^{n+1} =w_{ice,\, snl+1}^{n} +q_{grnd,\, ice} \Delta t

.. math::
   :label: 8.24

   \Delta z_{snl+1}^{n+1} =\Delta z_{snl+1}^{n} +\Delta z_{sno} .

In the second step, after surface fluxes and snow/soil temperatures have been determined (Chapters :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes` and :numref:`rst_Soil and Snow Temperatures`), :math:`w_{ice,\, snl+1}` is updated for frost or sublimation as

.. math::
   :label: 8.25

   w_{ice,\, snl+1}^{n+1} =w_{ice,\, snl+1}^{n} +f_{sno} \left(q_{frost} -q_{subl} \right)\Delta t.

If :math:`w_{ice,\, snl+1}^{n+1} <0` upon solution of equation :eq:`8.25`, the ice content is reset to zero and the liquid water content :math:`w_{liq,\, snl+1}` is reduced by the amount required to bring :math:`w_{ice,\, snl+1}^{n+1}` up to zero.

The snow water equivalent :math:`W_{sno}` is capped to not exceed 10,000 kg m\ :sup:`-2`. If the addition of :math:`q_{frost}` were to result in :math:`W_{sno} > 10,000` kg m\ :sup:`-2`, the frost term :math:`q_{frost}` is instead added to the ice runoff term :math:`q_{snwcp,\, ice}` (section :numref:`Runoff from glaciers and snow-capped surfaces`).

.. _Water Content:

Water Content
^^^^^^^^^^^^^^^^^^^

The conservation equation for mass of water in snow layers is

.. math::
   :label: 8.26

   \frac{\partial w_{liq,\, i} }{\partial t} =\left(q_{liq,\, i-1} -q_{liq,\, i} \right)+\frac{\left(\Delta w_{liq,\, i} \right)_{p} }{\Delta t}

where :math:`q_{liq,\, i-1}` is the flow of liquid water into layer :math:`i` from the layer above, :math:`q_{liq,\, i}` is the flow of water out of layer :math:`i` to the layer below, :math:`{\left(\Delta w_{liq,\, i} \right)_{p} \mathord{\left/ {\vphantom {\left(\Delta w_{liq,\, i} \right)_{p} \Delta t}} \right.} \Delta t}` is the change in liquid water due to phase change (melting rate) (section :numref:`Phase Change`). For the top snow layer only,

.. math::
   :label: 8.27

   q_{liq,\, i-1} =f_{sno} \left(q_{grnd,\, liq} +\left(q_{sdew} -q_{seva} \right)\right)

where :math:`q_{grnd,\, liq}` is the rate of liquid precipitation reaching the snow (section :numref:`Canopy Water`), :math:`q_{seva}` is the evaporation of liquid water and :math:`q_{sdew}` is the liquid dew (section :numref:`Update of Ground Sensible and Latent Heat Fluxes`). After surface fluxes and snow/soil temperatures have been determined (Chapters :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes` and :numref:`rst_Soil and Snow Temperatures`), :math:`w_{liq,\, snl+1}` is updated for the liquid precipitation reaching the ground and dew or evaporation as

.. math::
   :label: 8.28

   w_{liq,\, snl+1}^{n+1} =w_{liq,\, snl+1}^{n} +f_{sno} \left(q_{grnd,\, liq} +q_{sdew} -q_{seva} \right)\Delta t.

When the liquid water within a snow layer exceeds the layer's holding capacity, the excess water is added to the underlying layer, limited by the effective porosity (:math:`1-\theta _{ice}` ) of the layer. The flow of water is assumed to be zero (:math:`q_{liq,\, i} =0`) if the effective porosity of either of the two layers (:math:`1-\theta _{ice,\, i} {\rm \; and\; }1-\theta _{ice,\, i+1}` ) is less than :math:`\theta _{imp} =0.05`, the water impermeable volumetric water content. Thus, water flow between layers, :math:`q_{liq,\, i}`, for layers :math:`i=snl+1,\ldots,0`, is initially calculated as

.. math::
   :label: 8.29

   q_{liq,\, i} =\frac{\rho _{liq} \left[\theta _{liq,\, i} -S_{r} \left(1-\theta _{ice,\, i} \right)\right]f_{sno} \Delta z_{i} }{\Delta t} \ge 0

where the volumetric liquid water :math:`\theta _{liq,\, i}` and ice :math:`\theta _{ice,\, i}` contents are

.. math::
   :label: 8.30

   \theta _{ice,\, i} =\frac{w_{ice,\, i} }{f_{sno} \Delta z_{i} \rho _{ice} } \le 1

.. math::
   :label: 8.31

   \theta _{liq,\, i} =\frac{w_{liq,\, i} }{f_{sno} \Delta z_{i} \rho _{liq} } \le 1-\theta _{ice,\, i} ,

and :math:`S_{r} =0.033` is the irreducible water saturation (snow holds a certain amount of liquid water due to capillary retention after drainage has ceased (:ref:`Anderson (1976) <Anderson1976>`)). The water holding capacity of the underlying layer limits the flow of water :math:`q_{liq,\, i}` calculated in equation :eq:`8.29`, unless the underlying layer is the surface soil layer, as

.. math::
   :label: 8.32

   q_{liq,\, i} \le \frac{\rho _{liq} \left[1-\theta _{ice,\, i+1} -\theta _{liq,\, i+1} \right]\Delta z_{i+1} }{\Delta t} \qquad i=snl+1,\ldots ,-1.

The liquid water content :math:`w_{liq,\, i}`  is updated as

.. math::
   :label: 8.33

   w_{liq,\, i}^{n+1} =w_{liq,\, i}^{n} +\left(q_{i-1} -q_{i} \right)\Delta t.

Equations :eq:`8.29` - :eq:`8.33` are solved sequentially from top (:math:`i=snl+1`) to bottom (:math:`i=0`) snow layer in each time step. The total flow of liquid water reaching the soil surface is then :math:`q_{liq,\, 0}` which is used in the calculation of surface runoff and infiltration (sections :numref:`Surface Runoff` and :numref:`Infiltration`).

.. _Black and organic carbon and mineral dust within snow:

Black and organic carbon and mineral dust within snow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Particles within snow originate from atmospheric aerosol deposition (:math:`D_{sp}` in Table 2.3 (kg m\ :sup:`-2` s\ :sup:`-1`) and influence snow radiative transfer (sections :numref:`Snow Albedo`, :numref:`Snowpack Optical Properties`, and :numref:`Snow Aging`). Particle masses and mixing ratios are represented with a simple mass-conserving scheme. The model maintains masses of the following eight particle species within each snow layer: hydrophilic black carbon, hydrophobic black carbon, hydrophilic organic carbon, hydrophobic organic carbon, and four species of mineral dust with the following particle sizes: 0.1-1.0, 1.0-2.5, 2.5-5.0, and 5.0-10.0 :math:`\mu m`. Each of these species has unique optical properties (:numref:`Table Single-scatter albedo values used for snowpack impurities and ice`) and meltwater removal efficiencies (:numref:`Table Meltwater scavenging`).

The black carbon and organic carbon deposition rates described in Table 2.3 are combined into four categories as follows

.. math::
   :label: 8.34

   D_{bc,\, hphil} =D_{bc,\, dryhphil} +D_{bc,\, wethphil}

.. math::
   :label: 8.35

   D_{bc,\, hphob} =D_{bc,\, dryhphob}

.. math::
   :label: 8.36

   D_{oc,\, hphil} =D_{oc,\, dryhphil} +D_{oc,\, wethphil}

.. math::
   :label: 8.37

   D_{oc,\, hphob} =D_{oc,\, dryhphob}

Deposited particles are assumed to be instantly mixed (homogeneously) within the surface snow layer and are added after the inter-layer water fluxes are computed (section :numref:`Water Content`) so that some aerosol is in the top layer after deposition and is not immediately washed out before radiative calculations are done. Particle masses are then redistributed each time step based on meltwater drainage through the snow column (section :numref:`Water Content`) and snow layer combination and subdivision (section :numref:`Snow Layer Combination and Subdivision`). The change in mass of each of the particle species :math:`\Delta m_{sp,\, i}` (kg m\ :sup:`-2`) is

.. math::
   :label: 8.38

   \Delta m_{sp,\, i} =\left[k_{sp} \left(q_{liq,\, i-1} c_{sp,\, i-1} -q_{liq,\, i} c_{i} \right)+D_{sp} \right]\Delta t

where :math:`k_{sp}` is the meltwater scavenging efficiency that is unique for each species (:numref:`Table Meltwater scavenging`), :math:`q_{liq,\, i-1}` is the flow of liquid water into layer :math:`i` from the layer above, :math:`q_{liq,\, i}` is the flow of water out of layer :math:`i` into the layer below (kg m\ :sup:`-2` s\ :sup:`-1`) (section :numref:`Water Content`), :math:`c_{sp,\, i-1}` and :math:`c_{sp,\, i}` are the particle mass mixing ratios in layers :math:`i-1` and :math:`i` (kg kg\ :sup:`-1`), :math:`D_{sp}` is the atmospheric deposition rate (zero for all layers except layer :math:`snl+1`), and :math:`\Delta t` is the model time step (s). The particle mass mixing ratio is

.. math::
   :label: 8.39

   c_{i} =\frac{m_{sp,\, i} }{w_{liq,\, i} +w_{ice,\, i} } .

Values of :math:`k_{sp}` are partially derived from experiments published by :ref:`Conway et al. (1996) <Conwayetal1996>`. Particles masses are re-distributed proportionately with snow mass when layers are combined or divided, thus conserving particle mass within the snow column. The mass of particles carried out with meltwater through the bottom snow layer is assumed to be permanently lost from the snowpack, and is not maintained within the model.

.. _Table Meltwater scavenging:

.. table:: Meltwater scavenging efficiency for particles within snow

 +------------------------------------------+-------------------+
 | Species                                  | :math:`k_{sp}`    |
 +==========================================+===================+
 | Hydrophilic black carbon                 | 0.20              |
 +------------------------------------------+-------------------+
 | Hydrophobic black carbon                 | 0.03              |
 +------------------------------------------+-------------------+
 | Hydrophilic organic carbon               | 0.20              |
 +------------------------------------------+-------------------+
 | Hydrophobic organic carbon               | 0.03              |
 +------------------------------------------+-------------------+
 | Dust species 1 (0.1-1.0 :math:`\mu m`)   | 0.02              |
 +------------------------------------------+-------------------+
 | Dust species 2 (1.0-2.5 :math:`\mu m`)   | 0.02              |
 +------------------------------------------+-------------------+
 | Dust species 3 (2.5-5.0 :math:`\mu m`)   | 0.01              |
 +------------------------------------------+-------------------+
 | Dust species 4 (5.0-10.0 :math:`\mu m`)  | 0.01              |
 +------------------------------------------+-------------------+

.. _Initialization of snow layer:

Initialization of snow layer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If there are no existing snow layers (:math:`snl+1=1`) but :math:`z_{sno} \ge 0.01` m after accounting for solid precipitation :math:`q_{sno}`, then a snow layer is initialized (:math:`snl=-1`) as follows

.. math::
   :label: 8.40

   \begin{array}{lcr}
   \Delta z_{0} & = & z_{sno}  \\
   z_{o} & = & -0.5\Delta z_{0}  \\
   z_{h,\, -1} & = & -\Delta z_{0}  \\
   T_{0} & = & \min \left(T_{f} ,T_{atm} \right) \\
   w_{ice,\, 0} & = & W_{sno}  \\
   w_{liq,\, 0} & = & 0
   \end{array}.

.. _Snow Compaction:

Snow Compaction
^^^^^^^^^^^^^^^^^^^^^

Snow compaction is initiated after the soil hydrology calculations [surface runoff (section :numref:`Surface Runoff`), infiltration (section :numref:`Infiltration`), soil water (section :numref:`Soil Water`)] are complete. Currently, there are four processes included that lead to snow compaction:

  #. destructive metamorphism of new snow (crystal breakdown due to wind or thermodynamic stress)
  #. snow load or compaction by overburden pressure
  #. melting (changes in snow structure due to melt-freeze cycles plus changes in crystals due to liquid water)
  #. drifting snow compaction.

The total fractional compaction rate for each snow layer :math:`C_{R,\, i}` (s\ :sup:`-1`) is the sum of multiple compaction processes

.. math::
   :label: 8.41

   C_{R,\, i} =\frac{1}{\Delta z_{i} } \frac{\partial \Delta z_{i} }{\partial t} =C_{R1,\, i} +C_{R2,\, i} +C_{R3,\, i} +C_{R4,\, i} +C_{R5,\, i} .

Compaction is not allowed if the layer is saturated

.. math::
   :label: 8.42

   1-\left(\frac{w_{ice,\, i} }{f_{sno} \Delta z_{i} \rho _{ice} } +\frac{w_{liq,\, i} }{f_{sno} \Delta z_{i} \rho _{liq} } \right)\le 0.001

or if the ice content is below a minimum value (:math:`w_{ice,\, i} \le 0.1`).

The snow layer thickness after compaction is

.. math::
   :label: 8.42b

   \Delta z_{i}^{n+1} =\Delta z_{i}^{n} \left(1+C_{R,\, i} \Delta t\right).

.. _Destructive metamorphism:

Destructive metamorphism
''''''''''''''''''''''''

Compaction as a result of destructive metamorphism :math:`C_{R1,\; i}` (s\ :sup:`-1`) is temperature dependent (:ref:`Anderson (1976) <Anderson1976>`)

.. math::
   :label: 8.43

   C_{R1,\, i} =\left[\frac{1}{\Delta z_{i} } \frac{\partial \Delta z_{i} }{\partial t} \right]_{metamorphism} =-c_{3} c_{1} c_{2} \exp \left[-c_{4} \left(T_{f} -T_{i} \right)\right]

where :math:`c_{3} =2.777\times 10^{-6}`  (s\ :sup:`-1`) is the fractional compaction rate for :math:`T_{i} =T_{f}`, :math:`c_{4} =0.04` K\ :sup:`-1`, and

.. math::
   :label: 8.44

   \begin{array}{lr}
   c_{1}  = 1 & \qquad \frac{w_{ice,\, i} }{f_{sno} \Delta z_{i} } \le 175{\rm \; kg\; m}^{{\rm -3}}  \\
   c_{1}  = \exp \left[-0.046\left(\frac{w_{ice,\, i} }{f_{sno} \Delta z_{i} } -175\right)\right] & \qquad \frac{w_{ice,\, i} }{f_{sno} \Delta z_{i} } >175{\rm \; kg\; m}^{{\rm -3}} \\
   c_{2}  = 2 & \qquad \frac{w_{liq,\, i} }{f_{sno} \Delta z_{i} } >0.01 \\
   c_{2}  = 1 & \qquad \frac{w_{liq,\, i} }{f_{sno} \Delta z_{i} } \le 0.01
   \end{array}

..  upper limit (upplim_destruct_metamorph) used to be 100 but was changed to 175 for CLM5 (Van Kampenhout et al., 2017)

where :math:`{w_{ice,\, i} \mathord{\left/ {\vphantom {w_{ice,\, i}  \left(f_{sno} \Delta z_{i} \right)}} \right.} \left(f_{sno} \Delta z_{i} \right)}` and
:math:`{w_{liq,\, i} \mathord{\left/ {\vphantom {w_{liq,\, i}  \left(f_{sno} \Delta z_{i} \right)}} \right.} \left(f_{sno} \Delta z_{i} \right)}` are the bulk densities of liquid water and ice (kg m\ :sup:`-3`).

.. _Overburden pressure compaction:

Overburden pressure compaction
''''''''''''''''''''''''''''''

The compaction rate as a result of overburden :math:`C_{R2,\; i}` (s\ :sup:`-1`) is a linear function of the snow load pressure :math:`P_{s,\, i}` (kg m\ :sup:`-2`) (:ref:`Anderson (1976) <Anderson1976>`):

.. math::
   :label: 8.45

   C_{R2,\, i} =\left[\frac{1}{\Delta z_{i} } \frac{\partial \Delta z_{i} }{\partial t} \right]_{overburden} =-\frac{P_{s,\, i} }{\eta }

The snow load pressure :math:`P_{s,\, i}` is calculated for each layer as the sum of the ice :math:`w_{ice,\, i}` and liquid water contents :math:`w_{liq,\, i}` of the layers above plus half the ice and liquid water contents of the layer being compacted

.. math::
   :label: 8.47

   P_{s,\, i} =\frac{w_{ice,\, i} +w_{liq,\, i} }{2} +\sum _{j=snl+1}^{j=i-1}\left(w_{ice,\, j} +w_{liq,\, j} \right) .

Variable :math:`\eta` in :eq:`8.45` is a viscosity coefficient (kg s m\ :sup:`-2`) that varies with density and temperature as

.. math::
   :label: 8.46

   \eta = f_{1} f_{2} \eta_{0} \frac{\rho_{i}}{c_{\eta}} \exp \left[ a_{\eta} \left(T_{f} -T_{i} \right) + b_{\eta} \rho_{i} \right]

with constant factors :math:`\eta _{0} = 7.62237 \times 10^{6}` kg s\ :sup:`-1` m\ :sup:`-2`, :math:`a_{\eta} = 0.1` K\ :sup:`-1`, :math:`b_{\eta} = 0.023` m\ :sup:`-3` kg\ :sup:`-1`, and :math:`c_{\eta} = 450` kg m\ :sup:`-3` (:ref:`van Kampenhout et al. (2017) <vanKampenhoutetal2017>`). Further, factor :math:`f_1` accounts for the presence of liquid water (:ref:`Vionnet et al. (2012) <Vionnetetal2012>`):

.. math::
   :label: 8.46b

   f_{1} = \frac{1}{1+ 60 \frac{w_{\mathrm{liq},\, i}}{\rho_{\mathrm{liq}} \Delta z_{i} }}.

Factor :math:`f_2` originally accounts for the presence of angular grains, but since grain shape is not modelled :math:`f_2` is fixed to the value 4.

.. _Compaction by melt:

Compaction by melt
''''''''''''''''''
The compaction rate due to melting :math:`C_{R3,\; i}` (s\ :sup:`-1`) is taken to be the ratio of the change in snow ice mass after the melting to the mass before melting

.. math::
   :label: 8.48

   C_{R3,\, i} = \left[\frac{1}{\Delta z_{i} } \frac{\partial \Delta z_{i} }{\partial t} \right]_{melt}
   = -\frac{1}{\Delta t} \max \left(0,\frac{W_{sno,\, i}^{n} -W_{sno,\, i}^{n+1} }{W_{sno,\, i}^{n} } \right)

and melting is identified during the phase change calculations (section :numref:`Phase Change`). Because snow depth is defined as the average depth of the snow covered area, the snow depth must also be updated for changes in :math:`f_{sno}` when :math:`W_{sno}` has changed.

 .. math::
    :label: 8.49

    C_{R4,\, i} =\left[\frac{1}{\Delta z_{i} } \frac{\partial \Delta z_{i} }{\partial t} \right]_{fsno} =-\frac{1}{\Delta t} \max \left(0,\frac{f_{sno,\, i}^{n} -f_{sno,\, i}^{n+1} }{f_{sno,\, i}^{n} } \right)

.. _Compaction by drifting snow:

Compaction by drifting snow
'''''''''''''''''''''''''''
Crystal breaking by drifting snow leads to higher snow densities at the surface. This process is particularly important on ice sheets, where destructive metamorphism is slow due to low temperatures but high wind speeds (katabatic winds) are prevailing. Therefore a drifting snow compaction parametrization was introduced, based on (:ref:`Vionnet et al. (2012) <Vionnetetal2012>`).

.. math::
   :label: 8.50

   C_{R5,\, i} = \left[\frac{1}{\Delta z_{i} } \frac{\partial \Delta z_{i} }{\partial t} \right]_{drift} = - \frac{\rho_{\max} - \rho_i}{\tau_{i}}.

Here, :math:`\rho_{\max} = 350` kg m\ :sup:`-3` is the upper limit to which this process is active, and
:math:`\tau_{i}` is a timescale which is depth dependent:

.. math::
   :label: 8.50b

   \tau_i = \frac{\tau}{\Gamma_{\mathrm{drift}}^i} \quad \mathrm{,} \:\; \Gamma^i_\mathrm{drift} = \max\left[ 0, S_\mathrm{I}^i \exp(-z_i / 0.1) \right].

Here, :math:`\tau` is a characteristic time scale for drifting snow compaction and is empirically set to 48 h, and
:math:`z_i` is a pseudo-depth which takes into account previous hardening of snow layers above the current layer:
:math:`z_i = \sum_j \Delta z_j \cdot (3.25 - S_\mathrm{I}^j)`.
The driftability index :math:`S_\mathrm{I}` reflects how well snow can be drifted and depends on the mobility of the snow
as well as the 10 m wind speed:

.. math::
   :label: 8.50c

   \begin{array}{rcl}
   S_\mathrm{I} & = & -2.868 \exp(-0.085 U) + 1 + M_{\mathrm{O}} \\
   M_\mathrm{O} & = & -0.069 + 0.66 F(\rho)
   \end{array}

The latter equation (for the mobility index :math:`M_\mathrm{O}`) is a simplification from the original paper by removing the dependency on grain size and assuming spherical grains (see :ref:`van Kampenhout et al. (2017) <vanKampenhoutetal2017>`).

.. _Snow Layer Combination and Subdivision:

Snow Layer Combination and Subdivision
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the determination of snow temperature including phase change(Chapter :numref:`rst_Soil and Snow Temperatures`), snow hydrology (Chapter :numref:`rst_Snow Hydrology`), and the compaction calculations (section :numref:`Snow Compaction`), the number of snow layers is adjusted by either combining or subdividing layers. The combination and subdivision of snow layers is based on :ref:`Jordan (1991) <Jordan1991>`.

.. _Combination:

Combination
'''''''''''''''''''

If a snow layer has nearly melted or if its thickness :math:`\Delta z_{i}` is less than the prescribed minimum thickness :math:`\Delta z_{\min }` (:numref:`Table snow layer thickness`), the layer is combined with a neighboring layer. The overlying or underlying layer is selected as the neighboring layer according to the following rules

#. If the top layer is being removed, it is combined with the underlying layer

#. If the underlying layer is not snow (i.e., it is the top soil layer), the layer is combined with the overlying layer

#. If the layer is nearly completely melted, the layer is combined with the underlying layer

#. If none of the above rules apply, the layer is combined with the thinnest neighboring layer.

A first pass is made through all snow layers to determine if any layer is nearly melted (:math:`w_{ice,\, i} \le 0.1`). If so, the remaining liquid water and ice content of layer :math:`i` is combined with the underlying neighbor :math:`i+1` as

.. math::
   :label: 8.51

   w_{liq,\, i+1} =w_{liq,\, i+1} +w_{liq,\, i}

.. math::
   :label: 8.52

   w_{ice,\, i+1} =w_{ice,\, i+1} +w_{ice,\, i} .

This includes the snow layer directly above the top soil layer. In this case, the liquid water and ice content of the melted snow layer is added to the contents of the top soil layer. The layer properties, :math:`T_{i}`, :math:`w_{ice,\, i}`, :math:`w_{liq,\, i}`, :math:`\Delta z_{i}`, are then re-indexed so that the layers above the eliminated layer are shifted down by one and the number of snow layers is decremented accordingly.

At this point, if there are no explicit snow layers remaining (:math:`snl=0`), the snow water equivalent :math:`W_{sno}` and snow depth :math:`z_{sno}` are set to zero, otherwise, :math:`W_{sno}` and :math:`z_{sno}` are re-calculated as

.. math::
   :label: 8.53

   W_{sno} =\sum _{i=snl+1}^{i=0}\left(w_{ice,\, i} +w_{liq,\, i} \right)

.. math::
   :label: 8.54

   z_{sno} =\sum _{i=snl+1}^{i=0}\Delta z_{i}  .

If the snow depth :math:`0<z_{sno} <0.01` m or the snow density :math:`\frac{W_{sno} }{f_{sno} z_{sno} } <50` kg/m3, the number of snow layers is set to zero, the total ice content of the snowpack :math:`\sum _{i=snl+1}^{i=0}w_{ice,\; i}` is assigned to :math:`W_{sno}`, and the total liquid water :math:`\sum _{i=snl+1}^{i=0}w_{liq,\; i}` is assigned to the top soil layer. Otherwise, the layers are combined according to the rules above.

When two snow layers are combined (denoted here as 1 and 2), their thickness combination (:math:`c`) is

.. math::
   :label: 8.55

   \Delta z_{c} =\Delta z_{1} +\Delta z_{2} ,

their mass combination is

.. math::
   :label: 8.56

   w_{liq,\, c} =w_{liq,\, 1} +w_{liq,\, 2}

.. math::
   :label: 8.57

   w_{ice,\, c} =w_{ice,\, 1} +w_{ice,\, 2} ,

and their temperatures are combined as

.. math::
   :label: 8.58

   T_{c} =T_{f} +\frac{h_{c} -L_{f} w_{liq,\, c} }{C_{ice} w_{ice,\, c} +C_{liq} w_{liq,\, c} }

where :math:`h_{c} =h_{1} +h_{2}` is the combined enthalpy :math:`h_{i}` of the two layers where

.. math::
   :label: 8.59

   h_{i} =\left(C_{ice} w_{ice,\, i} +C_{liq} w_{liq,\, i} \right)\left(T_{i} -T_{f} \right)+L_{f} w_{liq,\, i} .

In these equations, :math:`L_{f}` is the latent heat of fusion (J kg\ :sup:`-1`) and :math:`C_{liq}` and :math:`C_{ice}` are the specific heat capacities (J kg\ :sup:`-1` K\ :sup:`-1`) of liquid water and ice, respectively (:numref:`Table Physical Constants`). After layer combination, the node depths and layer interfaces (:numref:`Figure three layer snow pack`) are recalculated from

.. math::
   :label: 8.60

   z_{i} =z_{h,\, i} -0.5\Delta z_{i} \qquad i=0,\ldots ,snl+1

.. math::
   :label: 8.61

   z_{h,\, i-1} =z_{h,\, i} -\Delta z_{i} \qquad i=0,\ldots ,snl+1

where :math:`\Delta z_{i}`  is the layer thickness.

.. _Table snow layer thickness:

.. table:: Minimum and maximum thickness of snow layers (m)

 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | Layer        | :math:`\Delta z_{\min }`    | :math:`N_{l}`    | :math:`N_{u}`    | :math:`\left(\Delta z_{\max } \right)_{l}`    | :math:`\left(\Delta z_{\max } \right)_{u}`              |
 +==============+=============================+==================+==================+===============================================+=========================================================+
 | 1 (top)      | 0.010                       | 1                | :math:`>`\ 1     | 0.03                                          | 0.02                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 2            | 0.015                       | 2                | :math:`>`\ 2     | 0.07                                          | 0.05                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 3            | 0.025                       | 3                | :math:`>`\ 3     | 0.18                                          | 0.11                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 4            | 0.055                       | 4                | :math:`>`\ 4     | 0.41                                          | 0.23                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 5            | 0.115                       | 5                | :math:`>`\ 5     | 0.88                                          | 0.47                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 6            | 0.235                       | 6                | :math:`>`\ 6     | 1.83                                          | 0.95                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 7            | 0.475                       | 7                | :math:`>`\ 7     | 3.74                                          | 1.91                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 8            | 0.955                       | 8                | :math:`>`\ 8     | 7.57                                          | 3.83                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 9            | 1.915                       | 9                | :math:`>`\ 9     | 15.24                                         | 7.67                                                    |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 10           | 3.835                       | 10               | :math:`>`\ 10    | 30.59                                         | 15.35                                                   |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 11           | 7.675                       | 11               | :math:`>`\ 11    | 61.30                                         | 30.71                                                   |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+
 | 12 (bottom)  | 15.355                      | 12               | -                | -                                             | -                                                       |
 +--------------+-----------------------------+------------------+------------------+-----------------------------------------------+---------------------------------------------------------+

The maximum snow layer thickness, :math:`\Delta z_{\max }`, depends on the number of layers, :math:`N_{l}` and :math:`N_{u}` (section :numref:`Subdivision`).

.. _Subdivision:

Subdivision
'''''''''''''''''''

The snow layers are subdivided when the layer thickness exceeds the prescribed maximum thickness :math:`\Delta z_{\max }` with lower and upper bounds that depend on the number of snow layers (:numref:`Table snow layer thickness`). For example, if there is only one layer, then the maximum thickness of that layer is 0.03 m, however, if there is more than one layer, then the maximum thickness of the top layer is 0.02 m. Layers are checked sequentially from top to bottom for this limit. If there is only one snow layer and its thickness is greater than 0.03 m (:numref:`Table snow layer thickness`), the layer is subdivided into two layers of equal thickness, liquid water and ice contents, and temperature. If there is an existing layer below the layer to be subdivided, the thickness :math:`\Delta z_{i}`, liquid water and ice contents, :math:`w_{liq,\; i}` and :math:`w_{ice,\; i}`, and temperature :math:`T_{i}` of the excess snow are combined with the underlying layer according to equations -. If there is no underlying layer after adjusting the layer for the excess snow, the layer is subdivided into two layers of equal thickness, liquid water and ice contents. The vertical snow temperature profile is maintained by calculating the slope between the layer above the splitting layer (:math:`T_{1}` ) and the splitting layer (:math:`T_{2}` ) and constraining the new temperatures (:math:`T_{2}^{n+1}`, :math:`T_{3}^{n+1}` ) to lie along this slope. The temperature of the lower layer is first evaluated from

.. math::
   :label: 8.62

   T'_{3} =T_{2}^{n} -\left(\frac{T_{1}^{n} -T_{2}^{n} }{{\left(\Delta z_{1}^{n} +\Delta z_{2}^{n} \right)\mathord{\left/ {\vphantom {\left(\Delta z_{1}^{n} +\Delta z_{2}^{n} \right) 2}} \right.} 2} } \right)\left(\frac{\Delta z_{2}^{n+1} }{2} \right),

then adjusted as,

.. math::
   :label: 8.63

   \begin{array}{lr}
   T_{3}^{n+1} = T_{2}^{n} & \qquad T'_{3} \ge T_{f}  \\
   T_{2}^{n+1} = T_{2}^{n} +\left(\frac{T_{1}^{n} -T_{2}^{n} }{{\left(\Delta z_{1} +\Delta z_{2}^{n} \right)\mathord{\left/ {\vphantom {\left(\Delta z_{1} +\Delta z_{2}^{n} \right) 2}} \right.} 2} } \right)\left(\frac{\Delta z_{2}^{n+1} }{2} \right) & \qquad T'_{3} <T_{f}
   \end{array}

where here the subscripts 1, 2, and 3 denote three layers numbered from top to bottom. After layer subdivision, the node depths and layer interfaces are recalculated from equations and.

