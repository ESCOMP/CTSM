.. _rst_Fire:

Fire
========

The fire parameterization in CLM contains four components: non-peat fires outside cropland and tropical closed forests, agricultural fires in cropland, deforestation fires in the tropical closed forests, and peat fires (see :ref:`Li et al. 2012a <Lietal2012a>`, :ref:`Li et al. 2012b <Lietal2012b>`, :ref:`Li et al. 2013 <Lietal2013a>`, :ref:`Li and Lawrence 2017 <LiLawrence2017>` for details). In this fire parameterization, burned area is affected by climate and weather conditions, vegetation composition and structure, and human activities. After burned area is calculated, we estimate the fire impact, including biomass and peat burning, fire-induced vegetation mortality, adjustment of the carbon and nitrogen (C/N) pools, and fire emissions.

.. _Non-peat fires outside cropland and tropical closed forest:

Non-peat fires outside cropland and tropical closed forest
---------------------------------------------------------------

Burned area in a grid cell, \ :math:`A_{b}` (km\ :sup:`2` s :sup:`-1`), is determined by

.. math::
   :label: 23.1

   A_{b} =N_{f} a

where :math:`N_{f}` (count s\ :sup:`-1`) is fire counts in the grid cell; :math:`a` (km\ :sup:`2`) is average fire spread area of a fire.

.. _Fire counts:

Fire counts
^^^^^^^^^^^^^^^^^^

Fire counts :math:`N_{f}`  is taken as

.. math::
   :label: 23.2

   N_{f} = N_{i} f_{b} f_{m} f_{se,o}

where :math:`N_{i}` ( count s\ :sup:`-1`) is the number of ignition sources due to natural causes and human activities; :math:`f_{b}` and :math:`f_{m}` (fractions) represent the availability and combustibility of fuel, respectively; :math:`f_{se,o}` is the fraction of anthropogenic and natural fires unsuppressed by humans and related to the socioeconomic conditions.

:math:`N_{i}`  (count s\ :sup:`-1`) is given as

.. math::
   :label: 23.3

   N_{i} = \left(I_{n} +I_{a} \right) A_{g}

where :math:`I_{n}` (count km\ :sup:`-2` s\ :sup:`-1`) and :math:`I_{a}` (count km\ :sup:`-2` s\ :sup:`-1`) are the number of natural and anthropogenic ignitions per km\ :sup:`2`, respectively; :math:`A_{g}` is the area of the grid cell (km\ :sup:`2`). :math:`I_{n}` is estimated by

.. math::
   :label: 23.4

   I_{n} = \gamma \psi I_{l}

where :math:`\gamma` \ =0.22 is ignition efficiency of cloud-to-ground lightning; :math:`\psi =\frac{1}{5.16+2.16\cos [3min(60,\lambda )]}` is the cloud-to-ground lightning fraction and depends on the latitude :math:`\lambda` (degrees); :math:`I_{l}` (flash km\ :sup:`-2` s\ :sup:`-1`) is the total lightning flashes. :math:`I_{a}` is modeled as a monotonic increasing function of population density:

.. math::
   :label: 23.5

   I_{a} =\frac{\alpha D_{P} k(D_{P} )}{n}

where :math:`\alpha =0.01` (count person\ :sup:`-1` mon\ :sup:`-1`) is the number of potential ignition sources by a person per month; :math:`D_{P}` (person km\ :sup:`-2`) is the population density; :math:`k(D_{P} )=6.8D_{P} ^{-0.6}` represents anthropogenic ignition potential as a function of human population density :math:`D_{P}`; *n* is the seconds in a month.

Fuel availability :math:`f_{b}` is given as

.. math::
   :label: 23.6

   f_{b} =\left\{\begin{array}{c}
   {0} \\ {\frac{B_{ag} -B_{low} }{B_{up} -B_{low} } } \\ {1} \end{array}
   \begin{array}{cc} {} & {} \end{array}\begin{array}{c} {B_{ag} <B_{low} } \\ {\begin{array}{cc} {} & {} \end{array}B_{low} \le B_{ag} \le B_{up} } \\ {B_{ag} >B_{up} }
   \end{array}\right\} \ ,

where :math:`B_{ag}` (g C m\ :sup:`-2`) is the biomass of combined leaf, stem, litter, and woody debris pools; :math:`B_{low}` = 105 g C m :sup:`-2` is the lower fuel threshold below which fire does not occur; :math:`B_{up}` = 1050 g C m\ :sup:`-2` is the upper fuel threshold above which fire occurrence is not limited by fuel availability.

Fuel combustibility :math:`f_{m}` is estimated by

.. math::
   :label: 23.7

   f_{m} = {f_{RH} f_{\beta}}, \qquad T_{17cm} > T_{f}

where :math:`f_{RH}` and :math:`f_{\beta }` represent the dependence of fuel combustibility on relative humidity :math:`RH` (%) and root-zone soil moisture limitation :math:`\beta` (fraction); :math:`T_{17cm}` is the temperature of the top 17 cm of soil (K) and :math:`T_{f}` is the freezing temperature. :math:`f_{RH}` is a weighted average of real time :math:`RH` (:math:`RH_{0}`) and 30-day running mean :math:`RH` (:math:`RH_{30d}`):

.. math::
   :label: 23.8

   f_{RH} = (1-w) l_{RH_{0}} + wl_{RH_{30d}}

where weight :math:`w=\max [0,\min (1,\frac{B_{ag}-2500}{2500})]`, :math:`l_{{RH}_{0}}=1-\max [0,\min (1,\frac{RH_{0}-30}{80-30})]`, and :math:`l_{{RH}_{30d}}=1-\max [0.75,\min (1,\frac{RH_{30d}}{90})]`. :math:`f_{\beta}` is given by

.. math::
   :label: 23.9

   f_{\beta } =\left\{\begin{array}{cccc}
   {1} & {} & {} & {\beta\le \beta_{low} } \\ {\frac{\beta_{up} -\beta}{\beta_{up} -\beta_{low} } } & {} & {} & {\beta_{low} <\beta<\beta_{up} } \\
   {0} & {} & {} & {\beta\ge \beta_{up} }
   \end{array}\right\} \ ,

where :math:`\beta _{low}` \ =0.85 and :math:`\beta _{up}` \ =0.98 are the lower and upper thresholds, respectively.

For scarcely populated regions (:math:`D_{p} \le 0.1` person km :sup:`-2`), we assume that anthropogenic suppression on fire occurrence is negligible, i.e., :math:`f_{se,o} =1.0`. In regions of :math:`D_{p} >0.1` person km\ :sup:`-2`, we parameterize the fraction of anthropogenic and natural fires unsuppressed by human activities as

.. math::
   :label: 23.10

   f_{se,o} =f_{d} f_{e}

where :math:`{f}_{d}` and :math:`{f}_{e}` are the effects of the demographic and economic conditions on fire occurrence. The demographic influence on fire occurrence is

.. math::
   :label: 23.11

   f_{d} =0.01 + 0.98 \exp (-0.025D_{P} ).

For shrub and grass PFTs, the economic influence on fire occurrence is parameterized as a function of Gross Domestic Product GDP (k 1995US$ capita\ :sup:`-1`):

.. math::
   :label: 23.12

   f_{e} =0.1+0.9\times \exp [-\pi (\frac{GDP}{8} )^{0.5} ]

which captures 73% of the observed MODIS fire counts with variable GDP in regions where shrub and grass PFTs are dominant (fractional coverage of shrub and grass PFTs :math:`>` 50%). In regions outside tropical closed forests and dominated by trees (fractional coverage of tree PFTs :math:`>` 50%), we use

.. math::
   :label: 23.13

   f_{e} =\left\{\begin{array}{c}
   {0.39} \\ {0.79} \\ {1} \end{array}
   \begin{array}{cc} {} & {} \end{array}\begin{array}{c} {GDP > 20 } \\
   { 8 < GDP \le 20 } \\  { GDP \le 8 }
   \end{array}\right\} \ ,

to reproduce the relationship between MODIS fire counts and GDP.

.. _Average spread area of a fire:

Average spread area of a fire
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fire fighting capacity depends on socioeconomic conditions and affects fire spread area. Due to a lack of observations, we consider the socioeconomic impact on the average burned area rather than separately on fire spread rate and fire duration:

.. math::
   :label: 23.14

   a=a^{*} F_{se}

where :math:`a^{*}` is the average burned area of a fire without anthropogenic suppression and :math:`F_{se}` is the socioeconomic effect on fire spread area.

Average burned area of a fire without anthropogenic suppression is assumed elliptical in shape with the wind direction along the major axis and the point of ignition at one of the foci. According to the area formula for an ellipse, average burned area of a fire can be represented as:

.. math::
   :label: 23.15

   a^{*} =\pi \frac{l}{2} \frac{w}{2} \times 10^{-6} =\frac{\pi u_{p}^{2} \tau ^{2} }{4L_{B} } (1+\frac{1}{H_{B} } )^{2} \times 10^{-6}

where :math:`u_{p}` (m s\ :sup:`-1`) is the fire spread rate in the downwind direction; :math:`\tau` (s) is average fire duration; :math:`L_{B}` and :math:`H_{B}` are length-to-breadth ratio and head-to-back ratio of the ellipse; 10 :sup:`-6` converts m :sup:`2` to km :sup:`2`.

According to :ref:`Arora and Boer (2005)<AroraBoer2005>`,

.. math::
   :label: 23.16

   L_{B} =1.0+10.0[1-\exp (-0.06W)]

where :math:`W`\ (m s\ :sup:`-1`) is the wind speed. According to the mathematical properties of the ellipse, the head-to-back ratio :math:`H_{B}` is

.. math::
   :label: 23.17

   H_{B} =\frac{u_{p} }{u_{b} } =\frac{L_{B} +(L_{B} ^{2} -1)^{0.5} }{L_{B} -(L_{B} ^{2} -1)^{0.5} } .

The fire spread rate in the downwind direction is represented as

.. math::
   :label: 23.18

   u_{p} =u_{\max } C_{m} g(W)

(:ref:`Arora and Boer, 2005<AroraBoer2005>`), where :math:`u_{\max }` (m s\ :sup:`-1`) is the PFT-dependent average maximum fire spread rate in natural vegetation regions; :math:`C_{m} =\sqrt{f_{m}}` and :math:`g(W)` represent the dependence of :math:`u_{p}` on fuel wetness and wind speed :math:`W`, respectively. :math:`u_{\max }` is set to 0.33 m s :sup:`-1`\ for grass PFTs, 0.28 m s :sup:`-1` for shrub PFTs, 0.26 m s\ :sup:`-1` for needleleaf tree PFTs, and 0.25 m s\ :sup:`-1` for other tree PFTs. :math:`g(W)` is derived from the mathematical properties of the ellipse and equation :eq:`23.16` and :eq:`23.17`.

.. math::
   :label: 23.19

   g(W)=\frac{2L_{B} }{1+\frac{1}{H_{B} } } g(0).

Since g(\ *W*)=1.0, and \ :math:`L_{B}` and :math:`H_{B}` are at their maxima \ :math:`L_{B} ^{\max } =11.0` and \ :math:`H_{B} ^{\max } =482.0` when :math:`W\to \infty`, g(0) can be derived as

.. math::
   :label: 23.20

   g(0)=\frac{1+\frac{1}{H_{B} ^{\max } } }{2L_{B} ^{\max } } =0.05.

In the absence of globally gridded data on barriers to fire (e.g. rivers, lakes, roads, firebreaks) and human fire-fighting efforts, average fire duration is simply assumed equal to 1 which is the observed 2001–2004 mean persistence of most fires in the world (:ref:`Giglio et al. 2006 <Giglioetal2006>`).

As with the socioeconomic influence on fire occurrence, we assume that the socioeconomic influence on fire spreading is negligible in regions of :math:`D_{p} \le 0.1` person km\ :sup:`-2`, i.e., :math:`F_{se} = 1.0`. In regions of :math:`D_{p} >0.1` person km\ :sup:`-2`, we parameterize such socioeconomic influence as:

.. math::
   :label: 23.21

   F_{se} =F_{d} F_{e}

where :math:`{F}_{d}` and :math:`{F}_{e}` are effects of the demographic and economic conditions on the average spread area of a fire, and are identified by maximizing the explained variability of the GFED3 burned area fraction with both socioeconomic indices in grid cells with various dominant vegetation types. For shrub and grass PFTs, the demographic impact factor is

.. math::
   :label: 23.22

   F_{d} =0.2+0.8\times \exp [-\pi (\frac{D_{p} }{450} )^{0.5} ]

and the economic impact factor is

.. math::
   :label: 23.23

   F_{e} =0.2+0.8\times \exp (-\pi \frac{GDP}{7} ).

For tree PFTs outside tropical closed forests, the demographic and economic impact factors are given as

.. math::
   :label: 23.24

   F_{d} =0.4+0.6\times \exp (-\pi \frac{D_{p} }{125} )

and

.. math::
   :label: 23.25

   F_{e} =\left\{\begin{array}{cc}
   {0.62,} & {GDP>20} \\ {0.83,} & {8<GDP\le 20} \\
   {1,} & {GDP\le 8}
   \end{array}\right. .

Equations :eq:`23.22` - :eq:`23.25` reflect that more developed and more densely populated regions have a higher fire fighting capability.

.. _Fire impact:

Fire impact
^^^^^^^^^^^^^^^^^^

In post-fire regions, we calculate PFT-level fire carbon emissions from biomass burning of the :math:`j`\ th PFT, :math:`{\phi}_{j}` (g C s\ :sup:`-1`), as

.. math::
   :label: 23.26

   \phi _{j} =A_{b,j} \mathbf{C}_{j} \bullet \mathbf{CC}_{j}

where :math:`A_{b,j}` (km\ :sup:`2` \s\ :sup:`-1`) is burned area for the :math:`j`\ th PFT; **C**\ :sub:`j` =(:math:`C_{leaf}`, :math:`C_{stem}`, :math:`C_{root}`, :math:`C_{ts}`) is a vector with carbon density (g C km :sup:`-2`) for leaf, stem (live and dead stem), root (fine, live coarse and dead coarse root), and transfer and storage carbon pools as elements; :math:`\mathbf{CC}_{j}` = (:math:`\mathbf{CC}_{leaf}`, :math:`\mathbf{CC}_{stem}`, :math:`\mathbf{CC}_{root}`, :math:`\mathbf{CC}_{ts}`) is the corresponding combustion completeness factor vector (:numref:`Table PFT-specific combustion completeness and fire mortality factors`). Moreover, we assume that 50% and 28% of column-level litter and coarse woody debris are burned and the corresponding carbon is transferred to atmosphere.

Tissue mortality due to fire leads to carbon transfers in two ways. First, carbon from uncombusted leaf, live stem, dead stem, root, and transfer and storage pools :math:`\mathbf{C^{'} _{j1}} ={(C_{{leaf}} (1-CC_{{leaf}} ),C_{{livestem}} (1-CC_{{stem}} ),C_{{deadstem}} (1-CC_{{stem}} ),C_{{root}} (1-CC_{{root}} ),C_{{ts}} (1-CC_{{ts}} ))}_{j}` (g C km\ :sup:`-2`) is transferred to litter as

.. math::
   :label: 23.27

   \Psi _{j1} =\frac{A_{b,j} }{f_{j} A_{g} } \mathbf{C^{'} _{j1}} \bullet M_{j1}

where :math:`M_{j1} =(M_{{leaf}},M_{{livestem,1}},M_{{deadstem}},M_{{root}},M_{{ts}} )_{j}` is the corresponding mortality factor vector (:numref:`Table PFT-specific combustion completeness and fire mortality factors`). Second, carbon from uncombusted live stems is transferred to dead stems as:

.. math::
   :label: 23.28

   \Psi _{j2} =\frac{A_{b,j} }{f_{j} A_{g} } C_{livestem} (1-CC_{stem} )M_{livestem,2}

where :math:`M_{livestem,2}` is the corresponding mortality factor (:numref:`Table PFT-specific combustion completeness and fire mortality factors`).

Fire nitrogen emissions and nitrogen transfers due to fire-induced mortality are calculated the same way as for carbon, using the same values for combustion completeness and mortality factors. With CLM's dynamic vegetation option enabled, the number of tree PFT individuals killed by fire per km\ :sup:`2` (individual km\ :sup:`-2` s\ :sup:`-1`) is given by

.. math::
   :label: 23.29

   P_{disturb,j} =\frac{A_{b,j} }{f_{j} A_{g} } P_{j} \xi _{j}

where :math:`P_{j}` (individual km\ :sup:`-2`) is the population density for the :math:`j` th tree PFT and :math:`\xi _{j}` is the whole-plant mortality factor (:numref:`Table PFT-specific combustion completeness and fire mortality factors`).

.. _Agricultural fires:

Agricultural fires
-----------------------

The burned area of cropland (km\ :sup:`2` s\ :sup:`-1`) is taken as :math:`{A}_{b}`:

.. math::
   :label: 23.30

   A_{b} =a_{1} f_{se} f_{t} f_{crop} A_{g}

where :math:`a_{1}` (s\ :sup:`-1`) is a constant; :math:`f_{se}` represents the socioeconomic effect on fires; :math:`f_{t}` determines the seasonality of agricultural fires; :math:`f_{crop}` is the fractional coverage of cropland. :math:`a_{1}` \ = 1.6x10\ :sup:`-4` \hr\ :sup:`-1`\ is estimated using an inverse method, by matching 1997-2004 simulations to the analysis of :ref:`van der Werf et al. (2010) <vanderWerfetal2010>` that shows the 2001-2009 average contribution of cropland fires is 4.7% of the total global burned area.

The socioeconomic factor :math:`f_{se}`  is given as follows:

.. math::
   :label: 23.31

   f_{se} =f_{d} f_{e} .

Here

.. math::
   :label: 23.32

   f_{d} =0.04+0.96\times \exp [-\pi (\frac{D_{p} }{350} )^{0.5} ]

and

.. math::
   :label: 23.33

   f_{e} =0.01+0.99\times \exp (-\pi \frac{GDP}{10} )

are the effects of population density and GDP on burned area, derived in a similar way to equation :eq:`23.32` and :eq:`23.33`. :math:`f_{t}` is set to 1 at the first time step during the climatological peak month for agricultural fires (:ref:`van der Werf et al. 2010 <vanderWerfetal2010>`); :math:`{f}_{t}` is set to 0 otherwise. Peak month in this dataset correlates with the month after harvesting or the month before planting. In CLM we use this dataset the same way whether the CROP option is active or not, without regard to the CROP option's simulated planting and harvesting dates.

In the post-fire region, fire impact is parameterized similar to section :numref:`Fire impact` but with combustion completeness factors and tissue mortality factors for crop PFTs (:numref:`Table PFT-specific combustion completeness and fire mortality factors`).

.. _Deforestation fires:

Deforestation fires
------------------------

CLM focuses on deforestation fires in tropical closed forests. Tropical closed forests are defined as grid cells with tropical tree (BET and BDT tropical) coverage :math:`>` 60% according to the FAO classification. Deforestation fires are defined as fires caused by deforestation, including escaped deforestation fires, termed degradation fires. Deforestation and degradation fires are assumed to occur outside of cropland areas in these grid cells. Burned area is controlled by the deforestation rate and climate:

.. math::
   :label: 23.34

   A_{b} = b \ f_{lu} f_{cli,d} f_{b} A_{g}

where :math:`b` (s\ :sup:`-1`) is a global constant; :math:`f_{lu}` (fraction) represents the effect of decreasing fractional coverage of tree PFTs derived from land use data; :math:`f_{cli,d}` (fraction) represents the effect of climate conditions on the burned area.

Constants :math:`b` and :math:`{f}_{lu}` are calibrated based on observations and reanalysis datasets in the Amazon rainforest (tropical closed forests within 15.5 °S :math:`\text{-}` 10.5 °N, 30.5 ° W :math:`\text{-}` 91 ° W). :math:`b` = 0.033 d\ :sup:`-1` and :math:`f_{lu}` is defined as

.. math::
   :label: 23.35

   f_{lu} = \max (0.0005,0.19D-0.001)

where :math:`D` (yr\ :sup:`-1`) is the annual loss of tree cover based on CLM land use and land cover change data.

The effect of climate on deforestation fires is parameterized as:

.. math::
   :label: 23.36

   \begin{array}{ll}
   f_{cli,d} \quad = & \quad \max \left[0,\min (1,\frac{b_{2} -P_{60d} }{b_{2} } )\right]^{0.5} \times \\
   & \quad \max \left[0,\min (1,\frac{b_{3} -P_{10d} }{b_{3} } )\right]^{0.5} \times \\
   & \quad \max \left[0,\min (1,\frac{0.25-P}{0.25} )\right]
   \end{array}

where :math:`P` (mm d :sup:`-1`) is instantaneous precipitation, while :math:`P_{60d}` (mm d\ :sup:`-1`) and :math:`P_{10d}` (mm d :sup:`-1`) are 60-day and 10-day running means of precipitation, respectively; :math:`b_{2}` (mm d :sup:`-1`) and :math:`b_{3}` (mm d :sup:`-1`) are the grid-cell dependent thresholds of :math:`P_{60d}` and :math:`P_{10d}`; 0.25 mm d :sup:`-1` is the maximum precipitation rate for drizzle. :ref:`Le Page et al. (2010) <LePageetal2010>` analyzed the relationship between large-scale deforestation fire counts and precipitation during 2003 :math:`\text{-}`\ 2006 in southern Amazonia where tropical evergreen trees (BET Tropical) are dominant. Figure 2 in :ref:`Le Page et al. (2010) <LePageetal2010>` showed that fires generally occurred if both :math:`P_{60d}` and :math:`P_{10d}` were less than about 4.0 mm d :sup:`-1`, and fires occurred more frequently in a drier environment. Based on the 30-yr (1985 to 2004) precipitation data in :ref:`Qian et al. (2006) <Qianetal2006>`. The climatological precipitation of dry months (P < 4.0 mm d :sup:`-1`) in a year over tropical deciduous tree (BDT Tropical) dominated regions is 46% of that over BET Tropical dominated regions, so we set the PFT-dependent thresholds of :math:`P_{60d}` and :math:`P_{10d}` as 4.0 mm d :sup:`-1` for BET Tropical and 1.8 mm d :sup:`-1` (= 4.0 mm d :sup:`-1` :math:`\times` 46%) for BDT Tropical, and :math:`b`\ :sub:`2` and :math:`b`\ :sub:`3` are the average of thresholds of BET Tropical and BDT Tropical weighted bytheir coverage.

The post-fire area due to deforestation is not limited to land-type conversion regions. In the tree-reduced region, the maximum fire carbon emissions are assumed to be 80% of the total conversion flux. According to the fraction of conversion flux for tropical trees in the tree-reduced region (60%) assigned by CLM4-CN, to reach the maximum fire carbon emissions in a conversion region requires burning this region about twice when we set PFT-dependent combustion completeness factors to about 0.3 for stem [the mean of 0.2\ :math:`{-}`\ 0.4 used in :ref:`van der Werf et al. (2010) <vanderWerfetal2010>`. Therefore, when the burned area calculated from equation :eq:`23.36` is no more than twice the tree-reduced area, we assume no escaped fires outside the land-type conversion region, and the fire-related fraction of the total conversion flux is estimated as :math:`\frac{A_{b} /A_{g} }{2D}`. Otherwise, 80% of the total conversion flux is assumed to be fire carbon emissions, and the biomass combustion and vegetation mortality outside the tree-reduced regions with an area fraction of :math:`\frac{A_{b} }{A_{g} } -2D` are set as in section :numref:`Fire impact`.

.. _Peat fires:

Peat fires
---------------

The burned area due to peat fires is given as :math:`{A}_{b}`:

.. math::
   :label: 23.37

   A_{b} = c \ f_{cli,p} f_{peat} (1 - f_{sat} ) A_{g}

where :math:`c` (s\ :sup:`-1`) is a constant; :math:`f_{cli,p}` represents the effect of climate on the burned area; :math:`f_{peat}` is the fractional coverage of peatland in the grid cell; and :math:`f_{sat}` is the fraction of the grid cell with a water table at the surface or higher. :math:`c` = 0.17 :math:`\times` 10 :sup:`-3` hr\ :sup:`-1` for tropical peat fires and :math:`c` = 0.9 :math:`\times` 10 :sup:`-5` hr :sup:`-1` for boreal peat fires are derived using an inverse method, by matching simulations to earlier studies: about 2.4 Mha peatland was burned over Indonesia in 1997 (:ref:`Page et al. 2002 <Pageetal2002>`) and the average burned area of peat fires in Western Canada was 0.2 Mha yr :sup:`-1` for 1980-1999 (:ref:`Turetsky et al. 2004 <Turetskyetal2004>`).

For tropical peat fires, :math:`f_{cli,p}` is set as a function of long-term precipitation :math:`P_{60d}` :

.. math::
   :label: 23.38

   f_{cli,p} = \ max \left[0,\min \left(1,\frac{4-P_{60d} }{4} \right)\right]^{2} .

For boreal peat fires, :math:`f_{cli,p}`  is set to

.. math::
   :label: 23.39

   f_{cli,p} = \exp (-\pi \frac{\theta _{17cm} }{0.3} )\cdot \max [0,\min (1,\frac{T_{17cm} -T_{f} }{10} )]

where :math:`\theta _{17cm}` is the wetness of the top 17 cm of soil.

Peat fires lead to peat burning and the combustion and mortality of vegetation over peatlands. For tropical peat fires, based on :ref:`Page et al. (2002) <Pageetal2002>`, about 6% of the peat carbon loss from stored carbon is caused by 33.9% of the peatland burned. Carbon emissions due to peat burning (g C m\ :sup:`-2` s\ :sup:`-1`) are therefore set as the product of 6%/33.9%, burned area fraction of peat fire (s\ :sup:`-1`), and soil organic carbon (g C m\ :sup:`-2`). For boreal peat fires, the carbon emissions due to peat burning are set as 2.2 kg C m\ :sup:`-2` \ peat fire area (:ref:`Turetsky et al. 2002 <Turetskyetal2002>`). Biomass combustion and vegetation mortality in post-fire peatlands are set the same as section :numref:`Fire impact` for non-crop PFTs and as section :numref:`Agricultural fires` for crops PFTs.

.. _Fire trace gas and aerosol emissions:

Fire trace gas and aerosol emissions
-------------------------------------
CESM2 is the first Earth system model that can model the full coupling among fire, fire emissions, land, and atmosphere. CLM5, as the land component of CESM2, calculates the surface trace gas and aerosol emissions due to fire and fire emission heights, as the inputs of atmospheric chemistry model and aerosol model.

Emissions for trace gas and aerosol species x and the j-th PFT, :math:`E_{x,j}` (g species s\ :sup:`-1`), are given by

.. math::
   :label: 23.40

   E_{x,j} = EF_{x,j}\frac{\phi _{j} }{[C]}.

Here, :math:`EF_{x,j}` (g species (g dm)\ :sup:`-1`) is PFT-dependent emission factor scaled from biome-level values (Li et al., in prep, also used for FireMIP fire emissions data) by Dr. Val Martin and Dr. Li. :math:`[C]` = 0.5 (g C (g dm)\ :sup:`-1`) is a conversion factor from dry matter to carbon.

Emission height is PFT-dependent: 4.3 km for needleleaf tree PFTs, 3 km for other boreal and temperate tree PFTs, 2.5 km for tropical tree PFTs, 2 km for shrub PFTs, and 1 km for grass and crop PFTs. These values are compiled from earlier studies by Dr. Val Martin.

.. _Table PFT-specific combustion completeness and fire mortality factors:

.. table:: PFT-specific combustion completeness and fire mortality factors.

 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | PFT                              | *CC*\ :sub:`leaf`         | *CC*\ :sub:`stem`         | *CC*\ :sub:`root`         | *CC*\ :sub:`ts`         | *M*\ :sub:`leaf`         | *M*\ :sub:`livestem,1`       | *M*\ :sub:`deadstem`         | *M*\ :sub:`root`         | *M*\ :sub:`ts`         | *M*\ :sub:`livestem,2`       | :math:`\xi`\ :sub:`j`           |
 +==================================+===========================+===========================+===========================+=========================+==========================+==============================+==============================+==========================+========================+==============================+=================================+
 | NET Temperate                    | 0.80                      | 0.30                      | 0.00                      | 0.50                    | 0.80                     | 0.15                         | 0.15                         | 0.15                     | 0.50                   | 0.35                         | 0.15                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | NET Boreal                       | 0.80                      | 0.30                      | 0.00                      | 0.50                    | 0.80                     | 0.15                         | 0.15                         | 0.15                     | 0.50                   | 0.35                         | 0.15                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | NDT Boreal                       | 0.80                      | 0.30                      | 0.00                      | 0.50                    | 0.80                     | 0.15                         | 0.15                         | 0.15                     | 0.50                   | 0.35                         | 0.15                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BET Tropical                     | 0.80                      | 0.27                      | 0.00                      | 0.45                    | 0.80                     | 0.13                         | 0.13                         | 0.13                     | 0.45                   | 0.32                         | 0.13                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BET Temperate                    | 0.80                      | 0.27                      | 0.00                      | 0.45                    | 0.80                     | 0.13                         | 0.13                         | 0.13                     | 0.45                   | 0.32                         | 0.13                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BDT Tropical                     | 0.80                      | 0.27                      | 0.00                      | 0.45                    | 0.80                     | 0.10                         | 0.10                         | 0.10                     | 0.35                   | 0.25                         | 0.10                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BDT Temperate                    | 0.80                      | 0.27                      | 0.00                      | 0.45                    | 0.80                     | 0.10                         | 0.10                         | 0.10                     | 0.35                   | 0.25                         | 0.10                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BDT Boreal                       | 0.80                      | 0.27                      | 0.00                      | 0.45                    | 0.80                     | 0.13                         | 0.13                         | 0.13                     | 0.45                   | 0.32                         | 0.13                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BES Temperate                    | 0.80                      | 0.35                      | 0.00                      | 0.55                    | 0.80                     | 0.17                         | 0.17                         | 0.17                     | 0.55                   | 0.38                         | 0.17                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BDS Temperate                    | 0.80                      | 0.35                      | 0.00                      | 0.55                    | 0.80                     | 0.17                         | 0.17                         | 0.17                     | 0.55                   | 0.38                         | 0.17                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | BDS Boreal                       | 0.80                      | 0.35                      | 0.00                      | 0.55                    | 0.80                     | 0.17                         | 0.17                         | 0.17                     | 0.55                   | 0.38                         | 0.17                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | C\ :sub:`3` Grass Arctic         | 0.80                      | 0.80                      | 0.00                      | 0.80                    | 0.80                     | 0.20                         | 0.20                         | 0.20                     | 0.80                   | 0.60                         | 0.20                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | C\ :sub:`3` Grass                | 0.80                      | 0.80                      | 0.00                      | 0.80                    | 0.80                     | 0.20                         | 0.20                         | 0.20                     | 0.80                   | 0.60                         | 0.20                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | C\ :sub:`4` Grass                | 0.80                      | 0.80                      | 0.00                      | 0.80                    | 0.80                     | 0.20                         | 0.20                         | 0.20                     | 0.80                   | 0.60                         | 0.20                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+
 | Crop                             | 0.80                      | 0.80                      | 0.00                      | 0.80                    | 0.80                     | 0.20                         | 0.20                         | 0.20                     | 0.80                   | 0.60                         | 0.20                            |
 +----------------------------------+---------------------------+---------------------------+---------------------------+-------------------------+--------------------------+------------------------------+------------------------------+--------------------------+------------------------+------------------------------+---------------------------------+

Leaves (:math:`CC_{leaf}` ), stems (:math:`CC_{stem}` ), roots (:math:`CC_{root}` ), and transfer and storage carbon (:math:`CC_{ts}` ); mortality factors for leaves (:math:`M_{leaf}` ), live stems (:math:`M_{livestem,1}` ), dead stems (:math:`M_{deadstem}` ), roots (:math:`M_{root}` ), and transfer and storage carbon (:math:`M_{ts}` ) related to the carbon transfers from these pools to litter pool; mortality factors for live stems (:math:`M_{livestem,2}` ) related to the carbon transfer from live stems to dead stems; whole-plant mortality factor (:math:`\xi _{j}` ).
