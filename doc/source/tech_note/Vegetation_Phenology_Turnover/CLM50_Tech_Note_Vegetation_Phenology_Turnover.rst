.. _rst_Vegetation Phenology and Turnover:

Vegetation Phenology and Turnover
=================================

The CLM phenology model consists of several algorithms controlling the transfer of stored carbon and nitrogen out of storage pools for the display of new growth and into litter pools for losses of displayed growth. PFTs are classified into three distinct phenological types that are represented by separate algorithms: an evergreen type, for which some fraction of annual leaf growth persists in the displayed pool for longer than one year; a seasonal-deciduous type with a single growing season per year, controlled mainly by temperature and daylength; and a stress-deciduous type with the potential for multiple growing seasons per year, controlled by temperature and soil moisture conditions.

The three phenology types share a common set of control variables. The calculation of the phenology fluxes is generalized, operating identically for all three phenology types, given a specification of the common control variables. The following sections describe first the general flux parameterization, followed by the algorithms for setting the control parameters for the three phenology types.

General Phenology Flux Parameterization
--------------------------------------------

Fluxes of carbon and nitrogen from storage pools and into displayed tissue pools pass through a special transfer pool (denoted *\_xfer*), maintained as a separate state variable for each tissue type. Storage (*\_stor*) and transfer (*\_xfer*) pools are maintained separately to reduce the complexity of accounting for transfers into and out of storage over the course of a single growing season.

.. _Figure annual phenology cycle:

.. figure:: image1.png

 Example of annual phenology cycle for seasonal deciduous.

14.1.1 Onset Periods
^^^^^^^^^^^^^^^^^^^^

The deciduous phenology algorithms specify the occurrence of onset growth periods (Figure 14.1). Carbon fluxes from the transfer pools into displayed growth are calculated during these periods as:

.. math::
   :label: 20.1)

   CF_{leaf\_ xfer,leaf} =r_{xfer\_ on} CS_{leaf\_ xfer}

.. math::
   :label: 20.2)

   CF_{froot\_ xfer,froot} =r_{xfer\_ on} CS_{froot\_ xfer}

.. math::
   :label: 20.3)

   CF_{livestem\_ xfer,livestem} =r_{xfer\_ on} CS_{livestem\_ xfer}

.. math::
   :label: 20.4)

   CF_{deadstem\_ xfer,deadstem} =r_{xfer\_ on} CS_{deadstem\_ xfer}

.. math::
   :label: 20.5)

   CF_{livecroot\_ xfer,livecroot} =r_{xfer\_ on} CS_{livecroot\_ xfer}

.. math::
   :label: 20.6)

   CF_{deadcroot\_ xfer,deadcroot} =r_{xfer\_ on} CS_{deadcroot\_ xfer} ,

with corresponding nitrogen fluxes:

.. math::
   :label: 20.7)

   NF_{leaf\_ xfer,leaf} =r_{xfer\_ on} NS_{leaf\_ xfer}

.. math::
   :label: 20.8)

   NF_{froot\_ xfer,froot} =r_{xfer\_ on} NS_{froot\_ xfer}

.. math::
   :label: 20.9)

   NF_{livestem\_ xfer,livestem} =r_{xfer\_ on} NS_{livestem\_ xfer}

.. math::
   :label: 20.10)

   NF_{deadstem\_ xfer,deadstem} =r_{xfer\_ on} NS_{deadstem\_ xfer}

.. math::
   :label: 20.11)

   NF_{livecroot\_ xfer,livecroot} =r_{xfer\_ on} NS_{livecroot\_ xfer}

.. math::
   :label: 20.12)

   NF_{deadcroot\_ xfer,deadcroot} =r_{xfer\_ on} NS_{deadcroot\_ xfer} ,

where CF is the carbon flux, CS is stored carbon, NF is the nitrogen flux, NS is stored nitrogen, :math:`{r}_{xfer\_on}` (s\ :sup:`-1`) is a time-varying rate coefficient controlling flux out of the transfer pool:

.. math::
   :label: ZEqnNum852972

   r_{xfer\_ on} =\left\{\begin{array}{l} {{2\mathord{\left/ {\vphantom {2 t_{onset} }} \right.} t_{onset} } \qquad {\rm for\; }t_{onset} \ne \Delta t} \\ {{1\mathord{\left/ {\vphantom {1 \Delta t}} \right.} \Delta t} \qquad {\rm for\; }t_{onset} =\Delta t} \end{array}\right.

and *t*\ :sub:`onset` (s) is the number of seconds remaining in the current phenology onset growth period (Figure 14.1). The form of Eq. :eq:`ZEqnNum852972` produces a flux from the transfer pool which declines linearly over the onset growth period, approaching zero flux in the final timestep.

14.1.2 Offset Periods
^^^^^^^^^^^^^^^^^^^^^

The deciduous phenology algorithms also specify the occurrence of litterfall during offset periods. In contrast to the onset periods, only leaf and fine root state variables are subject to litterfall fluxes. Carbon fluxes from display pools into litter are calculated during these periods as:

.. math::
   :label: 20.14)

   CF_{leaf,litter}^{n} =\left\{\begin{array}{l} {CF_{leaf,litter}^{n-1} + r_{xfer\_ off} \left(CS_{leaf} -CF_{leaf,litter}^{n-1} {\kern 1pt} t_{offset} \right)\qquad {\rm for\; }t_{offset} \ne \Delta t}
   \\ {\left({CS_{leaf} \mathord{\left/ {\vphantom {CS_{leaf}  \Delta t}} \right.} \Delta t} \right)
   \left( 1-biofuel\_harvfrac  \right)
   +CF_{alloc,leaf} \qquad {\rm for\; }t_{offset} =\Delta t} \end{array}\right.

.. math::
   :label: 20.15)

   CF_{froot,litter}^{n} =\left\{\begin{array}{l} {CF_{froot,litter}^{n-1} +
   r_{xfer\_ off} \left(CS_{froot} -CF_{froot,litter}^{n-1} {\kern 1pt} t_{offset} \right)\qquad {\rm for\; }t_{offset} \ne \Delta t} \\ {\left({CS_{froot} \mathord{\left/ {\vphantom {CS_{froot}  \Delta t}} \right.} \Delta t} \right)+CF_{alloc,\, froot} \qquad \qquad \qquad {\rm for\; }t_{offset} =\Delta t} \end{array}\right.

.. math::
   :label: 20.16)

   r_{xfer\_ off} =\frac{2\Delta t}{t_{offset} ^{2} }

where superscripts *n* and *n-1* refer to fluxes on the current and previous timesteps, respectively. The rate coefficient :math:`{r}_{xfer\_off}` varies with time to produce a linearly increasing litterfall rate throughout the offset period. The :math:`biofuel\_harvfrac` (:numref:`Harvest to food and seed`) is the harvested fraction of aboveground biomass (leaf & livestem) for bioenergy crops. The special case for fluxes in the final litterfall timestep (:math:`{t}_{offset}` = :math:`\Delta t`\ ) ensures that all of the displayed growth is sent to the litter pools or biofuel feedstock pools. The fraction (:math:`biofuel\_harvfrac`) of leaf biomass going to the biofuel feedstock pools (Equation :eq:`25.9`) is defined in Table 26.3 and is only non-zero for prognostic crops. The remaining fraction of leaf biomass (:math:`1-biofuel\_harvfrac`) for deciduous plant types is sent to the litter pools. Similar modifications made for livestem carbon pools for prognostic crops can be found in section :numref:`Harvest to food and seed` in Equations :eq:`25.9`-:eq:`25.14`.

Corresponding nitrogen fluxes during litterfall take into account retranslocation of nitrogen out of the displayed leaf pool prior to litterfall (:math:`{NF}_{leaf,retrans}`, gN m\ :sup:`-2` s\ :sup:`-1`). Retranslocation of nitrogen out of fine roots is assumed to be negligible. The fluxes are:

.. math::
   :label: 20.17)

   NF_{leaf,litter} ={CF_{leaf,litter} \mathord{\left/ {\vphantom {CF_{leaf,litter}  CN_{leaf\_ litter} }} \right.} CN_{leaf\_ litter} }

.. math::
   :label: 20.18)

   NF_{froot,litter} ={CF_{leaf,litter} \mathord{\left/ {\vphantom {CF_{leaf,litter}  CN_{froot} }} \right.} CN_{froot} }

.. math::
   :label: 20.19)

   NF_{leaf,retrans} =\left({CF_{leaf,litter} \mathord{\left/ {\vphantom {CF_{leaf,litter}  CN_{leaf} }} \right.} CN_{leaf} } \right)-NF_{leaf,litter} .

where CN is C:N.

14.1.3 Background Onset Growth
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The stress-deciduous phenology algorithm includes a provision for the case when stress signals are absent, and the vegetation shifts from a deciduous habit to an evergreen habit, until the next occurrence of an offset stress trigger. In that case, the regular onset flux mechanism is switched off and a background onset growth algorithm is invoked (:math:`{r}_{bgtr} > 0`). During this period, small fluxes of carbon and nitrogen from the storage pools into the associated transfer pools are calculated on each time step, and the entire contents of the transfer pool are added to the associated displayed growth pool on each time step. The carbon fluxes from transfer to display pools under these conditions are:

.. math::
   :label: 20.20)

   CF_{leaf\_ xfer,leaf} ={CS_{leaf\_ xfer} \mathord{\left/ {\vphantom {CS_{leaf\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.21)

   CF_{froot\_ xfer,froot} ={CS_{froot\_ xfer} \mathord{\left/ {\vphantom {CS_{froot\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.22)

   CF_{livestem\_ xfer,livestem} ={CS_{livestem\_ xfer} \mathord{\left/ {\vphantom {CS_{livestem\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.23)

   CF_{deadstem\_ xfer,deadstem} ={CS_{deadstem\_ xfer} \mathord{\left/ {\vphantom {CS_{deadstem\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.24)

   CF_{livecroot\_ xfer,livecroot} ={CS_{livecroot\_ xfer} \mathord{\left/ {\vphantom {CS_{livecroot\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.25)

   CF_{deadcroot\_ xfer,deadcroot} ={CS_{deadcroot\_ xfer} \mathord{\left/ {\vphantom {CS_{deadcroot\_ xfer}  \Delta t}} \right.} \Delta t} ,

and the corresponding nitrogen fluxes are:

.. math::
   :label: 20.26)

   NF_{leaf\_ xfer,leaf} ={NS_{leaf\_ xfer} \mathord{\left/ {\vphantom {NS_{leaf\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.27)

   NF_{froot\_ xfer,froot} ={NS_{froot\_ xfer} \mathord{\left/ {\vphantom {NS_{froot\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.28)

   NF_{livestem\_ xfer,livestem} ={NS_{livestem\_ xfer} \mathord{\left/ {\vphantom {NS_{livestem\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.29)

   NF_{deadstem\_ xfer,deadstem} ={NS_{deadstem\_ xfer} \mathord{\left/ {\vphantom {NS_{deadstem\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.30)

   NF_{livecroot\_ xfer,livecroot} ={NS_{livecroot\_ xfer} \mathord{\left/ {\vphantom {NS_{livecroot\_ xfer}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.31)

   NF_{deadcroot\_ xfer,deadcroot} ={NS_{deadcroot\_ xfer} \mathord{\left/ {\vphantom {NS_{deadcroot\_ xfer}  \Delta t}} \right.} \Delta t} .

14.1.4 Background Litterfall
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Both evergreen and stress-deciduous phenology algorithms can specify a litterfall flux that is not associated with a specific offset period, but which occurs instead at a slow rate over an extended period of time, referred to as background litterfall. For evergreen types the background litterfall is the only litterfall flux. For stress-deciduous types either the offset period litterfall or the background litterfall mechanism may be active, but not both at once. Given a specification of the background litterfall rate (:math:`{r}_{bglf}`, s\ :sup:`-1`), litterfall carbon fluxes are calculated as

.. math::
   :label: 20.32)

   CF_{leaf,litter} =r_{bglf} CS_{leaf}

.. math::
   :label: 20.33)

   CS_{froot,litter} =r_{bglf} CS_{froot} ,

with corresponding nitrogen litterfall and retranslocation fluxes:

.. math::
   :label: 20.34)

   NF_{leaf,litter} ={CF_{leaf,litter} \mathord{\left/ {\vphantom {CF_{leaf,litter}  CN_{leaf\_ litter} }} \right.} CN_{leaf\_ litter} }

.. math::
   :label: 20.35)

   NF_{froot,litter} ={CF_{froot,litter} \mathord{\left/ {\vphantom {CF_{froot,litter}  CN_{froot} }} \right.} CN_{froot} }

.. math::
   :label: 20.36)

   NF_{leaf,retrans} =\left({CF_{leaf,litter} \mathord{\left/ {\vphantom {CF_{leaf,litter}  CN_{leaf} }} \right.} CN_{leaf} } \right)-NF_{leaf,litter} .

14.1.5 Livewood Turnover
^^^^^^^^^^^^^^^^^^^^^^^^

The conceptualization of live wood vs. dead wood fractions for stem and coarse root pools is intended to capture the difference in maintenance respiration rates between these two physiologically distinct tissue types. Unlike displayed pools for leaf and fine root, which are lost to litterfall, live wood cells reaching the end of their lifespan are retained as a part of the dead woody structure of stems and coarse roots. A mechanism is therefore included in the phenology routine to effect the transfer of live wood to dead wood pools, which also takes into account the different nitrogen concentrations typical of these tissue types.

A live wood turnover rate (:math:`{r}_{lwt}`, s\ :sup:`-1`) is defined as

.. math::
   :label: 20.37)

   r_{lwt} ={p_{lwt} \mathord{\left/ {\vphantom {p_{lwt}  \left(365\cdot 86400\right)}} \right.} \left(365\cdot 86400\right)}

where :math:`{p}_{lwt} = 0.7` is the assumed annual live wood turnover fraction. Carbon fluxes from live to dead wood pools are:

.. math::
   :label: 20.38)

   CF_{livestem,deadstem} =CS_{livestem} r_{lwt}

.. math::
   :label: 20.39)

   CF_{livecroot,deadcroot} =CS_{livecroot} r_{lwt} ,

and the associated nitrogen fluxes, including retranslocation of nitrogen out of live wood during turnover, are:

.. math::
   :label: 20.40)

   NF_{livestem,deadstem} ={CF_{livestem,deadstem} \mathord{\left/ {\vphantom {CF_{livestem,deadstem}  CN_{dw} }} \right.} CN_{dw} }

.. math::
   :label: 20.41)

   NF_{livestem,retrans} =\left({CF_{livestem,deadstem} \mathord{\left/ {\vphantom {CF_{livestem,deadstem}  CN_{lw} }} \right.} CN_{lw} } \right)-NF_{livestem,deadstem}

.. math::
   :label: 20.42)

   NF_{livecroot,deadcroot} ={CF_{livecroot,deadcroot} \mathord{\left/ {\vphantom {CF_{livecroot,deadcroot}  CN_{dw} }} \right.} CN_{dw} }

.. math::
   :label: 20.43)

   NF_{livecroot,retrans} =\left({CF_{livecroot,deadcroot} \mathord{\left/ {\vphantom {CF_{livecroot,deadcroot}  CN_{lw} }} \right.} CN_{lw} } \right)-NF_{livecroot,deadcroot} .

Evergreen Phenology
------------------------

The evergreen phenology algorithm is by far the simplest of the three possible types. It is assumed for all evergreen types that all carbon and nitrogen allocated for new growth in the current timestep goes immediately to the displayed growth pools (i.e. f\ :math:`{f}_{cur} = 1.0` (Chapter 13)). As such, there is never an accumulation of carbon or nitrogen in the storage or transfer pools, and so the onset growth and background onset growth mechanisms are never invoked for this type. Litterfall is specified to occur only through the background litterfall mechanism – there are no distinct periods of litterfall for evergreen types, but rather a continuous (slow) shedding of foliage and fine roots. This is an obvious area for potential improvements in the model, since it is known, at least for evergreen needleleaf trees in the temperate and boreal zones, that there are distinct periods of higher and lower leaf litterfall (Ferrari, 1999; Gholz et al., 1985). The rate of background litterfall (:math:`{r}_{bglf}`, section 14.1.4) depends on the specified leaf longevity (:math:`\tau_{leaf}`\, y), as

.. math::
   :label: 20.44)

   r_{bglf} =\frac{1}{\tau _{leaf} \cdot 365\cdot 86400} .

Seasonal-Deciduous Phenology
---------------------------------

The seasonal-deciduous phenology algorithm derives directly from the treatment used in the offline model Biome-BGC v. 4.1.2, (Thornton et al., 2002), which in turn is based on the parameterizations for leaf onset and offset for temperate deciduous broadleaf forest from White et al. (1997). Initiation of leaf onset is triggered when a common degree-day summation exceeds a critical value, and leaf litterfall is initiated when daylength is shorter than a critical value. Because of the dependence on daylength, the seasonal deciduous phenology algorithm is only valid for latitudes outside of the tropical zone, defined here as :math:`\left|{\rm latitude}\right|>19.5{\rm {}^\circ }`. Neither the background onset nor background litterfall mechanism is invoked for the seasonal-deciduous phenology algorithm. The algorithm allows a maximum of one onset period and one offset period each year.

The algorithms for initiation of onset and offset periods use the winter and summer solstices as coordination signals. The period between winter and summer solstice is identified as :math:`{dayl}_{n} > {dayl}_{n-1}`, and the period between summer and winter solstice is identified as :math:`{dayl}_{n} < {dayl}_{n-1}`, where :math:`{dayl}_{n}` and :math:`{dayl}_{n-1}` are the day length(s) calculated for the current and previous timesteps, respectively, using

.. math::
   :label: 20.45)

   dayl=2\cdot 13750.9871\cdot acos\left(\frac{-\sin (lat)\sin (decl)}{\cos (lat)\cos (decl)} \right),

where *lat* and *decl* are the latitude and solar declination (radians), respectively, and the factor 13750.9871 is the number of seconds per radian of hour-angle.

14.3.1 Seasonal-Deciduous Onset Trigger
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The onset trigger for the seasonal-deciduous phenology algorithm is based on an accumulated growing-degree-day approach (White et al., 1997). The growing-degree-day summation (:math:`{GDD}_{sum}`) is initiated ( :math:`{GDD}_{sum} = 0`) when the phenological state is dormant and the model timestep crosses the winter solstice. Once these conditions are met, :math:`{GDD}_{sum}` is updated on each timestep as

.. math::
   :label: ZEqnNum510730

   GDD_{sum}^{n} =\left\{\begin{array}{l} {GDD_{sum}^{n-1} +\left(T_{s,3} -TKFRZ\right)f_{day} \qquad {\rm for\; }T_{s,3} >TKFRZ} \\ {GDD_{sum}^{n-1} \qquad \qquad \qquad {\rm for\; }T_{s,3} \le TKFRZ} \end{array}\right.

where :math:`{T}_{s,3}` (K) is the temperature of the third soil layer, and :math:`f_{day} ={\Delta t\mathord{\left/ {\vphantom {\Delta t 86400}} \right.} 86400}`. The onset period is initiated if :math:`GDD_{sum} >GDD_{sum\_ crit}`, where

.. math::
   :label: ZEqnNum598907

   GDD_{sum\_ crit} =\exp \left(4.8+0.13{\kern 1pt} \left(T_{2m,ann\_ avg} -TKFRZ\right)\right)

and where :math:`{T}_{2m,ann\_avg}` (K) is the annual average of the 2m air temperature, and TKFRZ is the freezing point of water (273.15 K). The following control variables are set when a new onset growth period is initiated:

.. math::
   :label: 20.48)

   GDD_{sum} =0

.. math::
   :label: 20.49)

   t_{onset} =86400\cdot n_{days\_ on} ,

where :math:`{n}_{days\_on}` is set to a constant value of 30 days. Fluxes from storage into transfer pools occur in the timestep when a new onset growth period is initiated. Carbon fluxes are:

.. math::
   :label: ZEqnNum904388

   CF_{leaf\_ stor,leaf\_ xfer} ={f_{stor,xfer} CS_{leaf\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{leaf\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.51)

   CF_{froot\_ stor,froot\_ xfer} ={f_{stor,xfer} CS_{froot\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{froot\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.52)

   CF_{livestem\_ stor,livestem\_ xfer} ={f_{stor,xfer} CS_{livestem\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{livestem\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.53)

   CF_{deadstem\_ stor,deadstem\_ xfer} ={f_{stor,xfer} CS_{deadstem\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{deadstem\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.54)

   CF_{livecroot\_ stor,livecroot\_ xfer} ={f_{stor,xfer} CS_{livecroot\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{livecroot\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.55)

   CF_{deadcroot\_ stor,deadcroot\_ xfer} ={f_{stor,xfer} CS_{deadcroot\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{deadcroot\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: ZEqnNum195642

   CF_{gresp\_ stor,gresp\_ xfer} ={f_{stor,xfer} CS_{gresp\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} CS_{gresp\_ stor}  \Delta t}} \right.} \Delta t}

and the associated nitrogen fluxes are:

.. math::
   :label: ZEqnNum812152

   NF_{leaf\_ stor,leaf\_ xfer} ={f_{stor,xfer} NS_{leaf\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} NS_{leaf\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.58)

   NF_{froot\_ stor,froot\_ xfer} ={f_{stor,xfer} NS_{froot\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} NS_{froot\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.59)

   NF_{livestem\_ stor,livestem\_ xfer} ={f_{stor,xfer} NS_{livestem\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} NS_{livestem\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.60)

   NF_{deadstem\_ stor,deadstem\_ xfer} ={f_{stor,xfer} NS_{deadstem\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} NS_{deadstem\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 20.61)

   NF_{livecroot\_ stor,livecroot\_ xfer} ={f_{stor,xfer} NS_{livecroot\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} NS_{livecroot\_ stor}  \Delta t}} \right.} \Delta t}

.. math::
   :label: ZEqnNum605338

   NF_{deadcroot\_ stor,deadcroot\_ xfer} ={f_{stor,xfer} NS_{deadcroot\_ stor} \mathord{\left/ {\vphantom {f_{stor,xfer} NS_{deadcroot\_ stor}  \Delta t}} \right.} \Delta t}

where :math:`{f}_{stor,xfer}` is the fraction of current storage pool moved into the transfer pool for display over the incipient onset period. This fraction is set to 0.5, based on the observation that seasonal deciduous trees are capable of replacing their canopies from storage reserves in the event of a severe early-season disturbance such as frost damage or defoliation due to insect herbivory.

If the onset criterion (:math:`{GDD}_{sum} > {GDD}_{sum\_crit}`) is not met before the summer solstice, then :math:`{GDD}_{sum}` is set to 0.0 and the growing-degree-day accumulation will not start again until the following winter solstice. This mechanism prevents the initiation of very short growing seasons late in the summer in cold climates. The onset counter is decremented on each time step after initiation of the onset period, until it reaches zero, signaling the end of the onset period:

.. math::
   :label: 20.63)

   t_{onfset}^{n} =t_{onfset}^{n-1} -\Delta t

14.3.2 Seasonal-Deciduous Offset Trigger
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the completion of an onset period, and once past the summer solstice, the offset (litterfall) period is triggered when daylength is shorter than 39300 s. The offset counter is set at the initiation of the offset period: :math:`t_{offset} =86400\cdot n_{days\_ off}`, where :math:`{n}_{days\_off}` is set to a constant value of 15 days. The offset counter is decremented on each time step after initiation of the offset period, until it reaches zero, signaling the end of the offset period:

.. math::
   :label: 20.64)

   t_{offset}^{n} =t_{offset}^{n-1} -\Delta t

Stress-Deciduous Phenology
-------------------------------

The stress-deciduous phenology algorithm was developed specifically for the CLM based in part on the grass phenology model proposed by White et al. (1997). The algorithm handles phenology for vegetation types such as grasses and tropical drought-deciduous trees that respond to both cold and drought-stress signals, and that can have multiple growing seasons per year. The algorithm also allows for the possibility that leaves might persist year-round in the absence of a suitable stress trigger. In that case the phenology switches to an evergreen habit, maintaining a marginally-deciduous leaf longevity (one year) until the occurrence of the next stress trigger.

14.4.1 Stress-Deciduous Onset Triggers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In climates that are warm year-round, onset triggering depends on soil water availability. At the beginning of a dormant period (end of previous offset period), an accumulated soil water index (:math:`{SWI}_{sum}`, d) is initialized (:math:`{SWI}_{sum} = 0`), with subsequent accumulation calculated as:

.. math::
   :label: ZEqnNum503826

   SWI_{sum}^{n} =\left\{\begin{array}{l} {SWI_{sum}^{n-1} +f_{day} \qquad {\rm for\; }\Psi _{s,3} \ge \Psi _{onset} } \\ {SWI_{sum}^{n-1} \qquad \qquad {\rm for\; }\Psi _{s,3} <\Psi _{onset} } \end{array}\right.

where :math:`\Psi`\ :sub:`s,3` is the soil water potential (MPa) in the third soil layer and :math:`{\Psi}_{onset} = -0.6 MPa` is the onset soil water potential threshold. Onset triggering is possible once :math:`{SWI}_{sum} > 15`. To avoid spurious onset triggering due to soil moisture in the third soil layer exceeding the threshold due only to soil water suction of water from deeper in the soil column, an additional precipitation trigger is included which requires at least 20 mm of rain over the previous 10 days :ref:`(Dahlin et al., 2015) <Dahlinetal2015>`. If the cold climate growing degree-day accumulator is not active at the time when the soil moisture and precipitation thresholds are reached (see below), and if the daylength is greater than 6 hours, then onset is triggered. Except as noted below, :math:`{SWI}_{sum}` continues to accumulate according to Eq. :eq:`ZEqnNum503826` during the dormant period if the daylength criterion prevents onset triggering, and onset is then triggered at the timestep when daylength exceeds 6 hours.

In climates with a cold season, onset triggering depends on both accumulated soil temperature summation and adequate soil moisture. At the beginning of a dormant period a freezing day accumulator (:math:`{FD}_{sum}`, d) is initialized (:math:`{FD}_{sum} = 0`), with subsequent accumulation calculated as:

.. math::
   :label: 20.66)

   FD_{sum}^{n} =\left\{\begin{array}{l} {FD_{sum}^{n-1} +f_{day} \qquad {\rm for\; }T_{s,3} >TKFRZ} \\ {FD_{sum}^{n-1} \qquad \qquad {\rm for\; }T_{s,3} \le TKFRZ} \end{array}\right. .

If :math:`{FD}_{sum} > 15` during the dormant period, then a cold-climate onset triggering criterion is introduced, following exactly the growing degree-day summation (:math:`{GDD}_{sum}`) logic of Eqs. :eq:`ZEqnNum510730` and :eq:`ZEqnNum598907`. At that time :math:`{SWI}_{sum}` is reset (:math:`{SWI}_{sum} = 0`). Onset triggering under these conditions depends on meeting all three of the following criteria: :math:`{SWI}_{sum} > 15`, :math:`{GDD}_{sum} > {GDD}_{sum\_crit}`, and daylength greater than 6 hrs.

The following control variables are set when a new onset growth period is initiated: :math:`{SWI}_{sum} = 0`, :math:`{FD}_{sum} = 0`, :math:`{GDD}_{sum} = 0`, :math:`{n}_{days\_active} = 0`, and :math:`t_{onset} = 86400\cdot n_{days\_ on}`, where :math:`{n}_{days\_on}` is set to a constant value of 30 days. Fluxes from storage into transfer pools occur in the timestep when a new onset growth period is initiated, and are handled identically to Eqs. :eq:`ZEqnNum904388` -:eq:`ZEqnNum195642` for carbon fluxes, and to Eqs. :eq:`ZEqnNum812152` - :eq:`ZEqnNum605338` for nitrogen fluxes. The onset counter is decremented on each time step after initiation of the onset period, until it reaches zero, signaling the end of the onset period:

.. math::
   :label: 20.67)

   t_{onfset}^{n} =t_{onfset}^{n-1} -\Delta t

14.4.2 Stress-Deciduous Offset Triggers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Any one of the following three conditions is sufficient to initiate an offset period for the stress-deciduous phenology algorithm: sustained period of dry soil, sustained period of cold temperature, or daylength shorter than 6 hours. Offset triggering due to dry soil or cold temperature conditions is only allowed once the most recent onset period is complete. Dry soil condition is evaluated with an offset soil water index accumulator (:math:`{OSWI}_{sum}`, d). To test for a sustained period of dry soils, this control variable can increase or decrease, as follows:

.. math::
   :label: 20.68)

   OSWI_{sum}^{n} =\left\{\begin{array}{l} {OSWI_{sum}^{n-1} +f_{day} \qquad \qquad \qquad {\rm for\; }\Psi _{s,3} \le \Psi _{offset} } \\ {{\rm max}\left(OSWI_{sum}^{n-1} -f_{day} ,0\right)\qquad {\rm for\; }\Psi _{s,3} >\Psi _{onset} } \end{array}\right.

where :math:`{\Psi}_{offset} = -0.8 MPa` is the offset soil water potential threshold. An offset period is triggered if the previous onset period is complete and :math:`{OSWI}_{sum}` :math:`\mathrm{\ge}` :math:`{OSWI}_{sum\_crit}`, where :math:`{OSWI}_{sum\_crit} = 15`.

The cold temperature trigger is calculated with an offset freezing day accumulator (:math:`{OFD}_{sum}`, d). To test for a sustained period of cold temperature, this variable can increase or decrease, as follows:

.. math::
   :label: 20.69)

   OFD_{sum}^{n} =\left\{\begin{array}{l} {OFD_{sum}^{n-1} +f_{day} \qquad \qquad \qquad {\rm for\; }T_{s,3} \le TKFRZ} \\ {{\rm max}\left(OFD_{sum}^{n-1} -f_{day} ,0\right)\qquad \qquad {\rm for\; }T_{s,3} >TKFRZ} \end{array}\right.

An offset period is triggered if the previous onset period is complete and :math:`{OFD}_{sum} > {OFD}_{sum\_crit}`, where :math:`{OFD}_{sum\_crit} = 15`.

The offset counter is set at the initiation of the offset period: :math:`t_{offset} =86400\cdot n_{days\_ off}`, where :math:`{n}_{days\_off}` is set to a constant value of 15 days. The offset counter is decremented on each time step after initiation of the offset period, until it reaches zero, signaling the end of the offset period:

.. math::
   :label: 20.70)

   t_{offset}^{n} =t_{offset}^{n-1} -\Delta t

14.4.3 Stress-Deciduous: Long Growing Season
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under conditions when the stress-deciduous conditions triggering offset are not met for one year or longer, the stress-deciduous algorithm shifts toward the evergreen behavior. This can happen in cases where a stress-deciduous vegetation type is assigned in a climate where suitably strong stresses occur less frequently than once per year. This condition is evaluated by tracking the number of days since the beginning of the most recent onset period (:math:`{n}_{days\_active}`, d). At the end of an offset period :math:`{n}_{days\_active}` is reset to 0. A long growing season control variable (*LGS*, range 0 to 1) is calculated as:

.. math::
   :label: 20.71)

   LGS=\left\{\begin{array}{l} {0\qquad \qquad \qquad {\rm for\; }n_{days\_ active} <365} \\ {\left({n_{days\_ active} \mathord{\left/ {\vphantom {n_{days\_ active}  365}} \right.} 365} \right)-1\qquad {\rm for\; }365\le n_{days\_ active} <730} \\ {1\qquad \qquad \qquad {\rm for\; }n_{days\_ active} \ge 730} \end{array}\right. .

The rate coefficient for background litterfall (:math:`{r}_{bglf}`, s\ :sup:`-1`) is calculated as a function of *LGS*:

.. math::
   :label: 20.72)

   r_{bglf} =\frac{LGS}{\tau _{leaf} \cdot 365\cdot 86400}

where :math:`{\tau}_{leaf}` is the leaf longevity. The result is a shift to continuous litterfall as :math:`{n}_{days\_active}` increases from 365 to 730. When a new offset period is triggered :math:`{r}_{bglf}` is set to 0.

The rate coefficient for background onset growth from the transfer pools ( :math:`{r}_{bgtr}`, s\ :sup:`-1`) also depends on *LGS*, as:

.. math::
   :label: 20.73)

   r_{bgtr} =\frac{LGS}{365\cdot 86400} .

On each timestep with :math:`{r}_{bgtr}` :math:`\neq` 0, carbon fluxes from storage to transfer pools are calculated as:

.. math::
   :label: 20.74)

   CF_{leaf\_ stor,leaf\_ xfer} =CS_{leaf\_ stor} r_{bgtr}

.. math::
   :label: 20.75)

   CF_{froot\_ stor,froot\_ xfer} =CS_{froot\_ stor} r_{bgtr}

.. math::
   :label: 20.76)

   CF_{livestem\_ stor,livestem\_ xfer} =CS_{livestem\_ stor} r_{bgtr}

.. math::
   :label: 20.77)

   CF_{deadstem\_ stor,deadstem\_ xfer} =CS_{deadstem\_ stor} r_{bgtr}

.. math::
   :label: 20.78)

   CF_{livecroot\_ stor,livecroot\_ xfer} =CS_{livecroot\_ stor} r_{bgtr}

.. math::
   :label: 20.79)

   CF_{deadcroot\_ stor,deadcroot\_ xfer} =CS_{deadcroot\_ stor} r_{bgtr} ,

with corresponding nitrogen fluxes:

.. math::
   :label: 20.80)

   NF_{leaf\_ stor,leaf\_ xfer} =NS_{leaf\_ stor} r_{bgtr}

.. math::
   :label: 20.81)

   NF_{froot\_ stor,froot\_ xfer} =NS_{froot\_ stor} r_{bgtr}

.. math::
   :label: 20.82)

   NF_{livestem\_ stor,livestem\_ xfer} =NS_{livestem\_ stor} r_{bgtr}

.. math::
   :label: 20.83)

   NF_{deadstem\_ stor,deadstem\_ xfer} =NS_{deadstem\_ stor} r_{bgtr}

.. math::
   :label: 20.84)

   NF_{livecroot\_ stor,livecroot\_ xfer} =NS_{livecroot\_ stor} r_{bgtr}

.. math::
   :label: 20.85)

   NF_{deadcroot\_ stor,deadcroot\_ xfer} =NS_{deadcroot\_ stor} r_{bgtr} .

The result, in conjunction with the treatment of background onset growth, is a shift to continuous transfer from storage to display pools at a rate that would result in complete turnover of the storage pools in one year at steady state, once *LGS* reaches 1 (i.e. after two years without stress-deciduous offset conditions). If and when conditions cause stress-deciduous triggering again, :math:`{r}_{bgtr}` is rest to 0.

Litterfall Fluxes Merged to the Column Level
-------------------------------------------------

CLM uses three litter pools, defined on the basis of commonly measured chemical fractionation of fresh litter into labile (LIT1 = hot water and alcohol soluble fraction), cellulose/hemicellulose (LIT2 = acid soluble fraction) and remaining material, referred to here for convenience as lignin (LIT3 = acid insoluble fraction) (Aber et al., 1990; Taylor et al., 1989). While multiple plant functional types can coexist on a single CLM soil column, each soil column includes a single instance of the litter pools. Fluxes entering the litter pools due to litterfall are calculated using a weighted average of the fluxes originating at the PFT level. Carbon fluxes are calculated as:

.. math::
   :label: 20.86)

   CF_{leaf,lit1} =\sum _{p=0}^{npfts}CF_{leaf,litter} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 20.87)

   CF_{leaf,lit2} =\sum _{p=0}^{npfts}CF_{leaf,litter} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 20.88)

   CF_{leaf,lit3} =\sum _{p=0}^{npfts}CF_{leaf,litter} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 20.89)

   CF_{froot,lit1} =\sum _{p=0}^{npfts}CF_{froot,litter} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 20.90)

   CF_{froot,lit2} =\sum _{p=0}^{npfts}CF_{froot,litter} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 20.91)

   CF_{froot,lit3} =\sum _{p=0}^{npfts}CF_{froot,litter} f_{lig\_ froot,p} wcol_{p}  ,

where :math:`{f}_{lab\_leaf,p}`, :math:`{f}_{cel\_leaf,p}`, and :math:`{f}_{lig\_leaf,p}` are the labile, cellulose/hemicellulose, and lignin fractions of leaf litter for PFT *p*, :math:`{f}_{lab\_froot,p}`, :math:`{f}_{cel\_froot,p}`, and :math:`{f}_{lig\_froot,p}` are the labile, cellulose/hemicellulose, and lignin fractions of fine root litter for PFT *p*, :math:`{wtcol}_{p}` is the weight relative to the column for PFT *p*, and *p* is an index through the plant functional types occurring on a column. Nitrogen fluxes to the litter pools are assumed to follow the C:N of the senescent tissue, and so are distributed using the same fractions used for carbon fluxes:

.. math::
   :label: 20.92)

   NF_{leaf,lit1} =\sum _{p=0}^{npfts}NF_{leaf,litter} f_{lab\_ leaf,p} wcol_{p}

.. math::
   :label: 20.93)

   NF_{leaf,lit2} =\sum _{p=0}^{npfts}NF_{leaf,litter} f_{cel\_ leaf,p} wcol_{p}

.. math::
   :label: 20.94)

   NF_{leaf,lit3} =\sum _{p=0}^{npfts}NF_{leaf,litter} f_{lig\_ leaf,p} wcol_{p}

.. math::
   :label: 20.95)

   NF_{froot,lit1} =\sum _{p=0}^{npfts}NF_{froot,litter} f_{lab\_ froot,p} wcol_{p}

.. math::
   :label: 20.96)

   NF_{froot,lit2} =\sum _{p=0}^{npfts}NF_{froot,litter} f_{cel\_ froot,p} wcol_{p}

.. math::
   :label: 20.97)

   NF_{froot,lit3} =\sum _{p=0}^{npfts}NF_{froot,litter} f_{lig\_ froot,p} wcol_{p}  .
