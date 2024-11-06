.. _rst_External Nitrogen Cycle:

External Nitrogen Cycle
===========================

.. _Summary of CLM5.0 updates relative to CLM4.5:

Summary of CLM5.0 updates relative to CLM4.5
-----------------------------------------------------

We describe external inputs to the nitrogen cycle in CLM5.0.  Much of the following information appeared in the CLM4.5 Technical Note (:ref:`Oleson et al. 2013 <Olesonetal2013>`) as well as :ref:`Koven et al. (2013) <Kovenetal2013>`.

CLM5.0 includes the following changes to terrestrial nitrogen inputs:

- Time varrying deposition of reactive nitrogen. In off-line runs this changes monthly. In coupled simulations N deposition is passed at the coupling timestep (e.g., half-hourly).

- Asymbiotic (or free living) N fixation is a function of evapotranspiration and is added to the inorganic nitrogen (NH\ :sub:`4`\ :sup:`+`) pool (described below).

- Symbiotic N fixation is handled by the FUN model (chapter :numref:`rst_FUN`) and is passed straight to the plant, not the mineral nitrogen pool.

Overview
-----------------------------------------------------

In addition to the relatively rapid cycling of nitrogen within the plant – litter – soil organic matter system, CLM also represents several processes which couple the internal nitrogen cycle to external sources and sinks. Inputs of new mineral nitrogen are from atmospheric deposition and biological nitrogen fixation. Losses of mineral nitrogen are due to nitrification, denitrification, leaching, and losses in fire. While the short-term dynamics of nitrogen limitation depend on the behavior of the internal nitrogen cycle, establishment of total ecosystem nitrogen stocks depends on the balance between sources and sinks in the external nitrogen cycle (:ref:`Thomas et al. 2015 <Thomasetal2015>`).

As with CLM4.5, CLM5.0 represents inorganic N transformations based on the Century N-gas model; this includes separate NH\ :sub:`4`\ :sup:`+` and NO\ :sub:`3`\ :sup:`-` pools, as well as environmentally controlled nitrification and denitrification rates that is described below.

Atmospheric Nitrogen Deposition
------------------------------------

CLM uses a single variable to represent the total deposition of mineral nitrogen onto the land surface, combining wet and dry deposition of NO\ :sub:`y` and NH\ :sub:`x` as a single flux (:math:`{NF}_{ndep\_sminn}`, gN m\ :sup:`-2` s\ :sup:`-1`). This flux is intended to represent total reactive nitrogen deposited to the land surface which originates from the following natural and anthropogenic sources (Galloway et al. 2004): formation of NO\ :sub:`x` during lightning, NO\ :math:`{}_{x }`\ and NH\ :sub:`3` emission from wildfire, NO\ :sub:`x` emission from natural soils, NH\ :sub:`3` emission from natural soils, vegetation, and wild animals, NO\ :sub:`x` and NH\ :sub:`3` emission during fossil fuel combustion (both thermal and fuel NO\ :sub:`x` production), NO\ :sub:`x` and NH\ :sub:`3` emission from other industrial processes, NO\ :sub:`x` and NH\ :sub:`3` emission from fire associated with deforestation, NO\ :sub:`x` and NH\ :sub:`3` emission from agricultural burning, NO\ :sub:`x` emission from agricultural soils, NH\ :sub:`3` emission from agricultural crops, NH\ :sub:`3` emission from agricultural animal waste, and NH\ :sub:`3` emission from human waste and waste water. The deposition flux is provided as a spatially and (potentially) temporally varying dataset (see section :numref:`Atmospheric Coupling` for a description of the default input dataset).

The nitrogen deposition flux is assumed to enter the NH\ :sub:`4`\ :sup:`+` pool, and is vertically distributed throughout the soil profile. Although N deposition inputs include both oxidized and reduced forms, CLM5 only reads in total N deposition. This approach is held over from CLM4.0, which only represented a single mineral nitrogen pool, however, real pathways for wet and dry nitrogen deposition can be more complex than currently represented in the CLM5.0, including release from melting snowpack and direct foliar uptake of deposited NO\ :sub:`y` (:ref:`Tye et al. 2005 <Tyeetal2005>`; :ref:`Vallano and Sparks, 2007 <VallanoSparks2007>`).

In offline (uncoupled) CLM5.0 simulations monthly estimates of N deposition are provided, as opposed to decadal files supplied with previous versions of the model. In coupled simulations, N depositions fluxes are passed to the land model at the frequency of the time step (every half hour) through the coupler.

Biological Nitrogen Fixation
---------------------------------

The fixation of new reactive nitrogen from atmospheric N\ :sub:`2` by soil microorganisms is an important component of both preindustrial and modern-day nitrogen budgets, but a mechanistic understanding of global-scale controls on biological nitrogen fixation (BNF) is still only poorly developed (:ref:`Cleveland et al. 1999 <Clevelandetal1999>`; :ref:`Galloway et al. 2004 <Gallowayetal2004>`). CLM5.0 uses the FUN model (chapter :numref:`rst_FUN`) to calculate the carbon cost and nitrogen acquired through symbotic nitrogen fixation. This nitrogen is immediately available to plants.

:ref:`Cleveland et al. (1999) <Clevelandetal1999>` suggested an empirical relationships that predicts BNF as a function of either evapotranspiration rate or net primary productivity for natural vegetation. CLM5.0 adopts the evapotranspiration approach to calculate asymbiotic, or free-living, N fixation. This function has been modified from the :ref:`Cleveland et al. (1999) <Clevelandetal1999>` estimates to provide lower estimate of free-living nitrogen fixation in CLM5.0 (:math:`{CF}_{ann\_ET}`, mm yr\ :sup:`-1`). This moves away from the NPP approach used in CLM4.0 and 4.5 and avoids unrealistically increasing freeliving rates of N fixation under global change scenarios (:ref:`Wieder et al. 2015 <Wiederetal2015>` The expression used is:

.. math::
   :label: 22.1)

   NF_{nfix,sminn} ={0.0006\left(0.0117+CF_{ann\_ ET}\right)\mathord{\left/ {\vphantom {0.0006\left(0.0117+ CF_{ann\_ ET}\right) \left(86400\cdot 365\right)}} \right.} \left(86400\cdot 365\right)}

Where :math:`{NF}_{nfix,sminn}` (gN m\ :sup:`-2` s\ :sup:`-1`) is the rate of free-living nitrogen fixation in :numref:`Figure Biological nitrogen fixation`.

.. _Figure Biological nitrogen fixation:

.. figure:: image1.png

 Free-living nitrogen fixation as a function of annual evapotranspiration. Results here show annual N inputs from free-living N fixations, but the model actually calculates inputs on a per second basis.

As with Atmospheric N deposition, free-living N inputs are added directly to the NH\ :sub:`4`\ :sup:`+` pool.

Nitrification and Denitrification Losses of Nitrogen
---------------------------------------------------------

Nitrification is an autotrophic process that converts less mobile ammonium ions into nitrate, that can more easily be lost from soil systems by leaching or denitrification. The process catalyzed by ammonia oxidizing archaea and bacteria that convert ammonium (NH\ :sub:`4`\ :sup:`+`) into nitrite, which is subsequently oxidized into nitrate (NO\ :sub:`3`\ :sup:`-`). Conditions favoring nitrification include high NH\ :sub:`4`\ :sup:`+` concentrations, well aerated soils, a neutral pH and warmer temperatures.

Under aerobic conditions in the soil oxygen is the preferred electron acceptor supporting the metabolism of heterotrophs, but anaerobic conditions favor the activity of soil heterotrophs which use nitrate as an electron acceptor (e.g. *Pseudomonas* and *Clostridium*) supporting respiration. This process, known as denitrification, results in the transformation of nitrate to gaseous N\ :sub:`2`, with smaller associated production of NO\ :sub:`x` and N\ :sub:`2`\ O. It is typically assumed that nitrogen fixation and denitrification were approximately balanced in the preindustrial biosphere ( :ref:`Galloway et al. 2004 <Gallowayetal2004>`). It is likely that denitrification can occur within anaerobic microsites within an otherwise aerobic soil environment, leading to large global denitrification fluxes even when fluxes per unit area are rather low (:ref:`Galloway et al. 2004 <Gallowayetal2004>`).

CLM includes a detailed representation of nitrification and denitrification based on the Century N model (:ref:`Parton et al. 1996 <Partonetal1996>`, :ref:`2001 <Partonetal2001>`; :ref:`del Grosso et al. 2000 <delGrossoetal2000>`). In this approach, nitrification of NH\ :sub:`4`\ :sup:`+` to NO\ :sub:`3`\ :sup:`-` is a function of temperature, moisture, and pH:

.. math::
   :label: 22.2)

   f_{nitr,p} =\left[NH_{4} \right]k_{nitr} f\left(T\right)f\left(H_{2} O\right)f\left(pH\right)

where :math:`{f}_{nitr,p}` is the potential nitrification rate (prior to competition for NH\ :sub:`4`\ :sup:`+` by plant uptake and N immobilization), :math:`{k}_{nitr}` is the maximum nitrification rate (10 % day\ :math:`\mathrm{-}`\ 1, (:ref:`Parton et al. 2001 <Partonetal2001>`), and *f(T)* and *f(H\)*\ :sub:`2`\ O) are rate modifiers for temperature and moisture content. CLM uses the same rate modifiers as are used in the decomposition routine. *f(pH)* is a rate modifier for pH; however, because CLM does not calculate pH, instead a fixed pH value of 6.5 is used in the pH function of :ref:`Parton et al. (1996) <Partonetal1996>`.

The potential denitrification rate is co-limited by NO\ :sup:`-3` concentration and C consumption rates, and occurs only in the anoxic fraction of soils:

.. math::
   :label: 22.3)

   f_{denitr,p} =\min \left(f(decomp),f\left(\left[NO_{3} ^{-} \right]\right)\right)frac_{anox}

where :math:`{f}_{denitr,p}` is the potential denitrification rate and *f(decomp)* and *f([NO*\ :sub:`3`\ :sup:`-` *])* are the carbon- and nitrate- limited denitrification rate functions, respectively, (:ref:`del Grosso et al. 2000 <delGrossoetal2000>`). Because the modified CLM includes explicit treatment of soil biogeochemical vertical profiles, including diffusion of the trace gases O\ :sub:`2` and CH\ :sub:`4` (:ref:`Riley et al. 2011a <Rileyetal2011a>`), the calculation of anoxic fraction :math:`{frac}_{anox}` uses this information following the anoxic microsite formulation of :ref:`Arah and Vinten (1995) <ArahVinten1995>`.

.. math::
   :label: 22.4)

   frac_{anox} =\exp \left(-aR_{\psi }^{-\alpha } V^{-\beta } C^{\gamma } \left[\theta +\chi \varepsilon \right]^{\delta } \right)

where *a*, :math:`\alpha`, :math:`\beta`, :math:`\gamma`, and :math:`\delta` are constants (equal to 1.5x10\ :sup:`-10`, 1.26, 0.6, 0.6, and 0.85, respectively), :math:`{R}_{\psi}` is the radius of a typical pore space at moisture content :math:`\psi`, *V* is the O\ :sub:`2` consumption rate, *C* is the O\ :sub:`2` concentration, :math:`\theta` is the water-filled pore space, :math:`\chi` is the ratio of diffusivity of oxygen in water to that in air, and :math:`\epsilon` is the air-filled pore space (:ref:`Arah and Vinten (1995) <ArahVinten1995>`). These parameters are all calculated separately at each layer to define a profile of anoxic porespace fraction in the soil.

The nitrification/denitrification models used here also predict fluxes of N\ :sub:`2`\ O via a "hole-in-the-pipe" approach (:ref:`Firestone and Davidson, 1989 <FirestoneDavidson1989>`). A constant fraction (6 * 10\ :math:`{}^{-4}`, :ref:`Li et al. 2000 <Lietal2000>`) of the nitrification flux is assumed to be N\ :sub:`2`\ O, while the fraction of denitrification going to N\ :sub:`2`\ O, \ :math:`{P}_{N2:N2O}`, is variable, following the Century (:ref:`del Grosso et al. 2000 <delGrossoetal2000>`) approach:

.. math::
   :label: 22.5)

   P_{N_{2} :N_{2} O} =\max \left(0.16k_{1} ,k_{1} \exp \left(-0.8P_{NO_{3} :CO_{2} } \right)\right)f_{WFPS}

where :math:`{P}_{NO3:CO2}` is the ratio of CO\ :sub:`2` production in a given soil layer to the NO\ :sub:`3`\ :sup:`-` concentration, :math:`{k}_{1}` is a function of :math:`{d}_{g}`, the gas diffusivity through the soil matrix:

.. math::
   :label: 22.6)

   k_{1} =\max \left(1.7,38.4-350*d_{g} \right)

and :math:`{f}_{WFPS}` is a function of the water filled pore space *WFPS:*

.. math::
   :label: 22.16)

   f_{WFPS} =\max \left(0.1,0.015\times WFPS-0.32\right)

Leaching Losses of Nitrogen
--------------------------------

Soil mineral nitrogen remaining after plant uptake, immobilization, and denitrification is subject to loss as a dissolved component of hydrologic outflow from the soil column (leaching). This leaching loss (:math:`{NF}_{leached}`, gN m\ :sup:`-2` s\ :sup:`-1`) depends on the concentration of dissolved mineral (inorganic) nitrogen in soil water solution (*DIN*, gN kgH\ :sub:`2`\ O), and the rate of hydrologic discharge from the soil column to streamflow (:math:`{Q}_{dis}`, kgH\ :sub:`2`\ O m\ :sup:`-2` s\ :sup:`-1`, section :numref:`Lateral Sub-surface Runoff`), as

.. math::
   :label: 22.17)

   NF_{leached} =DIN\cdot Q_{dis} .

*DIN* is calculated assuming that a constant fraction (*sf*, proportion) of the remaining soil mineral N pool is in soluble form, and that this entire fraction is dissolved in the total soil water. For the Century- based formulation in CLM5.0, the leaching acts only on the NO\ :sub:`3`\ :sup:`-` pool (which is assumed to be 100% soluble), while the NH\ :sub:`4`\ :sup:`+` pool is assumed to be 100% adsorbed onto mineral surfaces and unaffected by leaching. *DIN* is then given as

.. math::
   :label: 22.18)

   DIN=\frac{NS_{sminn} sf}{WS_{tot\_ soil} }

where :math:`{WS}_{tot\_soil}` (kgH\ :sub:`2`\ O m\ :sup:`-2`) is the total mass of soil water content integrated over the column. The total mineral nitrogen leaching flux is limited on each time step to not exceed the soluble fraction of :math:`{NS}_{sminn}`

.. math::
   :label: 22.19)

   NF_{leached} =\min \left(NF_{leached} ,\frac{NS_{sminn} sf}{\Delta t} \right).

Losses of Nitrogen Due to Fire
-----------------------------------

The final pathway for nitrogen loss is through combustion, also known as pyrodenitrification. Detailed equations are provided, together with the effects of fire on the carbon budget, in Chapter :numref:`rst_Fire`. It is assumed in CLM-CN that losses of N due to fire are restricted to vegetation and litter pools (including coarse woody debris). Loss rates of N are determined by the fraction of biomass lost to combustion, assuming that most of the nitrogen in the burned biomass is lost to the atmosphere (:ref:`Schlesinger, 1997 <Schlesinger1997>`; :ref:`Smith et al. 2005 <Smithetal2005>`). It is assumed that soil organic matter pools of carbon and nitrogen are not directly affected by fire (:ref:`Neff et al. 2005 <Neffetal2005>`).

