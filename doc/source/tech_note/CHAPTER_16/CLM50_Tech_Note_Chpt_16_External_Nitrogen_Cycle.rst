External Nitrogen Cycle
===========================

In addition to the relatively rapid cycling of nitrogen within the plant
– litter – soil organic matter system, CLM also represents several slow
processes which couple the internal nitrogen cycle to external sources
and sinks. Inputs of new mineral nitrogen are from atmospheric
deposition and biological nitrogen fixation. Losses of mineral nitrogen
are due to nitrification, denitrification, leaching, and losses in fire.
While the short-term dynamics of nitrogen limitation depend on the
behavior of the internal nitrogen cycle, establishment of total
ecosystem nitrogen stocks depends on the balance between sources and
sinks in the external nitrogen cycle.

CLM includes two separate treatments of the slow nitrogen cycle. The
first is the original CLM-CN formulations, which includes a single soil
mineral nitrogen pool to represent both
NO\ :sub:`3`\ :sup:`-` and NH\ :sub:`4`\ :sup:`+`, and with nitrogen losses via
denitrification calculated as a constant fraction of mineralization plus
a fixed first-order loss of unutilized mineral nitrogen at the end of
every timestep. The second is based on the Century N-gas model; this
includes separate NH\ :sub:`4`\ :sup:`+` and
NO\ :sub:`3`\ :sup:`-` pools, as well as
environmentally controlled nitrification and denitrification rates. Both
are described below.

Atmospheric Nitrogen Deposition
------------------------------------

CLM uses a single variable to represent the total deposition of mineral
nitrogen onto the land surface, combining wet and dry deposition of
NO\ :sub:`y` and NH\ :sub:`x` as a single flux
(:math:`{NF}_{ndep\_sminn}`, gN m\ :sup:`-2` s\ :sup:`-1`). This flux is intended to represent total reactive
nitrogen deposited to the land surface which originates from the
following natural and anthropogenic sources (Galloway et al. 2004):
formation of NO\ :sub:`x` during lightning,
NO\ :math:`{}_{x }`\ and NH\ :sub:`3` emission from wildfire,
NO\ :sub:`x` emission from natural soils, NH\ :sub:`3`
emission from natural soils, vegetation, and wild animals,
NO\ :sub:`x` and NH\ :sub:`3` emission during fossil fuel
combustion (both thermal and fuel NO\ :sub:`x` production),
NO\ :sub:`x` and NH\ :sub:`3` emission from other industrial
processes, NO\ :sub:`x` and NH\ :sub:`3` emission from fire
associated with deforestation, NO\ :sub:`x` and NH\ :sub:`3`
emission from agricultural burning, NO\ :sub:`x` emission from
agricultural soils, NH\ :sub:`3` emission from agricultural crops,
NH\ :sub:`3` emission from agricultural animal waste, and
NH\ :sub:`3` emission from human waste and waste water. The
deposition flux is provided as a spatially and (potentially) temporally
varying dataset (see section 2.2.3 for a description of the default
input dataset).

In the CLM-CN mineral N pool model, the nitrogen deposition flux is
assumed to enter the soil mineral nitrogen pool
(:math:`{NS}_{sminn}`) directly; while in the Century-based model,
all of the nitrogen is assumed to enter the
NH\ :sub:`4`\ :sup:`+` pool. Real pathways for wet and dry
nitrogen deposition can be more complex than currently represented in
the CLM-CN, including release from melting snowpack and direct foliar
uptake of deposited NO\ :sub:`y` (e.g. Tye et al. 2005; Vallano
and Sparks, 2007).

Biological Nitrogen Fixation
---------------------------------

The fixation of new reactive nitrogen from atmospheric N\ :sub:`2`
by soil microorganisms is an important component of both preindustrial
and modern-day nitrogen budgets, but a mechanistic understanding of
global-scale controls on biological nitrogen fixation (BNF) is still
only poorly developed (Cleveland et al. 1999; Galloway et al. 2004).
Cleveland et al. (1999) suggested empirical relationships that predict
BNF as a function of either evapotranspiration rate or net primary
productivity for natural vegetation. CLM assumes that BNF is a function
of net primary production (:math:`{CF}_{ann\_NPP}`, gC m\ :sup:`-2` y\ :sup:`-1`). The rationale for choosing net
primary production over evapotranspiration as the predictor is that the
two are well-correlated (Parton et al. 1993; Running et al. 1989), and
the use of primary production also introduces a known dependence of BNF
on the carbon supply to nitrogen fixing microorganisms (Cleveland et al.
1999). The expression used is:

.. math::
   :label: ZEqnNum802819: 

   NF_{nfix,sminn} ={1.8\left(1-\exp \left(-0.003{\kern 1pt} CF_{ann\_ NPP} \right)\right)\mathord{\left/ {\vphantom {1.8\left(1-\exp \left(-0.003{\kern 1pt} CF_{ann\_ NPP} \right)\right) \left(86400\cdot 365\right)}} \right. \kern-\nulldelimiterspace} \left(86400\cdot 365\right)}

where :math:`{NF}_{nfix,sminn}` (gN m\ :sup:`-2` s\ :sup:`-1`) is the rate of BNF. Eq. is plotted over a range of
annual NPP in Figure 16.1.

Figure 16.1. Biological nitrogen fixation as a function of annual net
primary production.

.. image:: image1.png

Because of the empirical nature of this NPP-BNF relationship, the
timescale for calculating NPP and thus BNF is unconstrained. Using
annual NPP, as in CLM4.0, introduces an error at high latitudes because
the aseasonal BNF inputs mean that much of the nitrogen is added when
the vegetation is dormant and may be lost before it is ever taken up by
vegetation. Thus an option was added to CLM to allow for an exponential
relaxation (with default e-folding time of 10 days) calculation of NPP,
and BNF calculated from that using equation 16.1.

As with Atmospheric N deposition, BNF N inputs are added directly to the
mineral N pools. In the CLM-CN N model, this is the single mineral N
pool; in the Century-based model, this is the
NH\ :sub:`4`\ :sup:`+` pool.

Nitrification and Denitrification Losses of Nitrogen
---------------------------------------------------------

In order to better understand the structural uncertainty in
biogeochemical responses to climate change, CLM includes two alternate
representations of the mineral N transformations and losses that define
the slow N cycle. Each of these is described below.

16.3.1 CLM-CN formulation
^^^^^^^^^^^^^^^^^^^^^^^^^

Under aerobic conditions in the soil oxygen is the preferred electron
acceptor supporting the metabolism of heterotrophs, but anaerobic
conditions favor the activity of soil heterotrophs which use nitrate as
an electron acceptor (e.g. *Pseudomonas* and *Clostridium*) supporting
respiration. This process, known as denitrification, results in the
transformation of nitrate to gaseous N\ :sub:`2`, with smaller
associated production of NO\ :sub:`x` and N\ :sub:`2`\ O. It
is typically assumed that nitrogen fixation and denitrification were
approximately balanced in the preindustrial biosphere (Galloway et al.
2004). It is likely that denitrification can occur within anaerobic
microsites within an otherwise aerobic soil environment, leading to
large global denitrification fluxes even when fluxes per unit area are
rather low (Galloway et al. 2004).

Because the vertical distribution of soil organic matter is not resolved
explicitly in CLM-CN, a simple denitrification parameterization is used
that treats denitrification as a constant fraction of gross nitrogen
mineralization. At each step in the decomposition cascade, if the
transformation from an upstream to a downstream pool is predicted to
mineralize (as opposed to immobilize) nitrogen, then a constant fraction
of the nitrogen mineralization flux is assumed to be lost via
denitrification. Due to large uncertainties in the mechanistic
understanding of the environmental controls on denitrification, no
modifications to the denitrification fraction are made for different
soil moisture conditions. This is identified as a high-priority area for
future model development.

Denitrification fluxes associated with gross mineralization in the
decomposition cascade are calculated as follows:

.. math::
   :label: 16.2) 

   NF_{denit,Lit1\to SOM1} =\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{pot\_ min,Lit1\to SOM1} >0} \\ {-NF_{pot\_ min,Lit1\to SOM1} {\kern 1pt} f_{denit} \qquad \qquad {\rm for\; }NF_{pot\_ min,Lit1\to SOM1} \le 0} \end{array}\right.

.. math::
   :label: 16.3) 

   NF_{denit,Lit2\to SOM2} =\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{pot\_ min,Lit2\to SOM2} >0} \\ {-NF_{pot\_ min,Lit2\to SOM2} {\kern 1pt} f_{denit} \qquad \qquad {\rm for\; }NF_{pot\_ min,Lit2\to SOM2} \le 0} \end{array}\right.

.. math::
   :label: 16.4) 

   NF_{denit,Lit3\to SOM3} =\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{pot\_ min,Lit3\to SOM3} >0} \\ {-NF_{pot\_ min,Lit3\to SOM3} {\kern 1pt} f_{denit} \qquad \qquad {\rm for\; }NF_{pot\_ min,Lit3\to SOM3} \le 0} \end{array}\right.

.. math::
   :label: 16.5) 

   NF_{denit,SOM1\to SOM2} =\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{pot\_ min,SOM1\to SOM2} >0} \\ {-NF_{pot\_ min,SOM1\to SOM2} {\kern 1pt} f_{denit} \qquad {\rm for\; }NF_{pot\_ min,SOM1\to SOM2} \le 0} \end{array}\right.

.. math::
   :label: 16.6) 

   NF_{denit,SOM2\to SOM3} =\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{pot\_ min,SOM2\to SOM3} >0} \\ {-NF_{pot\_ min,SOM2\to SOM3} {\kern 1pt} f_{denit} \qquad {\rm for\; }NF_{pot\_ min,SOM2\to SOM3} \le 0} \end{array}\right.

.. math::
   :label: 16.7) 

   NF_{denit,SOM3\to SOM4} =\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{pot\_ min,SOM3\to SOM4} >0} \\ {-NF_{pot\_ min,SOM3\to SOM4} {\kern 1pt} f_{denit} \qquad {\rm for\; }NF_{pot\_ min,SOM3\to SOM4} \le 0} \end{array}\right.

.. math::
   :label: 16.8) 

   NF_{denit,SOM4} =-NF_{pot\_ min,SOM4}

where :math:`{f}_{denit} = 0.01` is the constant denitrification
fraction of gross mineralization, and the denitrification fluxes are
assumed to be leaving the soil mineral nitrogen pool
(:math:`{NS}_{sminn}`) and entering the atmosphere. The speciation
of gaseous nitrogen fluxes entering the atmosphere (e.g.
N\ :sub:`2` vs. NO\ :sub:`x` or N\ :sub:`2`\ O) is not
specified. Providing an explicit speciation of these nitrogen losses is
another high-priority area for future model development.

The model includes one other denitrification pathway, intended to
represent the observed losses of mineral nitrogen in systems
experiencing nitrogen saturation. One reason this mechanism has been
included is in anticipation of an agricultural fertilization flux,
provided either through a prescribed dataset or through a prognostic
agricultural management routine. The model does not currently include an
explicit representation of the fertilization flux, but when it is
introduced, it will be necessary to account for the substantial
denitrification losses associated with high nitrate concentrations in
some heavily fertilized agricultural soils. Nitrogen saturation can also
occur in natural vegetation systems, especially under conditions of high
atmospheric nitrogen deposition, and so this mechanism plays a useful
role even prior to the introduction within the model of agricultural
fertilization.

For the purpose of this calculation, nitrogen saturation is evaluated on
each timestep, by comparing the total demand for new mineral nitrogen
from plants and immobilization with the available soil mineral nitrogen
pool. The denitrification of excess soil mineral nitrogen is non-zero
whenever the supply of mineral nitrogen exceeds the demand:

.. math::
   :label: 16.9) 

   NF_{sminn,denit} =\left\{\begin{array}{l} {\left(\frac{NS_{sminn} }{\Delta t} \right)-NF_{total\_ demand} f_{dnx} \qquad {\rm for\; }NF_{total\_ demand} \Delta t<NS_{sminn} } \\ {0\qquad \qquad \qquad \qquad {\rm for\; }NF_{total\_ demand} \Delta t\ge NS_{sminn} } \end{array}\right.

where :math:`{f}_{dnx}` (unitless) is the fraction of excess soil
mineral nitrogen subject to denitrification on each timestep. This
fraction is parameterized such that 50% of any excess soil mineral
nitrogen would be lost to denitrification per day:

.. math::
   :label: 16.10) 

   f_{dnx} =0.5\frac{\Delta t}{86400}

16.3.2 Century-based formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CLM includes a detailed representation of nitrification and
denitrification based on the Century N model (Parton et al. 1996, 2001;
del Grosso et al. 2000). In this approach, nitrification of
NH\ :sub:`4`\ :sup:`+` to NO\ :sub:`3`\ :sup:`-``
is a function of temperature, moisture, and pH:

.. math::
   :label: 16.11) 

   f_{nitr,p} =\left[NH_{4} \right]k_{nitr} f\left(T\right)f\left(H_{2} O\right)f\left(pH\right)

where :math:`{f}_{nitr,p}` is the potential nitrification rate
(prior to competition for NH\ :sub:`4`\ :sup:`+` by plant
uptake and N immobilization), :math:`{k}_{nitr}` is the maximum
nitrification rate (10 % day\ :math:`\mathrm{-}`\ 1, (Parton et al.
2001)), and *f(T)* and *f(H\)*\ :sub:`2`\ O) are rate modifiers for
temperature and moisture content. CLM uses the same rate modifiers as
are used in the decomposition routine. *f(pH)* is a rate modifier for
pH; however, because CLM does not calculate pH, instead a fixed pH value
of 6.5 is used in the pH function of Parton et al. (1996).

The potential denitrification rate is co-limited by
NO\ :sup:`-3` concentration and C consumption rates, and occurs only in the anoxic fraction of soils:

.. math::
   :label: 16.12) 

   f_{denitr,p} =\min \left(f(decomp),f\left(\left[NO_{3} ^{-} \right]\right)\right)frac_{anox}

where :math:`{f}_{denitr,p}` is the potential denitrification rate
and *f(decomp)* and *f([NO*\ :sub:`3`\ :sup:`-` *])*
are the carbon- and nitrate- limited denitrification rate functions,
respectively, (del Grosso et al. 2000). Because the modified CLM
includes explicit treatment of soil biogeochemical vertical profiles,
including diffusion of the trace gases O\ :sub:`2` and
CH\ :sub:`4` (Riley et al. 2011a), the calculation of anoxic
fraction  :math:`{frac}_{anox}` uses this information following the
anoxic microsite formulation of Arah and Vinten (1995):

.. math::
   :label: 16.13) 

   frac_{anox} =\exp \left(-aR_{\psi }^{-\alpha } V^{-\beta } C^{\gamma } \left[\theta +\chi \varepsilon \right]^{\delta } \right)

where *a*, :math:`\alpha`, :math:`\beta`, :math:`\gamma`, and :math:`\delta` are constants (equal to
1.5x10\ :sup:`-10`, 1.26, 0.6, 0.6, and 0.85, respectively), :math:`{R}_{\psi}` is the
radius of a typical pore space at moisture content :math:`\psi`, *V*
is the O\ :sub:`2` consumption rate, *C* is the O\ :sub:`2`
concentration, :math:`\theta` is the water-filled pore space,
:math:`\chi` is the ratio of diffusivity of oxygen in water to that in
air, and :math:`\epsilon` is the air-filled pore space (Arah and
Vinten, 1995). These parameters are all calculated separately at each
layer to define a profile of anoxic porespace fraction in the soil.

The nitrification/denitrification models used here also predict fluxes
of N\ :sub:`2`\ O via a “hole-in-the-pipe” approach (Firestone and
Davidson, 1989). A constant fraction (6 \* 10\ :math:`{}^{-4}`, Li et
al. 2000) of the nitrification flux is assumed to be
N\ :sub:`2`\ O, while the fraction of denitrification going to
N\ :sub:`2`\ O, P\ :math:`{P}_{N2:N2O}`, is variable, following
the Century (del Grosso et al. 2000) approach:

.. math::
   :label: 16.14) 

   P_{N_{2} :N_{2} O} =\max \left(0.16k_{1} ,k_{1} \exp \left(-0.8P_{NO_{3} :CO_{2} } \right)\right)f_{WFPS}

where :math:`{P}_{NO3:CO2}` is the ratio of CO\ :sub:`2`
production in a given soil layer to the
NO\ :sub:`3`\ :sup:`-`` concentration, :math:`{k}_{1}` is
a function of :math:`{d}_{g}`, the gas diffusivity through the soil
matrix:

.. math::
   :label: 16.15) 

   k_{1} =\max \left(1.7,38.4-350*d_{g} \right)

and :math:`{f}_{WFPS}` is a function of the water filled pore space *WFPS:*

.. math::
   :label: 16.16) 

   f_{WFPS} =\max \left(0.1,0.015\times WFPS-0.32\right)

Leaching Losses of Nitrogen
--------------------------------

Soil mineral nitrogen remaining after plant uptake, immobilization, and
denitrification is subject to loss as a dissolved component of
hydrologic outflow from the soil column (leaching). This leaching loss
(:math:`{NF}_{leached}`, gN m\ :sup:`-2` s\ :sup:`-1`)
depends on the concentration of dissolved mineral (inorganic) nitrogen
in soil water solution (*DIN*, gN kgH\ :sub:`2`\ O), and the rate
of hydrologic discharge from the soil column to streamflow
(:math:`{Q}_{dis}`, kgH\ :sub:`2`\ O m\ :sup:`-2`
s\ :sup:`-1`, section 7.6), as

.. math::
   :label: 16.17) 

   NF_{leached} =DIN\cdot Q_{dis} .

*DIN* is calculated assuming that a constant fraction (*sf*, proportion)
of the remaining soil mineral N pool is in soluble form, and that this
entire fraction is dissolved in the total soil water. For the CLM-CN
soil model, it is further assumed that *sf* = 0.1, representing an
estimated 10% of the total :math:`{NS}_{sminn}` pool as soluble
nitrate, with the remaining 90% as less soluble ammonia; for the
Century-based formulation, the leaching acts only on the
NO\ :sub:`3`\ :sup:`-`` pool (which is assumed to be 100%
soluble), while the NH\ :sub:`4`\ :sup:`+` pool is assumed
to be 100% adsorbed onto mineral surfaces and unaffected by leaching.
*DIN* is then given as

.. math::
   :label: 16.18) 

   DIN=\frac{NS_{sminn} sf}{WS_{tot\_ soil} }

where :math:`{WS}_{tot\_soil}` (kgH:sub:`2`\ O m\ :sup:`-2`) is the total mass of soil water content integrated
over the column. The total mineral nitrogen leaching flux is limited on
each time step to not exceed the soluble fraction of :math:`{NS}_{sminn}`

.. math::
   :label: 16.19) 

   NF_{leached} =\min \left(NF_{leached} ,\frac{NS_{sminn} sf}{\Delta t} \right).

The CLM-CN parameterization of the soluble fraction is poorly
constrained by observations. Fraction of total soil mineral N pool
present as nitrate will vary spatially and temporally, depending on
oxygen status of soils and rates of nitrification. A calibration of this
parameterization against observations of dissolved nitrate in headwater
streams might be an effective method for imposing better observational
constraints at broad spatial scales.

Losses of Nitrogen Due to Fire
-----------------------------------

The final pathway for nitrogen loss is through combustion, also known as
pyrodenitrification. Detailed equations are provided, together with the
effects of fire on the carbon budget, in Chapter 18. It is assumed in
CLM-CN that losses of N due to fire are restricted to vegetation and
litter pools (including coarse woody debris). Loss rates of N are
determined by the fraction of biomass lost to combustion, assuming that
most of the nitrogen in the burned biomass is lost to the atmosphere
(Schlesinger, 1997; Smith et al. 2005). It is assumed that soil organic
matter pools of carbon and nitrogen are not directly affected by fire
(Neff et al. 2005).

