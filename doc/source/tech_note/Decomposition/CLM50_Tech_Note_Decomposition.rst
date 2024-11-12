.. _rst_Decomposition:

Decomposition
=================

Decomposition of fresh litter material into progressively more recalcitrant forms of soil organic matter is represented in CLM is defined as a cascade of :math:`{k}_{tras}` transformations between :math:`{m}_{pool}` decomposing coarse woody debris (CWD), litter, and soil organic matter (SOM) pools, each defined at :math:`{n}_{lev}` vertical levels. CLM allows the user to define, at compile time, between 2 contrasting hypotheses of decomposition as embodied by two separate decomposition submodels: the CLM-CN pool structure used in CLM4.0, or a second pool structure, characterized by slower decomposition rates, based on the fCentury model (Parton et al 1988). In addition, the user can choose, at compile time, whether to allow :math:`{n}_{lev}` to equal 1, as in CLM4.0, or to equal the number of soil levels used for the soil hydrological and thermal calculations (see Section :numref:`Soil Layers` for soil layering).

.. _Figure Schematic of decomposition model in CLM:

.. figure:: CLM4_vertsoil_soilstruct_drawing.png

 Schematic of decomposition model in CLM.

Model is structured to allow different representations of the soil C and N decomposition cascade, as well as a vertically-explicit treatment of soil biogeochemistry.

For the single-level model structure, the fundamental equation for carbon balance of the decomposing pools is:

.. math::
   :label: 21.1)

   \frac{\partial C_{i} }{\partial t} =R_{i} +\sum _{j\ne i}\left(i-r_{j} \right)T_{ji} k_{j} C_{j} -k_{i} C_{i}

where :math:`{C}_{i}` is the carbon content of pool *i*, :math:`{R}_{i}` are the carbon inputs from plant tissues directly to pool *i* (only non-zero for CWD and litter pools), :math:`{k}_{i}` is the decay constant of pool *i*; :math:`{T}_{ji}` is the fraction of carbon directed from pool *j* to pool *i* with fraction :math:`{r}_{j}` lost as a respiration flux along the way.

Adding the vertical dimension to the decomposing pools changes the balance equation to the following:

.. math::
   :label: 21.2)

   \begin{array}{l} {\frac{\partial C_{i} (z)}{\partial t} =R_{i} (z)+\sum _{i\ne j}\left(1-r_{j} \right)T_{ji} k_{j} (z)C_{j} (z) -k_{i} (z)C_{i} (z)} \\ {+\frac{\partial }{\partial z} \left(D(z)\frac{\partial C_{i} }{\partial z} \right)+\frac{\partial }{\partial z} \left(A(z)C_{i} \right)} \end{array}

where :math:`{C}_{i}`\ (z) is now defined at each model level, and in volumetric (gC m\ :sup:`-3`) rather than areal (gC m\ :sup:`-2`) units, along with :math:`{R}_{i}`\ (z) and :math:`{k}_{j}`\ (z). In addition, vertical transport is handled by the last two terms, for diffusive and advective transport. In the base model, advective transport is set to zero, leaving only a diffusive flux with diffusivity *D(z)* defined for all decomposing carbon and nitrogen pools. Further discussion of the vertical distribution of carbon inputs :math:`{R}_{i}`\ (z), vertical turnover times :math:`{k}_{j}`\ (z), and vertical transport *D(z)* is below Discussion of the vertical model and analysis of both decomposition structures is in :ref:`Koven et al. (2013) <Kovenetal2013>`.

.. _Figure Pool structure:

.. figure:: soil_C_pools_CN_century.png

 Pool structure, transitions, respired fractions (numbers at
 end of arrows), and turnover times (numbers in boxes) for the 2
 alternate soil decomposition models included in CLM.

CLM-CN Pool Structure, Rate Constants and Parameters
---------------------------------------------------------

The CLM-CN structure in CLM45 uses three state variables for fresh litter and four state variables for soil organic matter (SOM). The masses of carbon and nitrogen in the live microbial community are not modeled explicitly, but the activity of these organisms is represented by decomposition fluxes transferring mass between the litter and SOM pools, and heterotrophic respiration losses associated with these transformations. The litter and SOM pools in CLM-CN are arranged as a converging cascade (Figure 15.2), derived directly from the implementation in Biome-BGC v4.1.2 (Thornton et al. 2002; Thornton and Rosenbloom, 2005).

Model parameters are estimated based on a synthesis of microcosm decomposition studies using radio-labeled substrates (Degens and Sparling, 1996; Ladd et al. 1992; Martin et al. 1980; Mary et al. 1993 Saggar et al. 1994; Sørensen, 1981; van Veen et al. 1984). Multiple exponential models are fitted to data from the microcosm studies to estimate exponential decay rates and respiration fractions (Thornton, 1998). The microcosm experiments used for parameterization were all conducted at constant temperature and under moist conditions with relatively high mineral nitrogen concentrations, and so the resulting rate constants are assumed not limited by the availability of water or mineral nitrogen. :numref:`Table Decomposition rate constants` lists the base decomposition rates for each litter and SOM pool, as well as a base rate for physical fragmentation for the coarse woody debris pool (CWD).

.. _Table Decomposition rate constants:

.. table:: Decomposition rate constants for litter and SOM pools, C:N ratios, and acceleration parameters for the CLM-CN decomposition pool structure.

 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 |                          | Biome-BGC                                      | CLM-CN                                        |               |                                         |
 +==========================+================================================+===============================================+===============+=========================================+
 |                          | :math:`{k}_{disc1}`\ (d\ :sup:`-1`)            | :math:`{k}_{disc2}` (hr\ :sup:`-1`)           | *C:N ratio*   | Acceleration term (:math:`{a}_{i}`)     |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{Lit1}`       | 0.7                                            | 0.04892                                       | -             | 1                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{Lit2}`       | 0.07                                           | 0.00302                                       | -             | 1                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{Lit3}`       | 0.014                                          | 0.00059                                       | -             | 1                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{SOM1}`       | 0.07                                           | 0.00302                                       | 12            | 1                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{SOM2}`       | 0.014                                          | 0.00059                                       | 12            | 1                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{SOM3}`       | 0.0014                                         | 0.00006                                       | 10            | 5                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{SOM4}`       | 0.0001                                         | 0.000004                                      | 10            | 70                                      |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+
 | :math:`{k}_{CWD}`        | 0.001                                          | 0.00004                                       | -             | 1                                       |
 +--------------------------+------------------------------------------------+-----------------------------------------------+---------------+-----------------------------------------+

The first column of :numref:`Table Decomposition rate constants` gives the rates as used for the Biome-BGC model, which uses a discrete-time model with a daily timestep. The second column of :numref:`Table Decomposition rate constants` shows the rates transformed for a one-hour discrete timestep typical of CLM-CN. The transformation is based on the conversion of the initial discrete-time value (:math:`{k}_{disc1}` first to a continuous time value (:math:`{k}_{cont}`), then to the new discrete-time value with a different timestep (:math:`{k}_{disc2}`), following Olson (1963):

.. math::
   :label: ZEqnNum608251

   k_{cont} =-\log \left(1-k_{disc1} \right)

.. math::
   :label: ZEqnNum772630

   k_{disc2} =1-\exp \left(-k_{cont} \frac{\Delta t_{2} }{\Delta t_{1} } \right)

where :math:`\Delta`\ :math:`{t}_{1}` (s) and :math:`\Delta`\ t\ :sub:`2` (s) are the time steps of the initial and new discrete-time models, respectively.

Respiration fractions are parameterized for decomposition fluxes out of each litter and SOM pool. The respiration fraction (*rf*, unitless) is the fraction of the decomposition carbon flux leaving one of the litter or SOM pools that is released as CO\ :sub:`2` due to heterotrophic respiration. Respiration fractions and exponential decay rates are estimated simultaneously from the results of microcosm decomposition experiments (Thornton, 1998). The same values are used in CLM-CN and Biome-BGC (:numref:`Table Respiration fractions for litter and SOM pools`).

.. _Table Respiration fractions for litter and SOM pools:

.. table:: Respiration fractions for litter and SOM pools

 +---------------------------+-----------------------+
 | Pool                      | *rf*                  |
 +===========================+=======================+
 |  :math:`{rf}_{Lit1}`      | 0.39                  |
 +---------------------------+-----------------------+
 |  :math:`{rf}_{Lit2}`      | 0.55                  |
 +---------------------------+-----------------------+
 |  :math:`{rf}_{Lit3}`      | 0.29                  |
 +---------------------------+-----------------------+
 |  :math:`{rf}_{SOM1}`      | 0.28                  |
 +---------------------------+-----------------------+
 |  :math:`{rf}_{SOM2}`      | 0.46                  |
 +---------------------------+-----------------------+
 |  :math:`{rf}_{SOM3}`      | 0.55                  |
 +---------------------------+-----------------------+
 |  :math:`{rf}_{SOM4}`      |  :math:`{1.0}^{a}`    |
 +---------------------------+-----------------------+

:sup:`a`:math:`{}^{a}` The respiration fraction for pool SOM4 is 1.0 by definition: since there is no pool downstream of SOM4, the entire carbon flux leaving this pool is assumed to be respired as CO\ :sub:`2`.

Century-based Pool Structure, Rate Constants and Parameters
----------------------------------------------------------------

The Century-based decomposition cascade is, like CLM-CN, a first-order decay model; the two structures differ in the number of pools, the connections between those pools, the turnover times of the pools, and the respired fraction during each transition (Figure 15.2). The turnover times are different for the Century-based pool structure, following those described in Parton et al. (1988) (:numref:`Table Turnover times`).

.. _Table Turnover times:

.. table:: Turnover times, C:N ratios, and acceleration parameters for the Century-based decomposition cascade.

 +------------+------------------------+-------------+-------------------------------------------+
 |            | Turnover time (year)   | C:N ratio   | Acceleration term (:math:`{a}_{i}`)       |
 +============+========================+=============+===========================================+
 | CWD        | 4.1                    | -           | 1                                         |
 +------------+------------------------+-------------+-------------------------------------------+
 | Litter 1   | 0.066                  | -           | 1                                         |
 +------------+------------------------+-------------+-------------------------------------------+
 | Litter 2   | 0.25                   | -           | 1                                         |
 +------------+------------------------+-------------+-------------------------------------------+
 | Litter 3   | 0.25                   | -           | 1                                         |
 +------------+------------------------+-------------+-------------------------------------------+
 | SOM 1      | 0.17                   | 8           | 1                                         |
 +------------+------------------------+-------------+-------------------------------------------+
 | SOM 2      | 6.1                    | 11          | 15                                        |
 +------------+------------------------+-------------+-------------------------------------------+
 | SOM 3      | 270                    | 11          | 675                                       |
 +------------+------------------------+-------------+-------------------------------------------+

Likewise, values for the respiration fraction of Century-based structure are in :numref:`Table Respiration fractions for Century-based structure`.

.. _Table Respiration fractions for Century-based structure:

.. table::  Respiration fractions for litter and SOM pools for Century-based structure

 +---------------------------+----------+
 | Pool                      | *rf*     |
 +===========================+==========+
 |  :math:`{rf}_{Lit1}`      | 0.55     |
 +---------------------------+----------+
 |  :math:`{rf}_{Lit2}`      | 0.5      |
 +---------------------------+----------+
 |  :math:`{rf}_{Lit3}`      | 0.5      |
 +---------------------------+----------+
 |  :math:`{rf}_{SOM1}`      | f(txt)   |
 +---------------------------+----------+
 |  :math:`{rf}_{SOM2}`      | 0.55     |
 +---------------------------+----------+
 |  :math:`{rf}_{SOM3}`      | 0.55     |
 +---------------------------+----------+

Environmental modifiers on decomposition rate
--------------------------------------------------

These base rates are modified on each timestep by functions of the current soil environment. For the single-level model, there are two rate modifiers, temperature (:math:`{r}_{tsoil}`, unitless) and moisture (:math:`{r}_{water}`, unitless), both of which are calculated using the average environmental conditions of the top five model levels (top 29 cm of soil column). For the vertically-resolved model, two additional environmental modifiers are calculated beyond the temperature and moisture limitations: an oxygen scalar (:math:`{r}_{oxygen}`, unitless), and a depth scalar (:math:`{r}_{depth}`, unitless).

The Temperature scalar :math:`{r}_{tsoil}` is calculated in CLM using a :math:`{Q}_{10}` approach, with :math:`{Q}_{10} = 1.5`.

.. math::
   :label: 21.5)

   r_{tsoil} =Q_{10} ^{\left(\frac{T_{soil,\, j} -T_{ref} }{10} \right)}

where *j* is the soil layer index, :math:`{T}_{soil,j}` (K) is the temperature of soil level *j*. The reference temperature :math:`{T}_{ref}` = 25C.

The rate scalar for soil water potential (:math:`{r}_{water}`, unitless) is calculated using a relationship from Andrén and Paustian (1987) and supported by additional data in Orchard and Cook (1983):

.. math::
   :label: 21.6)

   r_{water} =\sum _{j=1}^{5}\left\{\begin{array}{l} {0\qquad {\rm for\; }\Psi _{j} <\Psi _{\min } } \\ {\frac{\log \left({\Psi _{\min } \mathord{\left/ {\vphantom {\Psi _{\min }  \Psi _{j} }} \right.} \Psi _{j} } \right)}{\log \left({\Psi _{\min } \mathord{\left/ {\vphantom {\Psi _{\min }  \Psi _{\max } }} \right.} \Psi _{\max } } \right)} w_{soil,\, j} \qquad {\rm for\; }\Psi _{\min } \le \Psi _{j} \le \Psi _{\max } } \\ {1\qquad {\rm for\; }\Psi _{j} >\Psi _{\max } \qquad \qquad } \end{array}\right\}

where :math:`{\Psi}_{j}` is the soil water potential in layer *j*, :math:`{\Psi}_{min}` is a lower limit for soil water potential control on decomposition rate (in CLM5, this was changed from a default value of -10 MPa used in CLM4.5 and earlier to a default value of -2.5 MPa). :math:`{\Psi}_{max,j}` (MPa) is the soil moisture at which decomposition proceeds at a moisture-unlimited rate. The default value of :math:`{\Psi}_{max,j}` for CLM5 is updated from a saturated value used in CLM4.5 and earlier, to a value nominally at field capacity, with a value of -0.002 MPa For frozen soils, the bulk of the rapid dropoff in decomposition with decreasing temperature is due to the moisture limitation, since matric potential is limited by temperature in the supercooled water formulation of Niu and Yang (2006),

.. math::
   :label: 21.8)

   \psi \left(T\right)=-\frac{L_{f} \left(T-T_{f} \right)}{10^{3} T}

An additional frozen decomposition limitation can be specified using a ‘frozen Q\ :sub:`10`' following :ref:`Koven et al. (2011) <Kovenetal2011>`, however the default value of this is the same as the unfrozen Q\ :sub:`10` value, and therefore the basic hypothesis is that frozen respiration is limited by liquid water availability, and can be modeled following the same approach as thawed but dry soils.

An additional rate scalar, :math:`{r}_{oxygen}` is enabled when the CH\ :sub:`4` submodel is used (set equal to 1 for the single layer model or when the CH\ :sub:`4` submodel is disabled). This limits decomposition when there is insufficient molecular oxygen to satisfy stoichiometric demand (1 mol O\ :sub:`2` consumed per mol CO\ :sub:`2` produced) from heterotrophic decomposers, and supply from diffusion through soil layers (unsaturated and saturated) or aerenchyma (Chapter 19). A minimum value of :math:`{r}_{oxygen}` is set at 0.2, with the assumption that oxygen within organic tissues can supply the necessary stoichiometric demand at this rate. This value lies between estimates of 0.025–0.1 (Frolking et al. 2001), and 0.35 (Wania et al. 2009); the large range of these estimates poses a large unresolved uncertainty.

Lastly, a possible explicit depth dependence, :math:`{r}_{depth}`, (set equal to 1 for the single layer model) can be applied to soil C decomposition rates to account for processes other than temperature, moisture, and anoxia that can limit decomposition. This depth dependence of decomposition was shown by Jenkinson and Coleman (2008) to be an important term in fitting total C and 14C profiles, and implies that unresolved processes, such as priming effects, microscale anoxia, soil mineral surface and/or aggregate stabilization may be important in controlling the fate of carbon at depth :ref:`Koven et al. (2013) <Kovenetal2013>`. CLM includes these unresolved depth controls via an exponential decrease in the soil turnover time with depth:

.. math::
   :label: 21.9)

   r_{depth} =\exp \left(-\frac{z}{z_{\tau } } \right)

where :math:`{z}_{\tau}` is the e-folding depth for decomposition. For CLM4.5, the default value of this was 0.5m. For CLM5, this has been changed to a default value of 10m, which effectively means that intrinsic decomposition rates may proceed as quickly at depth as at the surface.

The combined decomposition rate scalar (:math:`{r}_{total}`,unitless) is:

.. math::
   :label: 21.10)

   r_{total} =r_{tsoil} r_{water} r_{oxygen} r_{depth} .

.. _decomp_mgmt_modifiers:

Management modifiers on decomposition rate
--------------------------------------------------

Tillage of cropland soil is represented as an additional rate scalar that depends on tillage intensity (default off), soil pool, and time since planting :ref:`(Graham et al., 2021) <Grahametal2021>`. The tillage enhancement is strongest in the first 14 days after planting (idpp < 15), weaker in the next 30 days (15 ≤ idpp < 45), weaker still in the next 30 days (45 ≤ idpp < 75), and nonexistent after that (idpp ≥ 75).

.. list-table:: Tillage decomposition rate scalars. Values in each cell represent enhancement in different periods of days past planting: [0, 14], [15, 44], [45, 74].
   :header-rows: 1

   * - \
     - low
     - high
   * - Litter 2 (cel_lit)
     - 1.5, 1.5, 1.1
     - 1.8, 1.5, 1.1
   * - Litter 3 (lig_lit)
     - 1.5, 1.5, 1.1
     - 1.8, 1.5, 1.1
   * - SOM 1 (act_som)
     - 1.0, 1.0, 1.0
     - 1.2, 1.0, 1.0
   * - SOM 2 (slo_som)
     - 3.0, 1.6, 1.3
     - 4.8, 3.5, 2.5
   * - SOM 3 (pas_som)
     - 3.0, 1.6, 1.3
     - 4.8, 3.5, 2.5

N-limitation of Decomposition Fluxes
-----------------------------------------

Decomposition rates can also be limited by the availability of mineral nitrogen, but calculation of this limitation depends on first estimating the potential rates of decomposition, assuming an unlimited mineral nitrogen supply. The general case is described here first, referring to a generic decomposition flux from an "upstream" pool (*u*) to a "downstream" pool (*d*), with an intervening loss due to respiration The potential carbon flux out of the upstream pool (:math:`{CF}_{pot,u}`, gC m\ :sup:`-2` s\ :sup:`-1`) is:

.. math::
   :label: 21.11)

   CF_{pot,\, u} =CS_{u} k_{u}

where :math:`{CS}_{u}` (gC m\ :sup:`-2`) is the initial mass in the upstream pool and :math:`{k}_{u}` is the decay rate constant (s\ :sup:`-1`) for the upstream pool, adjusted for temperature and moisture conditions. Depending on the C:N ratios of the upstream and downstream pools and the amount of carbon lost in the transformation due to respiration (the respiration fraction), the execution of this potential carbon flux can generate either a source or a sink of new mineral nitrogen (:math:`{NF}_{pot\_min,u}`\ :math:`{}_{\rightarrow}`\ :math:`{}_{d}`, gN m\ :sup:`-2` s\ :sup:`-1`). The governing equation (Thornton and Rosenbloom, 2005) is:

.. math::
   :label: 21.12)

   NF_{pot\_ min,\, u\to d} =\frac{CF_{pot,\, u} \left(1-rf_{u} -\frac{CN_{d} }{CN_{u} } \right)}{CN_{d} }

where :math:`{rf}_{u}` is the respiration fraction for fluxes leaving the upstream pool, :math:`{CN}_{u}` and :math:`{CN}_{d}` are the C:N ratios for upstream and downstream pools, respectively Negative values of :math:`{NF}_{pot\_min,u}`\ :math:`{}_{\rightarrow}`\ :math:`{}_{d}` indicate that the decomposition flux results in a source of new mineral nitrogen, while positive values indicate that the potential decomposition flux results in a sink (demand) for mineral nitrogen.

Following from the general case, potential carbon fluxes leaving individual pools in the decomposition cascade, for the example of the CLM-CN pool structure, are given as:

.. math::
   :label: 21.13)

   CF_{pot,\, Lit1} ={CS_{Lit1} k_{Lit1} r_{total} \mathord{\left/ {\vphantom {CS_{Lit1} k_{Lit1} r_{total}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 21.14)

   CF_{pot,\, Lit2} ={CS_{Lit2} k_{Lit2} r_{total} \mathord{\left/ {\vphantom {CS_{Lit2} k_{Lit2} r_{total}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 21.15)

   CF_{pot,\, Lit3} ={CS_{Lit3} k_{Lit3} r_{total} \mathord{\left/ {\vphantom {CS_{Lit3} k_{Lit3} r_{total}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 21.16)

   CF_{pot,\, SOM1} ={CS_{SOM1} k_{SOM1} r_{total} \mathord{\left/ {\vphantom {CS_{SOM1} k_{SOM1} r_{total}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 21.17)

   CF_{pot,\, SOM2} ={CS_{SOM2} k_{SOM2} r_{total} \mathord{\left/ {\vphantom {CS_{SOM2} k_{SOM2} r_{total}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 21.18)

   CF_{pot,\, SOM3} ={CS_{SOM3} k_{SOM3} r_{total} \mathord{\left/ {\vphantom {CS_{SOM3} k_{SOM3} r_{total}  \Delta t}} \right.} \Delta t}

.. math::
   :label: 21.19)

   CF_{pot,\, SOM4} ={CS_{SOM4} k_{SOM4} r_{total} \mathord{\left/ {\vphantom {CS_{SOM4} k_{SOM4} r_{total}  \Delta t}} \right.} \Delta t}

where the factor (1/:math:`\Delta`\ *t*) is included because the rate constant is calculated for the entire timestep (Eqs. and ), but the convention is to express all fluxes on a per-second basis. Potential mineral nitrogen fluxes associated with these decomposition steps are, again for the example of the CLM-CN pool structure (the CENTURY structure will be similar but without the different terminal step):

.. math::
   :label: ZEqnNum934998

   NF_{pot\_ min,\, Lit1\to SOM1} ={CF_{pot,\, Lit1} \left(1-rf_{Lit1} -\frac{CN_{SOM1} }{CN_{Lit1} } \right)\mathord{\left/ {\vphantom {CF_{pot,\, Lit1} \left(1-rf_{Lit1} -\frac{CN_{SOM1} }{CN_{Lit1} } \right) CN_{SOM1} }} \right.} CN_{SOM1} }

.. math::
   :label: 21.21)

   NF_{pot\_ min,\, Lit2\to SOM2} ={CF_{pot,\, Lit2} \left(1-rf_{Lit2} -\frac{CN_{SOM2} }{CN_{Lit2} } \right)\mathord{\left/ {\vphantom {CF_{pot,\, Lit2} \left(1-rf_{Lit2} -\frac{CN_{SOM2} }{CN_{Lit2} } \right) CN_{SOM2} }} \right.} CN_{SOM2} }

.. math::
   :label: 21.22)

   NF_{pot\_ min,\, Lit3\to SOM3} ={CF_{pot,\, Lit3} \left(1-rf_{Lit3} -\frac{CN_{SOM3} }{CN_{Lit3} } \right)\mathord{\left/ {\vphantom {CF_{pot,\, Lit3} \left(1-rf_{Lit3} -\frac{CN_{SOM3} }{CN_{Lit3} } \right) CN_{SOM3} }} \right.} CN_{SOM3} }

.. math::
   :label: 21.23)

   NF_{pot\_ min,\, SOM1\to SOM2} ={CF_{pot,\, SOM1} \left(1-rf_{SOM1} -\frac{CN_{SOM2} }{CN_{SOM1} } \right)\mathord{\left/ {\vphantom {CF_{pot,\, SOM1} \left(1-rf_{SOM1} -\frac{CN_{SOM2} }{CN_{SOM1} } \right) CN_{SOM2} }} \right.} CN_{SOM2} }

.. math::
   :label: 21.24)

   NF_{pot\_ min,\, SOM2\to SOM3} ={CF_{pot,\, SOM2} \left(1-rf_{SOM2} -\frac{CN_{SOM3} }{CN_{SOM2} } \right)\mathord{\left/ {\vphantom {CF_{pot,\, SOM2} \left(1-rf_{SOM2} -\frac{CN_{SOM3} }{CN_{SOM2} } \right) CN_{SOM3} }} \right.} CN_{SOM3} }

.. math::
   :label: 21.25)

   NF_{pot\_ min,\, SOM3\to SOM4} ={CF_{pot,\, SOM3} \left(1-rf_{SOM3} -\frac{CN_{SOM4} }{CN_{SOM3} } \right)\mathord{\left/ {\vphantom {CF_{pot,\, SOM3} \left(1-rf_{SOM3} -\frac{CN_{SOM4} }{CN_{SOM3} } \right) CN_{SOM4} }} \right.} CN_{SOM4} }

.. math::
   :label: ZEqnNum473594

   NF_{pot\_ min,\, SOM4} =-{CF_{pot,\, SOM4} \mathord{\left/ {\vphantom {CF_{pot,\, SOM4}  CN_{SOM4} }} \right.} CN_{SOM4} }

where the special form of Eq. arises because there is no SOM pool downstream of SOM4 in the converging cascade: all carbon fluxes leaving that pool are assumed to be in the form of respired CO\ :sub:`2`, and all nitrogen fluxes leaving that pool are assumed to be sources of new mineral nitrogen.

Steps in the decomposition cascade that result in release of new mineral nitrogen (mineralization fluxes) are allowed to proceed at their potential rates, without modification for nitrogen availability. Steps that result in an uptake of mineral nitrogen (immobilization fluxes) are subject to rate limitation, depending on the availability of mineral nitrogen, the total immobilization demand, and the total demand for soil mineral nitrogen to support new plant growth. The potential mineral nitrogen fluxes from Eqs. - are evaluated, summing all the positive fluxes to generate the total potential nitrogen immobilization flux (:math:`{NF}_{immob\_demand}`, gN m\ :sup:`-2` s\ :sup:`-1`), and summing absolute values of all the negative fluxes to generate the total nitrogen mineralization flux (:math:`{NF}_{gross\_nmin}`, gN m\ :sup:`-2` s\ :sup:`-1`). Since :math:`{NF}_{griss\_nmin}` is a source of new mineral nitrogen to the soil mineral nitrogen pool it is not limited by the availability of soil mineral nitrogen, and is therefore an actual as opposed to a potential flux.

N Competition between plant uptake and soil immobilization fluxes
----------------------------------------------------------------------

Once :math:`{NF}_{immob\_demand }` and :math:`{NF}_{nit\_demand }` for each layer *j* are known, the competition between plant and microbial nitrogen demand can be resolved. Mineral nitrogen in the soil pool (:math:`{NS}_{sminn}`, gN m\ :sup:`-2`) at the beginning of the timestep is considered the available supply.

Here, the :math:`{NF}_{plant\_demand}` is the theoretical maximum demand for nitrogen by plants to meet the entire carbon uptake given an N cost of zero (and therefore represents the upper bound on N requirements). N uptake costs that are :math:`>` 0 imply that the plant will take up less N that it demands, ultimately. However, given the heuristic nature of the N competition algorithm, this discrepancy is not explicitly resolved here.

The hypothetical plant nitrogen demand from the soil mineral pool is distributed between layers in proportion to the profile of available mineral N:

.. math::
   :label: 21.291

   NF_{plant\_ demand,j} =  NF_{plant\_ demand} NS_{sminn\_ j}  / \sum _{j=1}^{nj}NS_{sminn,j}

Plants first compete for ammonia (NH4). For each soil layer (*j*), we calculate the total NH4 demand as:

.. math::
   :label: 21.292

   NF_{total\_ demand_nh4,j}  = NF_{immob\_ demand,j}  + NF_{immob\_ demand,j} + NF_{nit\_ demand,j}

where If :math:`{NF}_{total\_demand,j}`\ :math:`\Delta`\ *t* :math:`<` :math:`{NS}_{sminn,j}`, then the available pool is large enough to meet both the maximum plant and microbial demand, then immobilization proceeds at the maximum rate.

.. math::
   :label: 21.29)

   f_{immob\_demand,j} = 1.0

where :math:`{f}_{immob\_demand,j}` is the fraction of potential immobilization demand that can be met given current supply of mineral nitrogen in this layer. We also set the actual nitrification flux to be the same as the potential flux (:math:`NF_{nit}` = :math:`NF_{nit\_ demand}`).

If :math:`{NF}_{total\_demand,j} \Delta t \mathrm{\ge} {NS}_{sminn,j}`, then there is not enough mineral nitrogen to meet the combined demands for plant growth and heterotrophic immobilization, immobilization is reduced proportional to the discrepancy, by :math:`f_{immob\_ demand,j}`, where

.. math::
   :label: 21.30)

   f_{immob\_ demand,j} = \frac{NS_{sminn,j} }{\Delta t\, NF_{total\_ demand,j} }

The N available to the FUN model for plant uptake (:math:`{NF}_ {plant\_ avail\_ sminn}` (gN m\ :sup:`-2`), which determines both the cost of N uptake, and the absolute limit on the N which is available for acquisition, is calculated as the total mineralized pool minus the actual immobilized flux:

.. math::
   :label: 21.311)

   NF_{plant\_ avail\_ sminn,j} = NS_{sminn,j} - f_{immob\_demand} NF_{immob\_ demand,j}

This treatment of competition for nitrogen as a limiting resource is referred to a demand-based competition, where the fraction of the available resource that eventually flows to a particular process depends on the demand from that process in comparison to the total demand from all processes. Processes expressing a greater demand acquire a larger vfraction of the available resource.

Final Decomposition Fluxes
-------------------------------

With :math:`{f}_{immob\_demand}` known, final decomposition fluxes can be calculated. Actual carbon fluxes leaving the individual litter and SOM pools, again for the example of the CLM-CN pool structure (the CENTURY structure will be similar but, again without the different terminal step), are calculated as:

.. math::
   :label: 21.32)

   CF_{Lit1} =\left\{\begin{array}{l} {CF_{pot,\, Lit1} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, Lit1\to SOM1} >0} \\ {CF_{pot,\, Lit1} \qquad {\rm for\; }NF_{pot\_ min,\, Lit1\to SOM1} \le 0} \end{array}\right\}

.. math::
   :label: 21.33)

   CF_{Lit2} =\left\{\begin{array}{l} {CF_{pot,\, Lit2} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, Lit2\to SOM2} >0} \\ {CF_{pot,\, Lit2} \qquad {\rm for\; }NF_{pot\_ min,\, Lit2\to SOM2} \le 0} \end{array}\right\}

.. math::
   :label: 21.34)

   CF_{Lit3} =\left\{\begin{array}{l} {CF_{pot,\, Lit3} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, Lit3\to SOM3} >0} \\ {CF_{pot,\, Lit3} \qquad {\rm for\; }NF_{pot\_ min,\, Lit3\to SOM3} \le 0} \end{array}\right\}

.. math::
   :label: 21.35)

   CF_{SOM1} =\left\{\begin{array}{l} {CF_{pot,\, SOM1} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, SOM1\to SOM2} >0} \\ {CF_{pot,\, SOM1} \qquad {\rm for\; }NF_{pot\_ min,\, SOM1\to SOM2} \le 0} \end{array}\right\}

.. math::
   :label: 21.36)

   CF_{SOM2} =\left\{\begin{array}{l} {CF_{pot,\, SOM2} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, SOM2\to SOM3} >0} \\ {CF_{pot,\, SOM2} \qquad {\rm for\; }NF_{pot\_ min,\, SOM2\to SOM3} \le 0} \end{array}\right\}

.. math::
   :label: 21.37)

   CF_{SOM3} =\left\{\begin{array}{l} {CF_{pot,\, SOM3} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, SOM3\to SOM4} >0} \\ {CF_{pot,\, SOM3} \qquad {\rm for\; }NF_{pot\_ min,\, SOM3\to SOM4} \le 0} \end{array}\right\}

.. math::
   :label: 21.38)

   CF_{SOM4} =CF_{pot,\, SOM4}

Heterotrophic respiration fluxes (losses of carbon as CO\ :sub:`2` to the atmosphere) are:

.. math::
   :label: 21.39)

   CF_{Lit1,\, HR} =CF_{Lit1} rf_{Lit1}

.. math::
   :label: 21.40)

   CF_{Lit2,\, HR} =CF_{Lit2} rf_{Lit2}

.. math::
   :label: 21.41)

   CF_{Lit3,\, HR} =CF_{Lit3} rf_{Lit3}

.. math::
   :label: 21.42)

   CF_{SOM1,\, HR} =CF_{SOM1} rf_{SOM1}

.. math::
   :label: 21.43)

   CF_{SOM2,\, HR} =CF_{SOM2} rf_{SOM2}

.. math::
   :label: 21.44)

   CF_{SOM3,\, HR} =CF_{SOM3} rf_{SOM3}

.. math::
   :label: 21.45)

   CF_{SOM4,\, HR} =CF_{SOM4} rf_{SOM4}

Transfers of carbon from upstream to downstream pools in the decomposition cascade are given as:

.. math::
   :label: 21.46)

   CF_{Lit1,\, SOM1} =CF_{Lit1} \left(1-rf_{Lit1} \right)

.. math::
   :label: 21.47)

   CF_{Lit2,\, SOM2} =CF_{Lit2} \left(1-rf_{Lit2} \right)

.. math::
   :label: 21.48)

   CF_{Lit3,\, SOM3} =CF_{Lit3} \left(1-rf_{Lit3} \right)

.. math::
   :label: 21.49)

   CF_{SOM1,\, SOM2} =CF_{SOM1} \left(1-rf_{SOM1} \right)

.. math::
   :label: 21.50)

   CF_{SOM2,\, SOM3} =CF_{SOM2} \left(1-rf_{SOM2} \right)

.. math::
   :label: 21.51)

   CF_{SOM3,\, SOM4} =CF_{SOM3} \left(1-rf_{SOM3} \right)

In accounting for the fluxes of nitrogen between pools in the decomposition cascade and associated fluxes to or from the soil mineral nitrogen pool, the model first calculates a flux of nitrogen from an upstream pool to a downstream pool, then calculates a flux either from the soil mineral nitrogen pool to the downstream pool (immobilization or from the downstream pool to the soil mineral nitrogen pool (mineralization). Transfers of nitrogen from upstream to downstream pools in the decomposition cascade are given as:

.. math::
   :label: 21.52)

   NF_{Lit1,\, SOM1} ={CF_{Lit1} \mathord{\left/ {\vphantom {CF_{Lit1}  CN_{Lit1} }} \right.} CN_{Lit1} }

.. math::
   :label: 21.53)

   NF_{Lit2,\, SOM2} ={CF_{Lit2} \mathord{\left/ {\vphantom {CF_{Lit2}  CN_{Lit2} }} \right.} CN_{Lit2} }

.. math::
   :label: 21.54)

   NF_{Lit3,\, SOM3} ={CF_{Lit3} \mathord{\left/ {\vphantom {CF_{Lit3}  CN_{Lit3} }} \right.} CN_{Lit3} }

.. math::
   :label: 21.55)

   NF_{SOM1,\, SOM2} ={CF_{SOM1} \mathord{\left/ {\vphantom {CF_{SOM1}  CN_{SOM1} }} \right.} CN_{SOM1} }

.. math::
   :label: 21.56)

   NF_{SOM2,\, SOM3} ={CF_{SOM2} \mathord{\left/ {\vphantom {CF_{SOM2}  CN_{SOM2} }} \right.} CN_{SOM2} }

.. math::
   :label: 21.57)

   NF_{SOM3,\, SOM4} ={CF_{SOM3} \mathord{\left/ {\vphantom {CF_{SOM3}  CN_{SOM3} }} \right.} CN_{SOM3} }

Corresponding fluxes to or from the soil mineral nitrogen pool depend on whether the decomposition step is an immobilization flux or a mineralization flux:

.. math::
   :label: 21.58)

   NF_{sminn,\, Lit1\to SOM1} =\left\{\begin{array}{l} {NF_{pot\_ min,\, Lit1\to SOM1} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, Lit1\to SOM1} >0} \\ {NF_{pot\_ min,\, Lit1\to SOM1} \qquad {\rm for\; }NF_{pot\_ min,\, Lit1\to SOM1} \le 0} \end{array}\right\}

.. math::
   :label: 21.59)

   NF_{sminn,\, Lit2\to SOM2} =\left\{\begin{array}{l} {NF_{pot\_ min,\, Lit2\to SOM2} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, Lit2\to SOM2} >0} \\ {NF_{pot\_ min,\, Lit2\to SOM2} \qquad {\rm for\; }NF_{pot\_ min,\, Lit2\to SOM2} \le 0} \end{array}\right\}

.. math::
   :label: 21.60)

   NF_{sminn,\, Lit3\to SOM3} =\left\{\begin{array}{l} {NF_{pot\_ min,\, Lit3\to SOM3} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, Lit3\to SOM3} >0} \\ {NF_{pot\_ min,\, Lit3\to SOM3} \qquad {\rm for\; }NF_{pot\_ min,\, Lit3\to SOM3} \le 0} \end{array}\right\}

.. math::
   :label: 21.61)

   NF_{sminn,SOM1\to SOM2} =\left\{\begin{array}{l} {NF_{pot\_ min,\, SOM1\to SOM2} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, SOM1\to SOM2} >0} \\ {NF_{pot\_ min,\, SOM1\to SOM2} \qquad {\rm for\; }NF_{pot\_ min,\, SOM1\to SOM2} \le 0} \end{array}\right\}

.. math::
   :label: 21.62)

   NF_{sminn,SOM2\to SOM3} =\left\{\begin{array}{l} {NF_{pot\_ min,\, SOM2\to SOM3} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, SOM2\to SOM3} >0} \\ {NF_{pot\_ min,\, SOM2\to SOM3} \qquad {\rm for\; }NF_{pot\_ min,\, SOM2\to SOM3} \le 0} \end{array}\right\}

.. math::
   :label: 21.63)

   NF_{sminn,SOM3\to SOM4} =\left\{\begin{array}{l} {NF_{pot\_ min,\, SOM3\to SOM4} f_{immob\_ demand} \qquad {\rm for\; }NF_{pot\_ min,\, SOM3\to SOM4} >0} \\ {NF_{pot\_ min,\, SOM3\to SOM4} \qquad {\rm for\; }NF_{pot\_ min,\, SOM3\to SOM4} \le 0} \end{array}\right\}

.. math::
   :label: 21.64)

   NF_{sminn,\, SOM4} =NF_{pot\_ min,\, SOM4}

Vertical Distribution and Transport of Decomposing C and N pools
---------------------------------------------------------------------

Additional terms are needed to calculate the vertically-resolved soil C and N budget: the initial vertical distribution of C and N from PFTs delivered to the litter and CWD pools, and the vertical transport of C and N pools.

For initial vertical inputs, CLM uses separate profiles for aboveground (leaf, stem) and belowground (root) inputs. Aboveground inputs are given a single exponential with default e-folding depth = 0.1m. Belowground inputs are distributed according to rooting profiles with default values based on the Jackson et al. (1996) exponential parameterization.

Vertical mixing is accomplished by an advection-diffusion equation. The goal of this is to consider slow, soild- and adsorbed-phase transport due to bioturbation, cryoturbation, and erosion. Faster aqueous-phase transport is not included in CLM, but has been developed as part of the CLM-BeTR suite of parameterizations (Tang and Riley 2013). The default value of the advection term is 0 cm/yr, such that transport is purely diffusive. Diffusive transport differs in rate between permafrost soils (where cryoturbation is the dominant transport term) and non-permafrost soils (where bioturbation dominates). For permafrost soils, a parameterization based on that of :ref:`Koven et al. (2009) <Kovenetal2009>` is used: the diffusivity parameter is constant through the active layer, and decreases linearly from the base of the active layer to zero at a set depth (default 3m); the default permafrost diffusivity is 5 cm\ :sup:`2`/yr. For non-permafrost soils, the default diffusivity is 1 cm\ :sup:`2`/yr.

Model Equilibration and its Acceleration
-----------------------------------------
For transient experiments, it is usually assumed that the carbon cycle is starting from a point of relatively close equilibrium, i.e. that productivity is balanced by ecosystem carbon losses through respiratory and disturbance pathways. In order to satisfy this assumption, the model is generally run until the productivity and loss terms find a stable long-term equilibrium; at this point the model is considered 'spun up'.

Because of the coupling between the slowest SOM pools and productivity through N downregulation of photosynthesis, equilibration of the model for initialization purposes will take an extremely long time in the standard mode. This is particularly true for the CENTURY-based decomposition cascade, which includes a passive pool. In order to rapidly equilibrate the model, a modified version of the "accelerated decomposition" :ref:`(Thornton and Rosenbloon, 2005) <ThorntonRosenbloom2005>` is used. The fundamental idea of this approach is to allow fluxes between the various pools (both turnover-defined and vertically-defined fluxes) adjust rapidly, while keeping the pool sizes themselves small so that they can fill quickly To do this, the base decomposition rate :math:`{k}_{i}` for each pool *i* is accelerated by a term :math:`{a}_{i}` such that the slow pools are collapsed onto an approximately annual timescale :ref:`Koven et al. (2013) <Kovenetal2013>`. Accelerating the pools beyond this timescale distorts the seasonal and/or diurnal cycles of decomposition and N mineralization, thus leading to a substantially different ecosystem productivity than the full model. For the vertical model, the vertical transport terms are also accelerated by the same term :math:`{a}_{i}`, as is the radioactive decay when :math:`{}^{14}`\ C is enabled, following the same principle of keeping fluxes between pools (or fluxes lost to decay close to the full model while keeping the pools sizes small. When leaving the accelerated decomposition mode, the concentration of C and N in pools that had been accelerated are multiplied by the same term :math:`{a}_{i}`, to bring the model into approximate equilibrium Note that in CLM, the model can also transition into accelerated decomposition mode from the standard mode (by dividing the pools by :math:`{a}_{i}`), and that the transitions into and out of accelerated decomposition mode are handled automatically by CLM upon loading from restart files (which preserve information about the mode of the model when restart files were written).

The base acceleration terms for the two decomposition cascades are shown in Tables 15.1 and 15.3. In addition to the base terms, CLM5 also includes a geographic term to the acceleration in order to apply larger values to high-latitude systems, where decomposition rates are particularly slow and thus equilibration can take significantly longer than in temperate or tropical climates. This geographic term takes the form of a logistic equation, where :math:`{a}_{i}` is equal to the product of the base acceleration term and :math:`{a}_{l}` below:

.. math::
   :label: 21.65)

    a_l = 1 + 50 / \left ( 1 + exp \left (-0.1 * (abs(latitude) -
    60 ) \right ) \right )

