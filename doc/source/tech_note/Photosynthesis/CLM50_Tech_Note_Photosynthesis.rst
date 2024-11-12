.. _rst_Stomatal Resistance and Photosynthesis:

Stomatal Resistance and Photosynthesis
=========================================

Summary of CLM5.0 updates relative to the CLM4.5
-----------------------------------------------------

We describe here the complete photosynthesis and stomatal conductance parameterizations that appear in CLM5.0. Corresponding information for CLM4.5 appeared in the CLM4.5 Technical Note (:ref:`Oleson et al. 2013 <Olesonetal2013>`).

CLM5 includes the following new changes to photosynthesis and stomatal conductance:

- Default stomatal conductance calculation uses the Medlyn conductance model

- :math:`V_{c,max}` and :math:`J_{max}` at 25 :sup:`\o`\ C: are now prognostic, and predicted via optimality by the LUNA model (Chapter :numref:`rst_Photosynthetic Capacity`)

- Leaf N concentration and the fraction of leaf N in Rubisco used to calculate :math:`V_{cmax25}` are determined by the LUNA model (Chapter :numref:`rst_Photosynthetic Capacity`)

- Water stress is applied by the hydraulic conductance model (Chapter :numref:`rst_Plant Hydraulics`)

Introduction
-----------------------

Leaf stomatal resistance, which is needed for the water vapor flux (Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`), is coupled to leaf photosynthesis similar to Collatz et al. (:ref:`1991 <Collatzetal1991>`, :ref:`1992 <Collatzetal1992>`). These equations are solved separately for sunlit and shaded leaves using average absorbed photosynthetically active radiation for sunlit and shaded leaves [:math:`\phi ^{sun}`,\ :math:`\phi ^{sha}` W m\ :sup:`-2` (section :numref:`Solar Fluxes`)] to give sunlit and shaded stomatal resistance (:math:`r_{s}^{sun}`,\ :math:`r_{s}^{sha}` s m\ :sup:`-1`) and photosynthesis (:math:`A^{sun}`,\ :math:`A^{sha}` µmol CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`). Canopy photosynthesis is :math:`A^{sun} L^{sun} +A^{sha} L^{sha}`, where :math:`L^{sun}` and :math:`L^{sha}` are the sunlit and shaded leaf area indices (section :numref:`Solar Fluxes`). Canopy conductance is :math:`\frac{1}{r_{b} +r_{s}^{sun} } L^{sun} +\frac{1}{r_{b} +r_{s}^{sha} } L^{sha}`, where :math:`r_{b}` is the leaf boundary layer resistance (section :numref:`Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces`).

.. _Stomatal resistance:

Stomatal resistance
-----------------------

CLM5 calculates stomatal conductance using the Medlyn stomatal conductance model (:ref:`Medlyn et al. 2011<Medlynetal2011>`). Previous versions of CLM calculated leaf stomatal resistance using the Ball-Berry conductance model as described by :ref:`Collatz et al. (1991)<Collatzetal1991>` and implemented in global climate models (:ref:`Sellers et al. 1996<Sellersetal1996>`). The Medlyn model calculates stomatal conductance (i.e., the inverse of resistance) based on net leaf photosynthesis, the leaf-to-air vapor pressure difference, and the CO\ :sub:`2` concentration at the leaf surface. Leaf stomatal resistance is:

.. math::
   :label: 9.1

   \frac{1}{r_{s} } =g_{s} = g_{o} + 1.6(1 + \frac{g_{1} }{\sqrt{D_{s}}}) \frac{A_{n} }{{c_{s} \mathord{\left/ {\vphantom {c_{s}  P_{atm} }} \right.} P_{atm} } }

where :math:`r_{s}` is leaf stomatal resistance (s m\ :sup:`2` :math:`\mu`\ mol\ :sup:`-1`), :math:`g_{o}` is the minimum stomatal conductance (:math:`\mu` mol m :sup:`-2` s\ :sup:`-1`), :math:`A_{n}` is leaf net photosynthesis (:math:`\mu`\ mol CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`), :math:`c_{s}` is the CO\ :sub:`2` partial pressure at the leaf surface (Pa), :math:`P_{atm}` is the atmospheric pressure (Pa), and :math:`D_{s}=(e_{i}-e{_s})/1000` is the leaf-to-air vapor pressure difference at the leaf surface (kPa) where :math:`e_{i}` is the saturation vapor pressure (Pa) evaluated at the leaf temperature :math:`T_{v}`, and :math:`e_{s}` is the vapor pressure at the leaf surface (Pa). :math:`g_{1}` is a plant functional type dependent parameter (:numref:`Table Plant functional type (PFT) stomatal conductance parameters`) and are the same as those used in the CABLE model (:ref:`de Kauwe et al. 2015 <deKauwe2015>`).

The value for :math:`g_{o}=100` :math:`\mu` mol m :sup:`-2` s\ :sup:`-1` for C\ :sub:`3` and C\ :sub:`4` plants. Photosynthesis is calculated for sunlit (:math:`A^{sun}`) and shaded (:math:`A^{sha}`) leaves to give :math:`r_{s}^{sun}` and :math:`r_{s}^{sha}`. Additionally, soil water influences stomatal resistance through plant hydraulic stress, detailed in the :ref:`rst_Plant Hydraulics` chapter.

Resistance is converted from units of s m\ :sup:`2` :math:`\mu` mol\ :sup:`-1` to s m\ :sup:`-1` as: 1 s m\ :sup:`-1` = :math:`1\times 10^{-9} R_{gas} \frac{\theta _{atm} }{P_{atm} }` :math:`\mu` mol\ :sup:`-1` m\ :sup:`2` s, where :math:`R_{gas}` is the universal gas constant (J K\ :sup:`-1` kmol\ :sup:`-1`) (:numref:`Table Physical constants`) and :math:`\theta _{atm}` is the atmospheric potential temperature (K).

.. _Table Plant functional type (PFT) stomatal conductance parameters:

.. table:: Plant functional type (PFT) stomatal conductance parameters.

 +----------------------------------+-------------------+
 | PFT                              |  g\ :sub:`1`      |
 +==================================+===================+
 | NET Temperate                    |        2.35       |
 +----------------------------------+-------------------+
 | NET Boreal                       |        2.35       |
 +----------------------------------+-------------------+
 | NDT Boreal                       |        2.35       |
 +----------------------------------+-------------------+
 | BET Tropical                     |        4.12       |
 +----------------------------------+-------------------+
 | BET temperate                    |        4.12       |
 +----------------------------------+-------------------+
 | BDT tropical                     |        4.45       |
 +----------------------------------+-------------------+
 | BDT temperate                    |        4.45       |
 +----------------------------------+-------------------+
 | BDT boreal                       |        4.45       |
 +----------------------------------+-------------------+
 | BES temperate                    |        4.70       |
 +----------------------------------+-------------------+
 | BDS temperate                    |        4.70       |
 +----------------------------------+-------------------+
 | BDS boreal                       |        4.70       |
 +----------------------------------+-------------------+
 | C\ :sub:`3` arctic grass         |        2.22       |
 +----------------------------------+-------------------+
 | C\ :sub:`3` grass                |        5.25       |
 +----------------------------------+-------------------+
 | C\ :sub:`4` grass                |        1.62       |
 +----------------------------------+-------------------+
 | Temperate Corn                   |        1.79       |
 +----------------------------------+-------------------+
 | Spring Wheat                     |        5.79       |
 +----------------------------------+-------------------+
 | Temperate Soybean                |        5.79       |
 +----------------------------------+-------------------+
 | Cotton                           |        5.79       |
 +----------------------------------+-------------------+
 | Rice                             |        5.79       |
 +----------------------------------+-------------------+
 | Sugarcane                        |        1.79       |
 +----------------------------------+-------------------+
 | Tropical Corn                    |        1.79       |
 +----------------------------------+-------------------+
 | Tropical Soybean                 |        5.79       |
 +----------------------------------+-------------------+
 | Miscanthus                       |        1.79       |
 +----------------------------------+-------------------+
 | Switchgrass                      |        1.79       |
 +----------------------------------+-------------------+

.. _Photosynthesis:

Photosynthesis
------------------

Photosynthesis in C\ :sub:`3` plants is based on the model of :ref:`Farquhar et al. (1980)<Farquharetal1980>`. Photosynthesis in C\ :sub:`4` plants is based on the model of :ref:`Collatz et al. (1992)<Collatzetal1992>`. :ref:`Bonan et al. (2011)<Bonanetal2011>` describe the implementation, modified here. In its simplest form, leaf net photosynthesis after accounting for respiration (:math:`R_{d}` ) is

.. math::
   :label: 9.2

   A_{n} =\min \left(A_{c} ,A_{j} ,A_{p} \right)-R_{d} .

The RuBP carboxylase (Rubisco) limited rate of carboxylation :math:`A_{c}` (:math:`\mu` \ mol CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 9.3

   A_{c} =\left\{\begin{array}{l} {\frac{V_{c\max } \left(c_{i} -\Gamma _{*} \right)}{c_{i} +K_{c} \left(1+{o_{i} \mathord{\left/ {\vphantom {o_{i}  K_{o} }} \right.} K_{o} } \right)} \qquad {\rm for\; C}_{{\rm 3}} {\rm \; plants}} \\ {V_{c\max } \qquad \qquad \qquad {\rm for\; C}_{{\rm 4}} {\rm \; plants}} \end{array}\right\}\qquad \qquad c_{i} -\Gamma _{*} \ge 0.

The maximum rate of carboxylation allowed by the capacity to regenerate RuBP (i.e., the light-limited rate) :math:`A_{j}` (:math:`\mu` \ mol CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 9.4

   A_{j} =\left\{\begin{array}{l} {\frac{J_{x}\left(c_{i} -\Gamma _{*} \right)}{4c_{i} +8\Gamma _{*} } \qquad \qquad {\rm for\; C}_{{\rm 3}} {\rm \; plants}} \\ {\alpha (4.6\phi )\qquad \qquad {\rm for\; C}_{{\rm 4}} {\rm \; plants}} \end{array}\right\}\qquad \qquad c_{i} -\Gamma _{*} \ge 0.

The product-limited rate of carboxylation for C\ :sub:`3` plants and the PEP carboxylase-limited rate of carboxylation for C\ :sub:`4` plants :math:`A_{p}` (:math:`\mu` \ mol CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 9.5

   A_{p} =\left\{\begin{array}{l} {3T_{p\qquad } \qquad \qquad {\rm for\; C}_{{\rm 3}} {\rm \; plants}} \\ {k_{p} \frac{c_{i} }{P_{atm} } \qquad \qquad \qquad {\rm for\; C}_{{\rm 4}} {\rm \; plants}} \end{array}\right\}.

In these equations, :math:`c_{i}` is the internal leaf CO\ :sub:`2` partial pressure (Pa) and :math:`o_{i} =0.20P_{atm}` is the O\ :sub:`2` partial pressure (Pa). :math:`K_{c}` and :math:`K_{o}` are the Michaelis-Menten constants (Pa) for CO\ :sub:`2` and O\ :sub:`2`. :math:`\Gamma _{*}` (Pa) is the CO\ :sub:`2` compensation point. :math:`V_{c\max }` is the maximum rate of carboxylation (µmol m\ :sup:`-2` s\ :sup:`-1`, Chapter :numref:`rst_Photosynthetic Capacity`) and :math:`J_{x}` is the electron transport rate (µmol m\ :sup:`-2` s\ :sup:`-1`). :math:`T_{p}` is the triose phosphate utilization rate (µmol m\ :sup:`-2` s\ :sup:`-1`), taken as :math:`T_{p} =0.167V_{c\max }` so that :math:`A_{p} =0.5V_{c\max }` for C\ :sub:`3` plants (as in :ref:`Collatz et al. 1992 <Collatzetal1992>`). For C\ :sub:`4` plants, the light-limited rate :math:`A_{j}` varies with :math:`\phi` in relation to the quantum efficiency (:math:`\alpha =0.05` mol CO\ :sub:`2` mol\ :sup:`-1` photon). :math:`\phi` is the absorbed photosynthetically active radiation (W m\ :sup:`-2`) (section :numref:`Solar Fluxes`), which is converted to photosynthetic photon flux assuming 4.6 :math:`\mu` \ mol photons per joule. :math:`k_{p}` is the initial slope of C\ :sub:`4` CO\ :sub:`2` response curve.

For C\ :sub:`3` plants, the electron transport rate depends on the photosynthetically active radiation absorbed by the leaf. A common expression is the smaller of the two roots of the equation

.. math::
   :label: 9.6

   \Theta _{PSII} J_{x}^{2} -\left(I_{PSII} +J_{\max } \right)J_{x}+I_{PSII} J_{\max } =0

where :math:`J_{\max }` is the maximum potential rate of electron transport (:math:`\mu`\ mol m\ :sup:`-2` s\ :sup:`-1`, Chapter :numref:`rst_Photosynthetic Capacity`), :math:`I_{PSII}` is the light utilized in electron transport by photosystem II (µmol m\ :sup:`-2` s\ :sup:`-1`), and :math:`\Theta _{PSII}` is a curvature parameter. For a given amount of photosynthetically active radiation absorbed by a leaf (:math:`\phi`, W m\ :sup:`-2`), converted to photosynthetic photon flux density with 4.6 :math:`\mu`\ mol J\ :sup:`-1`, the light utilized in electron transport is

.. math::
   :label: 9.7

   I_{PSII} =0.5\Phi _{PSII} (4.6\phi )

where :math:`\Phi _{PSII}` is the quantum yield of photosystem II, and the term 0.5 arises because one photon is absorbed by each of the two photosystems to move one electron. Parameter values are :math:`\Theta _{PSII}` \ = 0.7 and :math:`\Phi _{PSII}` \ = 0.85. In calculating :math:`A_{j}` (for both C\ :sub:`3` and C\ :sub:`4` plants), :math:`\phi =\phi ^{sun}` for sunlit leaves and :math:`\phi =\phi ^{sha}` for shaded leaves.

The model uses co-limitation as described by :ref:`Collatz et al. (1991, 1992) <Collatzetal1991>`. The actual gross photosynthesis rate, :math:`A`, is given by the smaller root of the equations

.. math::
   :label: 9.8

   \begin{array}{rcl} {\Theta _{cj} A_{i}^{2} -\left(A_{c} +A_{j} \right)A_{i} +A_{c} A_{j} } & {=} & {0} \\ {\Theta _{ip} A^{2} -\left(A_{i} +A_{p} \right)A+A_{i} A_{p} } & {=} & {0} \end{array} .

Values are :math:`\Theta _{cj} =0.98` and :math:`\Theta _{ip} =0.95` for C\ :sub:`3` plants; and :math:`\Theta _{cj} =0.80`\ and :math:`\Theta _{ip} =0.95` for C\ :sub:`4` plants. :math:`A_{i}` is the intermediate co-limited photosynthesis. :math:`A_{n} =A-R_{d}`.

The parameters :math:`K_{c}`, :math:`K_{o}`, and :math:`\Gamma` depend on temperature. Values at 25 °C are :math:`K_{c25} ={\rm 4}0{\rm 4}.{\rm 9}\times 10^{-6} P_{atm}`, :math:`K_{o25} =278.4\times 10^{-3} P_{atm}`, and :math:`\Gamma _{25} {\rm =42}.75\times 10^{-6} P_{atm}`. :math:`V_{c\max }`, :math:`J_{\max }`, :math:`T_{p}`, :math:`k_{p}`, and :math:`R_{d}` also vary with temperature.

:math:`J_{\max 25}`  at 25 :sup:`\o`\ C: is calculated by the LUNA model (Chapter :numref:`rst_Photosynthetic Capacity`)

Parameter values at 25 :sup:`\o`\ C are calculated from :math:`V_{c\max }` \ at 25
:sup:`\o`\ C:, including:
:math:`T_{p25} =0.167V_{c\max 25}`, and
:math:`R_{d25} =0.015V_{c\max 25}` (C\ :sub:`3`) and
:math:`R_{d25} =0.025V_{c\max 25}` (C\ :sub:`4`).

For C\ :sub:`4` plants, :math:`k_{p25} =20000\; V_{c\max 25}`.

However, when the biogeochemistry is active (the default mode), :math:`R_{d25}` is calculated from leaf nitrogen as described in (Chapter :numref:`rst_Plant Respiration`)

The parameters :math:`V_{c\max 25}`, :math:`J_{\max 25}`, :math:`T_{p25}`, :math:`k_{p25}`, and :math:`R_{d25}` are scaled over the canopy for sunlit and shaded leaves (section :numref:`Canopy scaling`). In C\ :sub:`3` plants, these are adjusted for leaf temperature, :math:`T_{v}` (K), as:

.. math::
   :label: 9.9

   \begin{array}{rcl} {V_{c\max } } & {=} & {V_{c\max 25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {J_{\max } } & {=} & {J_{\max 25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {T_{p} } & {=} & {T_{p25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {R_{d} } & {=} & {R_{d25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {K_{c} } & {=} & {K_{c25} \; f\left(T_{v} \right)} \\ {K_{o} } & {=} & {K_{o25} \; f\left(T_{v} \right)} \\ {\Gamma } & {=} & {\Gamma _{25} \; f\left(T_{v} \right)} \end{array}

.. math::
   :label: 9.10

   f\left(T_{v} \right)=\; \exp \left[\frac{\Delta H_{a} }{298.15\times 0.001R_{gas} } \left(1-\frac{298.15}{T_{v} } \right)\right]

and

.. math::
   :label: 9.11

   f_{H} \left(T_{v} \right)=\frac{1+\exp \left(\frac{298.15\Delta S-\Delta H_{d} }{298.15\times 0.001R_{gas} } \right)}{1+\exp \left(\frac{\Delta ST_{v} -\Delta H_{d} }{0.001R_{gas} T_{v} } \right)}  .

:numref:`Table Temperature dependence parameters for C3 photosynthesis` lists parameter values for :math:`\Delta H_{a}` and :math:`\Delta H_{d}`. :math:`\Delta S` is calculated separately for :math:`V_{c\max }` and :math:`J_{max }` to allow for temperature acclimation of photosynthesis (see equation :eq:`9.16`), and :math:`\Delta S` is 490 J mol :sup:`-1` K :sup:`-1` for :math:`R_d` (:ref:`Bonan et al. 2011<Bonanetal2011>`, :ref:`Lombardozzi et al. 2015<Lombardozzietal2015>`). Because :math:`T_{p}` as implemented here varies with :math:`V_{c\max }`, :math:`T_{p}` uses the same temperature parameters as :math:`V_{c\max}`. For C\ :sub:`4` plants,

.. math::
   :label: 9.12

   \begin{array}{l} {V_{c\max } =V_{c\max 25} \left[\frac{Q_{10} ^{(T_{v} -298.15)/10} }{f_{H} \left(T_{v} \right)f_{L} \left(T_{v} \right)} \right]} \\ {f_{H} \left(T_{v} \right)=1+\exp \left[s_{1} \left(T_{v} -s_{2} \right)\right]} \\ {f_{L} \left(T_{v} \right)=1+\exp \left[s_{3} \left(s_{4} -T_{v} \right)\right]} \end{array}

with :math:`Q_{10} =2`,
:math:`s_{1} =0.3`\ K\ :sup:`-1`
:math:`s_{2} =313.15` K,
:math:`s_{3} =0.2`\ K\ :sup:`-1`, and
:math:`s_{4} =288.15` K. Additionally,

.. math::
   :label: 9.13

   R_{d} =R_{d25} \left\{\frac{Q_{10} ^{(T_{v} -298.15)/10} }{1+\exp \left[s_{5} \left(T_{v} -s_{6} \right)\right]} \right\}

with :math:`Q_{10} =2`, :math:`s_{5} =1.3` K\ :sup:`-1` and :math:`s_{6} =328.15`\ K, and

.. math::
   :label: 9.14

   k_{p} =k_{p25} \, Q_{10} ^{(T_{v} -298.15)/10}

with :math:`Q_{10} =2`.

.. _Table Temperature dependence parameters for C3 photosynthesis:

.. table:: Temperature dependence parameters for C3 photosynthesis.

 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | Parameter              | :math:`\Delta H_{a}`  (J mol\ :sup:`-1`)                        | :math:`\Delta H_{d}`  (J mol\ :sup:`-1`)                        |
 +========================+=================================================================+=================================================================+
 | :math:`V_{c\max }`     | 72000                                                           | 200000                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | :math:`J_{\max }`      | 50000                                                           | 200000                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | :math:`T_{p}`          | 72000                                                           | 200000                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | :math:`R_{d}`          | 46390                                                           | 150650                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | :math:`K_{c}`          | 79430                                                           | –                                                               |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | :math:`K_{o}`          | 36380                                                           | –                                                               |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+
 | :math:`\Gamma _{*}`    | 37830                                                           | –                                                               |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+

In the model, acclimation is implemented as in :ref:`Kattge and Knorr (2007) <KattgeKnorr2007>`. In this parameterization, :math:`V_{c\max }` and :math:`J_{\max }` vary with the plant growth temperature. This is achieved by allowing :math:`\Delta S`\ to vary with growth temperature according to

.. math::
   :label: 9.15

   \begin{array}{l} {\Delta S=668.39-1.07(T_{10} -T_{f} )\qquad \qquad {\rm for\; }V_{c\max } } \\ {\Delta S=659.70-0.75(T_{10} -T_{f} )\qquad \qquad {\rm for\; }J_{\max } } \end{array}

The effect is to cause the temperature optimum of :math:`V_{c\max }` and :math:`J_{\max }` to increase with warmer temperatures. Additionally, the ratio :math:`J_{\max 25} /V_{c\max 25}` at 25 °C decreases with growth temperature as

.. math::
   :label: 9.16

   J_{\max 25} /V_{c\max 25} =2.59-0.035(T_{10} -T_{f} ).

In these acclimation functions, :math:`T_{10}` is the 10-day mean air temperature (K) and :math:`T_{f}` is the freezing point of water (K). For lack of data, :math:`T_{p}` acclimates similar to :math:`V_{c\max }`. Acclimation is restricted over the temperature range :math:`T_{10} -T_{f} \ge` 11°C and :math:`T_{10} -T_{f} \le` 35°C.

.. _Canopy scaling:

Canopy scaling
--------------------------------------------

When LUNA is on, the :math:`V_{c\max 25}` for sun leaves is scaled to the shaded leaves :math:`J_{\max 25}`, :math:`T_{p25}`, :math:`k_{p25}`, and :math:`R_{d25}` scale similarly.

.. math::
   :label: 9.17

   \begin{array}{rcl}
   {V_{c\max 25 sha}} & {=} & {V_{c\max 25 sha} \frac{i_{v,sha}}{i_{v,sun}}}  \\
   {J_{\max 25 sha}}  & {=} & {J_{\max 25 sun}  \frac{i_{v,sha}}{i_{v,sun}}}  \\
   {T_{p sha}}        & {=} & {T_{p sun}        \frac{i_{v,sha}}{i_{v,sun}}}  \end{array}

Where :math:`i_{v,sun}` and :math:`i_{v,sha}` are the leaf-to-canopy scaling coefficients of the twostream radiation model, calculated as

.. math::
   :label: 9.18

   i_{v,sun} = \frac{(1 - e^{-(k_{n,ext}+k_{b,ext})*lai_e)} / (k_{n,ext}+k_{b,ext})}{f_{sun}*lai_e}\\
   i_{v,sha} = \frac{(1 - e^{-(k_{n,ext}+k_{b,ext})*lai_e)} / (k_{n,ext}+k_{b,ext})}{(1 - f_{sun})*lai_e}

k_{n,ext} is the extinction coefficient for N through the canopy (0.3).  k_{b,ext} is the direct beam extinction coefficient calculated in the surface albedo routine, and :math:`f_{sun}` is the fraction of sunlit leaves, both derived from Chapter :numref:`rst_Surface Albedos`.

When LUNA is off, scaling defaults to the mechanism used in CLM4.5.

.. _Numerical implementation photosynthesis:

Numerical implementation
----------------------------

The CO\ :sub:`2` partial pressure at the leaf surface, :math:`c_{s}` (Pa), and the vapor pressure at the leaf surface, :math:`e_{s}` (Pa), needed for the stomatal resistance model in equation :eq:`9.1`, and the internal leaf CO\ :sub:`2` partial pressure :math:`c_{i}` (Pa), needed for the photosynthesis model in equations :eq:`9.3`-:eq:`9.5`, are calculated assuming there is negligible capacity to store CO\ :sub:`2` and water vapor at the leaf surface so that

.. math::
   :label: 9.19

   A_{n} =\frac{c_{a} -c_{i} }{\left(1.4r_{b} +1.6r_{s} \right)P_{atm} } =\frac{c_{a} -c_{s} }{1.4r_{b} P_{atm} } =\frac{c_{s} -c_{i} }{1.6r_{s} P_{atm} }

and the transpiration fluxes are related as

.. math::
   :label: 9.20

   \frac{e_{a} -e_{i} }{r_{b} +r_{s} } =\frac{e_{a} -e_{s} }{r_{b} } =\frac{e_{s} -e_{i} }{r_{s} }

where :math:`r_{b}` is leaf boundary layer resistance (s m\ :sup:`2` :math:`\mu` \ mol\ :sup:`-1`) (section :numref:`Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces`), the terms 1.4 and 1.6 are the ratios of diffusivity of CO\ :sub:`2` to H\ :sub:`2`\ O for the leaf boundary layer resistance and stomatal resistance, :math:`c_{a} ={\rm CO}_{{\rm 2}} \left({\rm mol\; mol}^{{\rm -1}} \right)`, :math:`P_{atm}` is the atmospheric pressure (Pa), :math:`e_{i}` is the saturation vapor pressure (Pa) evaluated at the leaf temperature :math:`T_{v}`, and :math:`e_{a}` is the vapor pressure of air (Pa). The vapor pressure of air in the plant canopy :math:`e_{a}` (Pa) is determined from

.. math::
   :label: 9.21

   e_{a} =\frac{P_{atm} q_{s} }{0.622}

where :math:`q_{s}` is the specific humidity of canopy air (kg kg\ :sup:`-1`, section :numref:`Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces`). Equations :eq:`9.19` and :eq:`9.20` are solved for :math:`c_{s}` and :math:`e_{s}`

.. math::
   :label: 9.34

   c_{s} =c_{a} -1.4r_{b} P_{atm} A_{n}

.. math::
   :label: 9.35

   e_{s} =\frac{e_{a} r_{s} +e_{i} r_{b} }{r_{b} +r_{s} }

In terms of conductance with :math:`g_{s} =1/r_{s}` and :math:`g_{b} =1/r_{b}`

.. math::
   :label: 9.36

   e_{s} =\frac{e_{a} g_{b} +e_{i} g_{s} }{g_{b} +g_{s} } .

Substitution of equation :eq:`9.36` into equation :eq:`9.1` gives an expression for the stomatal resistance (:math:`r_{s}`) as a function of photosynthesis (:math:`A_{n}` )

.. math::
   :label: 9.37

   ag_{s}^{2} + bg_{s} + c = 0

where

.. math::
   :label: 9.38

   \begin{array}{l}  a = 1 \\
   b = -[2(g_{o} * 10^{-6} + d) + \frac{(g_{1}d)^{2}}{g_{b}*10^{-6}D_{l}}] \\
   c = (g_{o}*10^{-6})^{2} + [2g_{o}*10^{-6} + d (1-\frac{g_{1}^{2}} {D_{l}})]d \end{array}

and

.. math::
   :label: 9.39

   d = \frac {1.6 A_{n}} {c_{s} / P_{atm} * 10^{6}}

   D_{l} = \frac {max(e_{i} - e_{a},50)} {1000}

Stomatal conductance, as solved by equation :eq:`9.36` (mol m :sup:`-2` s :sup:`-1`), is the larger of the two roots that satisfy the quadratic equation. Values for :math:`c_{i}` are given by

.. math::
   :label: 9.40

   c_{i} =c_{a} -\left(1.4r_{b} +1.6r_{s} \right)P_{atm} A{}_{n}

The equations for :math:`c_{i}`, :math:`c_{s}`, :math:`r_{s}`, and :math:`A_{n}` are solved iteratively until :math:`c_{i}` converges. :ref:`Sun et al. (2012)<Sunetal2012>` pointed out that the CLM4 numerical approach does not always converge. Therefore, the model uses a hybrid algorithm that combines the secant method and Brent's method to solve for :math:`c_{i}`. The equation set is solved separately for sunlit (:math:`A_{n}^{sun}`, :math:`r_{s}^{sun}` ) and shaded (:math:`A_{n}^{sha}`, :math:`r_{s}^{sha}` ) leaves.

