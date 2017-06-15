.. _rst_Stomatal Resistance and Photosynthesis:

Stomatal Resistance and Photosynthesis
=========================================

Summary of CLM5.0 updates relative to the CLM4.5
-----------------------------------------------------

We describe here the complete photosynthesis and stomatal conductance parameterizations that
appear in CLM5.0. Corresponding information for CLM4.5 appeared in the
CLM4.5 Technical Note (:ref:`Oleson et al. 2013 <Olesonetal2013>`).

CLM5 includes the following new changes to photosynthesis and stomatal conductance:

- Default stomatal conductance calculation uses the Medlyn conductance model

- :math:`V_{cmax}` is predicted by the LUNA model (Chapter :numref:`rst_Photosynthetic Capacity`)

- Water stress is applied by the hydraulic conductance model (Chapter :numref:`rst_Plant Hydraulics`) 


Introduction
-----------------------

Leaf stomatal resistance, which is needed for the water vapor flux
(Chapter :numref:`rst_Momentum, Sensible Heat, and Latent Heat Fluxes`), 
is coupled to leaf photosynthesis similar to Collatz et al.
(:ref:`1991 <Collatzetal1991>`, :ref:`1992 <Collatzetal1992>`). These equations are solved separately for sunlit and
shaded leaves using average absorbed photosynthetically active radiation
for sunlit and shaded leaves
[:math:`\phi ^{sun}` ,\ :math:`\phi ^{sha}`  W m\ :sup:`-2`
(section :numref:`Solar Fluxes`)] to give sunlit and shaded stomatal resistance
(:math:`r_{s}^{sun}` ,\ :math:`r_{s}^{sha}` s m\ :sup:`-1`) and
photosynthesis (:math:`A^{sun}` ,\ :math:`A^{sha}`  µmol CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`). Canopy
photosynthesis is :math:`A^{sun} L^{sun} +A^{sha} L^{sha}` , where
:math:`L^{sun}`  and :math:`L^{sha}`  are the sunlit and shaded leaf
area indices (section :numref:`Solar Fluxes`). Canopy conductance is
:math:`\frac{1}{r_{b} +r_{s}^{sun} } L^{sun} +\frac{1}{r_{b} +r_{s}^{sha} } L^{sha}` ,
where :math:`r_{b}`  is the leaf boundary layer resistance (section
:numref:`Sensible and Latent Heat Fluxes and Temperature for Vegetated Surfaces`). 
The implementation is described by Bonan et al. (:ref:`2011 <Bonanetal2011>`), though different 
methods of calculating stomatal conductance and :math:`V_{cmax}` are used in CLM5.

.. _Stomatal resistance:

Stomatal resistance
-----------------------

Leaf stomatal resistance is calculated from the Ball-Berry conductance
model as described by Collatz et al. (1991) and implemented in global
climate models (Sellers et al. 1996). The model relates stomatal
conductance (i.e., the inverse of resistance) to net leaf
photosynthesis, scaled by the relative humidity at the leaf surface and
the CO\ :sub:`1` concentration at the leaf surface. Leaf stomatal
resistance is

.. math::
   :label: ZEqnNum230514 

   \frac{1}{r_{s} } =g_{s} =m\frac{A_{n} }{{c_{s} \mathord{\left/ {\vphantom {c_{s}  P_{atm} }} \right. \kern-\nulldelimiterspace} P_{atm} } } h_{s} +b\, \beta _{t}

where :math:`r_{s}` is leaf stomatal resistance (s m\ :sup:`2`
:math:`\mu`\ mol\ :sup:`-1`), :math:`m` is a plant functional
type dependent parameter (Table 8.1), :math:`A_{n}` is leaf net
photosynthesis (:math:`\mu`\ mol CO\ :sub:`2` m\ :sup:`-2`
s\ :sup:`-1`), :math:`c_{s}` is the CO\ :sub:`2` partial
pressure at the leaf surface (Pa), :math:`P_{atm}` is the atmospheric
pressure (Pa), :math:`h_{s} =e_{s} /e_{i}` is the leaf surface humidity
with :math:`e_{s}` the vapor pressure at the leaf surface (Pa) and
:math:`e_{i}` the saturation vapor pressure (Pa) inside the leaf at the
vegetation temperature\ :math:`T_{v}`, and :math:`b` is the minimum
stomatal conductance (:math:`\mu` mol m :sup:`-2`
s\ :sup:`-1`). Parameter values are :math:`m=9` for
C\ :sub:`3` plants and :math:`m=4` for C\ :sub:`4` plants
(Collatz et al. 1991, 1992, Sellers et al. 1996). Sellers et al. (1996)
used :math:`b=10000` for C\ :sub:`3` plants and
:math:`b=40000` for C\ :sub:`4` plants, also used here.
Photosynthesis is calculated for sunlit (:math:`A^{sun}`) and shaded
(:math:`A^{sha}`) leaves to give :math:`r_{s}^{sun}` and
:math:`r_{s}^{sha}`. Additionally, soil water influences stomatal
resistance directly by multiplying the minimum conductance by a soil
water stress function :math:`\beta _{t}` (which ranges from 0 to 1) and
also indirectly through :math:`A_{n}`, as in (Sellers et al. 1996).

Resistance is converted from units of 
s m\ :sup:`2` :math:`\mu` mol\ :sub:`-1` to  s m\ :sup:`-1` as: 
1 s m\ sup:`-1` = :math:`1\times 10^{-9} R_{gas} \frac{\theta _{atm} }{P_{atm} }`
:math:`\mu` mol\ :sup:`-1` m\ :sup:`2` s, where :math:`R_{gas}` is the universal gas constant (J K\ :sup:`-1`
kmol\ :sup:`-1`) (Table 2.6) and :math:`\theta _{atm}` is the
atmospheric potential temperature (K).

.. _Table Plant functional type (PFT) photosynthetic parameters:

.. table:: Plant functional type (PFT) photosynthetic parameters.

 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | PFT                              | m   | :math:`\alpha`     | :math:`CN_{L}`    | :math:`F_{LNR}`    | :math:`SLA_{0}`    | :math:`\psi _{o}`    | :math:`\psi _{c}`    | :math:`{V}_{cmax25}`      |
 +==================================+=====+====================+===================+====================+====================+======================+======================+===========================+
 | NET Temperate                    | 9   | –                  | 35                | 0.0509             | 0.010              | -66000               | -255000              | 62.5                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | NET Boreal                       | 9   | –                  | 40                | 0.0466             | 0.008              | -66000               | -255000              | 62.6                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | NDT Boreal                       | 9   | –                  | 25                | 0.0546             | 0.024              | -66000               | -255000              | 39.1                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BET Tropical                     | 9   | –                  | 30                | 0.0461             | 0.012              | -66000               | -255000              | 55.0                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BET temperate                    | 9   | –                  | 30                | 0.0515             | 0.012              | -66000               | -255000              | 61.5                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BDT tropical                     | 9   | –                  | 25                | 0.0716             | 0.030              | -35000               | -224000              | 41.0                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BDT temperate                    | 9   | –                  | 25                | 0.1007             | 0.030              | -35000               | -224000              | 57.7                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BDT boreal                       | 9   | –                  | 25                | 0.1007             | 0.030              | -35000               | -224000              | 57.7                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BES temperate                    | 9   | –                  | 30                | 0.0517             | 0.012              | -83000               | -428000              | 61.7                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BDS temperate                    | 9   | –                  | 25                | 0.0943             | 0.030              | -83000               | -428000              | 54.0                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | BDS boreal                       | 9   | –                  | 25                | 0.0943             | 0.030              | -83000               | -428000              | 54.0                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | C\ :sub:`3` arctic grass         | 9   | –                  | 25                | 0.1365             | 0.030              | -74000               | -275000              | 78.2                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | C\ :sub:`3` grass                | 9   | –                  | 25                | 0.1365             | 0.030              | -74000               | -275000              | 78.2                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | C\ :sub:`4` grass                | 4   | 0.05               | 25                | 0.0900             | 0.030              | -74000               | -275000              | 51.6                      |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Crop R                           | 9   | –                  | 25                | 0.1758             | 0.030              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Crop I                           | 9   | –                  | 25                | 0.1758             | 0.030              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Corn R                           | 4   | 0.05               | 25                | 0.2930             | 0.050              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Corn I                           | 4   | 0.05               | 25                | 0.2930             | 0.050              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Temp Cereal R                    | 9   | –                  | 25                | 0.4102             | 0.070              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Temp Cereal I                    | 9   | –                  | 25                | 0.4102             | 0.070              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Winter Cereal R                  | 9   | –                  | 25                | 0.4102             | 0.070              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Winter Cereal I                  | 9   | –                  | 25                | 0.4102             | 0.070              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Soybean R                        | 9   | –                  | 25                | 0.4102             | 0.070              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+
 | Soybean I                        | 9   | –                  | 25                | 0.4102             | 0.070              | -74000               | -275000              | 100.7                     |
 +----------------------------------+-----+--------------------+-------------------+--------------------+--------------------+----------------------+----------------------+---------------------------+

:math:`\alpha` (mol CO\ :sub:`2` mol\ :sup:`-1` photon);
:math:`CN_{L}`  (g C g\ :sup:`-1` N); :math:`F_{LNR}`  (g N Rubisco g\ :sup:`-1` N); :math:`SLA_{0}`  (m\ :sup:`2` g\ :sup:`-1` C);
:math:`\psi _{o}`  and :math:`\psi _{c}`  (mm); 
V\ :sub:`cmax25` (:math:`\mu`\ mol m\ :sup:`-2` s\ :sup:`-1`, calculated from equation for canopy top).

.. _Photosynthesis:

Photosynthesis
------------------

Photosynthesis in C\ :sub:`3` plants is based on the model of
Farquhar et al. (1980). Photosynthesis in C\ :sub:`4` plants is
based on the model of Collatz et al. (1992). Bonan et al. (2011)
describe the implementation, modified here. In its simplest form, leaf
net photosynthesis after accounting for respiration (:math:`R_{d}` ) is

.. math::
   :label: 9.2) 

   A_{n} =\min \left(A_{c} ,A_{j} ,A_{p} \right)-R_{d} .

The RuBP carboxylase (Rubisco) limited rate of carboxylation
:math:`A_{c}`  (:math:`\mu` \ mol CO\ :sub:`2` m\ :sup:`-2`
s\ :sup:`-1`) is

.. math::
   :label: ZEqnNum141081 

   \label{ZEqnNum141081} 
   A_{c} =\left\{\begin{array}{l} {\frac{V_{c\max } \left(c_{i} -\Gamma _{\*} \right)}{c_{i} +K_{c} \left(1+{o_{i} \mathord{\left/ {\vphantom {o_{i}  K_{o} }} \right. \kern-\nulldelimiterspace} K_{o} } \right)} \qquad {\rm for\; C}_{{\rm 3}} {\rm \; plants}} \\ {V_{c\max } \qquad \qquad \qquad {\rm for\; C}_{{\rm 4}} {\rm \; plants}} \end{array}\right\}\qquad \qquad c_{i} -\Gamma _{\*} \ge 0.

The maximum rate of carboxylation allowed by the capacity to regenerate
RuBP (i.e., the light-limited rate) :math:`A_{j}`  (:math:`\mu` \ mol
CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 9.4) 

   A_{j} =\left\{\begin{array}{l} {\frac{J\left(c_{i} -\Gamma _{\*} \right)}{4c_{i} +8\Gamma _{\*} } \qquad \qquad {\rm for\; C}_{{\rm 3}} {\rm \; plants}} \\ {\alpha (4.6\phi )\qquad \qquad {\rm for\; C}_{{\rm 4}} {\rm \; plants}} \end{array}\right\}\qquad \qquad c_{i} -\Gamma _{\*} \ge 0.

The product-limited rate of carboxylation for C\ :sub:`3` plants
and the PEP carboxylase-limited rate of carboxylation for
C\ :sub:`4` plants :math:`A_{p}`  (:math:`\mu` \ mol
CO\ :sub:`2` m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: ZEqnNum104028 

   A_{p} =\left\{\begin{array}{l} {3T_{p\qquad } \qquad \qquad {\rm for\; C}_{{\rm 3}} {\rm \; plants}} \\ {k_{p} \frac{c_{i} }{P_{atm} } \qquad \qquad \qquad {\rm for\; C}_{{\rm 4}} {\rm \; plants}} \end{array}\right\}.

In these equations, :math:`c_{i}`  is the internal leaf
CO\ :sub:`2` partial pressure (Pa) and :math:`o_{i} =0.20P_{atm}` 
is the O\ :sub:`2` partial pressure (Pa). :math:`K_{c}`  and
:math:`K_{o}`  are the Michaelis-Menten constants (Pa) for
CO\ :sub:`2` and O\ :sub:`2`. :math:`\Gamma _{\*}`  (Pa) is
the CO\ :sub:`2` compensation point. :math:`V_{c\max }`  is the
maximum rate of carboxylation (µmol m\ :sup:`-2`
s\ :sup:`-1`) and :math:`J` is the electron transport rate (µmol
m\ :sup:`-2` s\ :sup:`-1`). :math:`T_{p}`  is the triose
phosphate utilization rate (µmol m\ :sup:`-2` s\ :sup:`-1`),
taken as :math:`T_{p} =0.167V_{c\max }`  so that
:math:`A_{p} =0.5V_{c\max }`  for C\ :sub:`3` plants (as in
Collatz et al. 1991). For C\ :sub:`4` plants, the light-limited
rate :math:`A_{j}`  varies with :math:`\phi`  in relation to the quantum
efficiency (:math:`\alpha =0.05` mol CO\ :sub:`2`
mol\ :sup:`-1` photon). :math:`\phi`  is the absorbed
photosynthetically active radiation (W m\ :sup:`-2`) (section
4.1), which is converted to photosynthetic photon flux assuming 4.6
:math:`\mu` \ mol photons per joule. :math:`k_{p}`  is the initial slope
of C\ :sub:`4` CO\ :sub:`2` response curve.

For C\ :sub:`3` plants, the electron transport rate depends on the
photosynthetically active radiation absorbed by the leaf. A common
expression is the smaller of the two roots of the equation

.. math::
   :label: 9.6) 

   \Theta _{PSII} J^{2} -\left(I_{PSII} +J_{\max } \right)J+I_{PSII} J_{\max } =0

where :math:`J_{\max }`  is the maximum potential rate of electron
transport (:math:`\mu`\ mol m\ :sup:`-2` s\ :sup:`-1`),
:math:`I_{PSII}`  is the light utilized in electron transport by
photosystem II (µmol m\ :sup:`-2` s\ :sup:`-1`), and
:math:`\Theta _{PSII}`  is a curvature parameter. For a given amount of
photosynthetically active radiation absorbed by a leaf :math:`\phi`  (W
m\ :sup:`-2`), converted to photosynthetic photon flux density
with 4.6 :math:`\mu`\ mol J\ :sup:`-1`, the light utilized in
electron transport is

.. math::
   :label: 9.7) 

   I_{PSII} =0.5\Phi _{PSII} (4.6\phi )

where :math:`\Phi _{PSII}`  is the quantum yield of photosystem II, and
the term 0.5 arises because one photon is absorbed by each of the two
photosystems to move one electron. Parameter values are
:math:`\Theta _{PSII}` \ = 0.7 and :math:`\Phi _{PSII}` \ = 0.85. In
calculating :math:`A_{j}`  (for both C\ :sub:`3` and
C\ :sub:`4` plants), :math:`\phi =\phi ^{sun}`  for sunlit leaves
and :math:`\phi =\phi ^{sha}`  for shaded leaves.

The model uses co-limitation as described by Collatz et al. (1991,
1992). The actual gross photosynthesis rate, :math:`A`, is given by the
smaller root of the equations

.. math::
   :label: 9.8) 

   \begin{array}{rcl} {\Theta _{cj} A_{i}^{2} -\left(A_{c} +A_{j} \right)A_{i} +A_{c} A_{j} } & {=} & {0} \\ {\Theta _{ip} A^{2} -\left(A_{i} +A_{p} \right)A+A_{i} A_{p} } & {=} & {0} \end{array} .

Values are :math:`\Theta _{cj} =0.98` and :math:`\Theta _{ip} =0.95` for
C\ :sub:`3` plants; and :math:`\Theta _{cj} =0.80`\ and
:math:`\Theta _{ip} =0.95` for C\ :sub:`4` plants.
:math:`A_{n} =A-R_{d}` .

The parameters :math:`K_{c}`, :math:`K_{o}` , and :math:`\Gamma _{*}` 
depend on temperature. Values at 25 :sup:`o` \ C are
:math:`K_{c25} ={\rm 4}0{\rm 4}.{\rm 9}\times 10^{-6} P_{atm}` ,
:math:`K_{o25} =278.4\times 10^{-3} P_{atm}` , and
:math:`\Gamma _{*25} {\rm =42}.75\times 10^{-6} P_{atm}` .
:math:`V_{c\max }` , :math:`J_{\max }` , :math:`T_{p}` , :math:`k_{p}` ,
and :math:`R_{d}`  also vary with temperature. Parameter values at 25
:math:`\circ`\ C are calculated from :math:`V_{c\max }` \ at 25
:math:`\circ`\ C: :math:`J_{\max 25} =1.97V_{c\max 25}` ,
:math:`T_{p25} =0.167V_{c\max 25}` , and
:math:`R_{d25} =0.015V_{c\max 25}`  (C\ :sub:`3`) and
:math:`R_{d25} =0.025V_{c\max 25}`  (C\ :sub:`4`). For
C\ :sub:`4` plants, :math:`k_{p25} =20000\; V_{c\max 25}` .
However, when the biogeochemistry is active, :math:`R_{d25}`  is
calculated from leaf nitrogen as :math:`R_{d25} =0.2577N_{a}` , where
:math:`N_{a}`  is the area-based leaf nitrogen concentration (g N
m\ :sup:`-2` leaf area, equation ) and 0.2577 :math:`\mu`\ mol
CO\ :sub:`2` g\ :sup:`-1` N s\ :sup:`-1` the base
respiration rate. The parameters :math:`V_{c\max 25}` ,
:math:`J_{\max 25}` , :math:`T_{p25}` , :math:`k_{p25}` , and
:math:`R_{d25}`  are scaled over the canopy for sunlit and shaded leaves
(section 8.3). In C3 plants, these are adjusted for leaf temperature
:math:`T_{v}`  (K) as:

.. math::
   :label: 9.9) 

   \begin{array}{rcl} {V_{c\max } } & {=} & {V_{c\max 25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {J_{\max } } & {=} & {J_{\max 25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {T_{p} } & {=} & {T_{p25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {R_{d} } & {=} & {R_{d25} \; f\left(T_{v} \right)f_{H} \left(T_{v} \right)} \\ {K_{c} } & {=} & {K_{c25} \; f\left(T_{v} \right)} \\ {K_{o} } & {=} & {K_{o25} \; f\left(T_{v} \right)} \\ {\Gamma _{*} } & {=} & {\Gamma _{*25} \; f\left(T_{v} \right)} \end{array}

with

.. math::
   :label: 9.10) 

   f\left(T_{v} \right)=\; \exp \left[\frac{\Delta H_{a} }{298.15\times 0.001R_{gas} } \left(1-\frac{298.15}{T_{v} } \right)\right]

and

.. math::
   :label: 9.11) 

   f_{H} \left(T_{v} \right)=\frac{1+\exp \left(\frac{298.15\Delta S-\Delta H_{d} }{298.15\times 0.001R_{gas} } \right)}{1+\exp \left(\frac{\Delta ST_{v} -\Delta H_{d} }{0.001R_{gas} T_{v} } \right)}  .

Table 8.2 list parameter values for :math:`\Delta H_{a}` ,
:math:`\Delta H_{d}` , and :math:`\Delta S`, from Bonan et al. (2011).
Because :math:`T_{p}`  as implemented here varies with
:math:`V_{c\max }` , the same temperature parameters are used for
:math:`T_{p}` . For C\ :sub:`4` plants,

.. math::
   :label: 9.12) 

   \begin{array}{l} {V_{c\max } =V_{c\max 25} \left[\frac{Q_{10} ^{(T_{v} -298.15)/10} }{f_{H} \left(T_{v} \right)f_{L} \left(T_{v} \right)} \right]} \\ {f_{H} \left(T_{v} \right)=1+\exp \left[s_{1} \left(T_{v} -s_{2} \right)\right]} \\ {f_{L} \left(T_{v} \right)=1+\exp \left[s_{3} \left(s_{4} -T_{v} \right)\right]} \end{array}

with :math:`Q_{10} =2`,
:math:`s_{1} =0.3`\ K\ :sup:`-1`
:math:`s_{2} =313.15` K,
:math:`s_{3} =0.2`\ K\ :sup:`-1`, and :math:`s_{4} =288.15` K. 
Additionally,

.. math::
   :label: 9.13) 

   R_{d} =R_{d25} \left\{\frac{Q_{10} ^{(T_{v} -298.15)/10} }{1+\exp \left[s_{5} \left(T_{v} -s_{6} \right)\right]} \right\}

with :math:`Q_{10} =2`, :math:`s_{5} =1.3`
K\ :sup:`-1` and :math:`s_{6} =328.15`\ K, and

.. math::
   :label: 9.14) 

   k_{p} =k_{p25} \, Q_{10} ^{(T_{v} -298.15)/10}

with :math:`Q_{10} =2`.

.. _Table Temperature dependence parameters for C3 photosynthesis:

.. table:: Temperature dependence parameters for C3 photosynthesis.

 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | Parameter              | :math:`\Delta H_{a}`  (J mol\ :sup:`-1`)                        | :math:`\Delta H_{d}`  (J mol\ :sup:`-1`)                        | :math:`\Delta S` (J mol\ :sup:`-1` K\ :sup:`-1`)                                             |
 +========================+=================================================================+=================================================================+==============================================================================================+
 | :math:`V_{c\max }`     | 65330                                                           | 149250                                                          | 485                                                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | :math:`J_{\max }`      | 43540                                                           | 152040                                                          | 495                                                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | :math:`T_{p}`          | 65330                                                           | 149250                                                          | 485                                                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | :math:`R_{d}`          | 46390                                                           | 150650                                                          | 490                                                                                          |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | :math:`K_{c}`          | 79430                                                           | –                                                               | –                                                                                            |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | :math:`K_{o}`          | 36380                                                           | –                                                               | –                                                                                            |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+
 | :math:`\Gamma _{\*}`   | 37830                                                           | –                                                               | –                                                                                            |
 +------------------------+-----------------------------------------------------------------+-----------------------------------------------------------------+----------------------------------------------------------------------------------------------+

The parameters in numref:`Table Temperature dependence parameters for C3 photosynthesis` 
do not allow for temperature acclimation of photosynthesis. In the model, acclimation is 
implemented as in :ref:`Kattge and Knorr (2007) <KattgeKnorr2007>`. In this parameterization, 
:math:`V_{c\max }` and :math:`J_{\max }`  vary with the plant growth temperature. This is
achieved by allowing :math:`\Delta S`\ to vary with growth temperature
according to

.. math::
   :label: 9.15) 

   \begin{array}{l} {\Delta S=668.39-1.07(T_{10} -T_{f} )\qquad \qquad {\rm for\; }V_{c\max } } \\ {\Delta S=659.70-0.75(T_{10} -T_{f} )\qquad \qquad {\rm for\; }J_{\max } } \end{array}

The effect is to cause the temperature optimum of :math:`V_{c\max }` 
and :math:`J_{\max }`  to increase with warmer temperature. In this
parameterization, :math:`\Delta H_{d}` \ = 200000,
:math:`\Delta H_{a}` \ = 72000 for :math:`V_{c\max }` , and
:math:`\Delta H_{a}` \ = 50000 for :math:`J_{\max }` . Additionally, the
ratio :math:`J_{\max 25} /V_{c\max 25}`  at 25 :sup:`o`\ C decreases with growth temperature as

.. math::
   :label: 9.16) 

   J_{\max 25} /V_{c\max 25} =2.59-0.035(T_{10} -T_{f} ).

In these acclimation functions, :math:`T_{10}`  is the 10-day mean air
temperature (K) and :math:`T_{f}`  is the freezing point of water (K).
For lack of data, :math:`T_{p}`  acclimates similar to V :sub:`cmax` . Acclimation is restricted over the temperature
range :math:`T_{10} -T_{f} \ge 11`\ :sup:`o`\ C and :math:`T_{10} -T_{f} \le 35`\ :sup:`o`\ C.

.. _V\ :sub:`cmax25` and Canopy scaling:

V\ :sub:`cmax25` and Canopy scaling
--------------------------------------------

The maximum rate of carboxylation at 25 :sup:`o`\ C varies with
foliage nitrogen concentration and specific leaf area and is calculated
as in Thornton and Zimmermann (2007). At 25ºC,

.. math::
   :label: ZEqnNum217783 

   V_{c\max 25} =N_{a} F_{LNR} F_{NR} a_{R25}

where :math:`N_{a}`  is the area-based leaf nitrogen concentration (g N
m\ :sup:`-2` leaf area), :math:`F_{LNR}`  is the fraction of leaf
nitrogen in Rubisco (g N in Rubisco g\ :sup:`-1` N),
:math:`F_{NR} =7.16` is the mass ratio of total Rubisco molecular mass
to nitrogen in Rubisco (g Rubisco g\ :sup:`-1` N in Rubisco), and
:math:`a_{R25} =60` is the specific activity of Rubisco (µmol
CO\ :sub:`2` g\ :sup:`-1` Rubisco s\ :sup:`-1`).
:math:`N_{a}`  is calculated from mass-based leaf N concentration and
specific leaf area

.. math::
   :label: ZEqnNum561340 

   N_{a} =\frac{1}{CN_{L} \; SLA_{0} }

where :math:`CN_{L}`  is the leaf carbon-to-nitrogen ratio (g C
g\ :sup:`-1` N) and :math:`SLA_{0}`  is specific leaf area at the
canopy top ( m\ :sup:`2` leaf area g\ :sup:`-1` C ). Table 8.1
lists values of :math:`F_{LNR}` , :math:`CN_{L}` , and :math:`SLA_{0}` 
for each plant functional type. :math:`F_{LNR}`  was chosen to give
:math:`V_{c\max 25}`  consistent with Kattge et al. (2009), as discussed
by Bonan et al. (2011, 2012). Table 8.1 lists derived values for
:math:`V_{c\max 25}`  at the top of the canopy using :math:`SLA_{0}` .
Tropical broadleaf evergreen trees are an exception, and a higher
:math:`V_{c\max 25}`  is used to alleviate model biases (Bonan et al.
2012).

:math:`V_{c\max 25}`  is calculated separately for sunlit and shaded
leaves using an exponential profile to area-based leaf nitrogen
(:math:`N_{a}` ), as in Bonan et al. (2011). :math:`V_{c\max 25}`  at
cumulative leaf area index :math:`x` from the canopy top scales directly
with :math:`N_{a}` , which decreases exponentially with greater
cumulative leaf area, so that

.. math::
   :label: ZEqnNum745439 

   V_{c\; \max 25}^{} \left(x\right)=V_{c\; \max 25}^{} \left(0\right)e^{-K_{n} x}

where :math:`V_{c\; \max 25}^{} \left(0\right)` is defined at the top of
the canopy using :math:`SLA_{0}` , and :math:`K_{n}`  is the decay
coefficient for nitrogen. The canopy integrated value for sunlit and
shaded leaves is

.. math::
   :label: 9.20) 

   \begin{array}{rcl} {V_{c\; \max 25}^{sun} } & {=} & {\int _{0}^{L}V_{c\; \max 25}^{} \left(x\right)f_{sun} \left(x\right)\,  dx} \\ {} & {=} & {V_{c\; \max 25}^{} \left(0\right)\left[1-e^{-\left(K_{n} +K\right)L} \right]\frac{1}{K_{n} +K} } \end{array}

.. math::
   :label: 9.21) 

   \begin{array}{rcl} {V_{c\; \max 25}^{sha} } & {=} & {\int _{0}^{L}V_{c\; \max 25}^{} \left(x\right)\left[1-f_{sun} \left(x\right)\right] \, dx} \\ {} & {=} & {V_{c\; \max 25}^{} \left(0\right)\left\{\left[1-e^{-K_{n} L} \right]\frac{1}{K_{n} } -\left[1-e^{-\left(K_{n} +K\right)L} \right]\frac{1}{K_{n} +K} \right\}} \end{array}

and the average value for the sunlit and shaded leaves is

.. math::
   :label: 9.22) 

   \bar{V}_{c\; \max 25}^{sun} ={V_{c\; \max 25}^{sun} \mathord{\left/ {\vphantom {V_{c\; \max 25}^{sun}  L^{sun} }} \right. \kern-\nulldelimiterspace} L^{sun} }

.. math::
   :label: 9.23) 

   \bar{V}_{c\; \max 25}^{sha} ={V_{c\; \max 25}^{sha} \mathord{\left/ {\vphantom {V_{c\; \max 25}^{sha}  L^{sha} }} \right. \kern-\nulldelimiterspace} L^{sha} } .

This integration is over all leaf area (:math:`L`) with
:math:`f_{sun} (x)=\exp \left(-Kx\right)` and :math:`K` the direct beam
extinction coefficient (equation 4.9). Photosynthetic parameters
:math:`J_{\max 25}` , :math:`T_{p25}` , :math:`k_{p25}` , and
:math:`R_{d25}`  scale similarly.

The value :math:`K_{n} = 0.11` chosen by Bonan et al. (2011) is
consistent with observationally-derived estimates for forests, mostly
tropical, and provides a gradient in V\ :sub:`cmax` similar to
the original CLM4 specific leaf area scaling. However, Bonan et al.
(2012) showed that the sunlit/shaded canopy parameterization does not
match an explicit multi-layer canopy parameterization. The discrepancy
arises from absorption of scattered radiation by shaded leaves and can
be tuned out with higher :math:`K_{n}` . The model uses
:math:`K_{n} =0.30` to match an explicit multi-layer canopy.

:math:`V_{c\max 25}`  additionally varies with daylength (:math:`DYL`)
using the function :math:`f(DYL)`, which introduces seasonal variation
to :math:`V_{c\max }` 

.. math::
   :label: 9.24) 

   f\left(DYL\right)=\frac{\left(DYL\right)^{2} }{\left(DYL_{\max } \right)^{2} }

with :math:`0.01\le f\left(DYL\right)\le 1`. Daylength (seconds) is
given by

.. math::
   :label: 9.25) 

   DYL=2\times 13750.9871\cos ^{-1} \left[\frac{-\sin \left(lat\right)\sin \left(decl\right)}{\cos \left(lat\right)\cos \left(decl\right)} \right]

where :math:`lat` (latitude) and :math:`decl` (declination angle) are
from section 3.3. Maximum daylength (:math:`DYL_{\max }` ) is calculated
similarly but using the maximum declination angle for present-day
orbital geometry (:math:`\pm`\ 23.4667º [:math:`\pm`\ 0.409571 radians],
positive for Northern Hemisphere latitudes and negative for Southern
Hemisphere).

.. _Soil water stress:

Soil water stress
---------------------

Soil water influences stomatal conductance directly by multiplying the
minimum conductance by a soil water stress function :math:`\beta _{t}` 
and also indirectly through :math:`A_{n}`  in the C\ :sub:`3` and
C\ :sub:`4` photosynthesis models, as in Sellers et al. (1996).
The latter effect is achieved by multiplying :math:`V_{c\; \max }`  and
:math:`R_{d}`  by :math:`\beta _{t}` .

The function :math:`\beta _{t}`  ranges from one when the soil is wet to
near zero when the soil is dry and depends on the soil water potential
of each soil layer, the root distribution of the plant functional type,
and a plant-dependent response to soil water stress

.. math::
   :label: 9.26) 

   \beta _{t} =\sum _{i}w_{i}  r_{i}

where :math:`w_{i}`  is a plant wilting factor for layer :math:`i` and
:math:`r_{i}`  is the fraction of roots in layer :math:`i`. The plant
wilting factor :math:`w_{i}`  is

.. math::
   :label: 9.27) 

   w_{i} =\left\{\begin{array}{l} {\frac{\psi _{c} -\psi _{i} }{\psi _{c} -\psi _{o} } \left[\frac{\theta _{sat,i} -\theta _{ice,i} }{\theta _{sat,i} } \right]\le 1\qquad {\rm for\; }T_{i} >T_{f} -2{\rm \; and\; }\theta _{liq,\, i} >0} \\ {0\qquad \qquad \qquad \qquad {\rm for\; }T_{i} \le T_{f} -2{\rm \; or\; }\theta _{liq,\, i} \le 0} \end{array}\right\}

where :math:`\psi _{i}`  is the soil water matric potential (mm) and
:math:`\psi _{c}`  and :math:`\psi _{o}`  are the soil water potential
(mm) when stomata are fully closed or fully open (respectively) (Table
8.1). The term in brackets scales :math:`w_{i}`  by the ratio of the
effective porosity (accounting for the ice fraction;
:math:`\theta _{sat} -\theta _{ice}` ) relative to the total porosity.

Here, the soil water matric potential :math:`\psi _{i}` \ is defined as

.. math::
   :label: 9.28) 

   \psi _{i} =\psi _{sat,i} \, s_{i}^{-B_{i} } \ge \psi _{c}

where :math:`s_{i}`  is the soil wetness for layer :math:`i` with
respect to the effective porosity and :math:`\psi _{sat,i}`  and
:math:`B_{i}`  are the saturated soil matric potential (mm) and the
Clapp and Hornberger (1978) parameter (section 7.4.1). The soil wetness
:math:`s_{i}`  is

.. math::
   :label: 9.29) 

   s_{i} =\frac{\theta _{liq,\, i} }{\theta _{sat,\, i} -\theta _{ice,\, i} } \ge 0.01

where
:math:`\theta _{ice,\, i} ={w_{ice,\, i} \mathord{\left/ {\vphantom {w_{ice,\, i}  \left(\rho _{ice} \Delta z_{i} \right)}} \right. \kern-\nulldelimiterspace} \left(\rho _{ice} \Delta z_{i} \right)} \le \theta _{sat,\, i}` \ and
:math:`\theta _{liq,\, i} ={w_{liq,\, i} \mathord{\left/ {\vphantom {w_{liq,\, i}  \left(\rho _{liq} \Delta z_{i} \right)}} \right. \kern-\nulldelimiterspace} \left(\rho _{liq} \Delta z_{i} \right)} \le \theta _{sat,\, i} -\theta _{ice,\, i}` .
:math:`w_{ice,\, i}`  and :math:`w_{liq,\, i}`  are the ice and liquid
water contents (kg m\ :sup:`-2`) (Chapter 7),
:math:`\theta _{sat,\, i}`  is the saturated volumetric water content
(section 7.4.1), :math:`\rho _{ice}`  and :math:`\rho _{liq}`  are the
densities of ice and liquid water (kg m\ :sup:`-3`) (Table 2.6),
and :math:`\Delta z_{i}`  is the soil layer thickness (m) (section 6.1).

.. _Numerical implementation photosynthesis:

Numerical implementation
----------------------------

The CO\ :sub:`2` partial pressure at the leaf surface
:math:`c_{s}`  (Pa) and the vapor pressure at the leaf surface
:math:`e_{s}`  (Pa), needed for the stomatal resistance model in
equation , and the internal leaf CO\ :sub:`2` partial pressure
:math:`c_{i}`  (Pa), needed for the photosynthesis model in equations -,
are calculated assuming there is negligible capacity to store
CO\ :sub:`2` and water vapor at the leaf surface so that

.. math::
   :label: ZEqnNum581596 

   A_{n} =\frac{c_{a} -c_{i} }{\left(1.4r_{b} +1.6r_{s} \right)P_{atm} } =\frac{c_{a} -c_{s} }{1.4r_{b} P_{atm} } =\frac{c_{s} -c_{i} }{1.6r_{s} P_{atm} }

and the transpiration fluxes are related as

.. math::
   :label: ZEqnNum323660 

   \frac{e_{a} -e_{i} }{r_{b} +r_{s} } =\frac{e_{a} -e_{s} }{r_{b} } =\frac{e_{s} -e_{i} }{r_{s} }

where :math:`r_{b}`  is leaf boundary layer resistance (s
m\ :sup:`2` :math:`\mu` \ mol\ :sup:`-1`) (section 5.3), the
terms 1.4 and 1.6 are the ratios of diffusivity of CO\ :sub:`2` to
H\ :sub:`2`\ O for the leaf boundary layer resistance and stomatal
resistance,
:math:`c_{a} ={\rm CO}_{{\rm 2}} \left({\rm mol\; mol}^{{\rm -1}} \right)P_{atm}` 
is the atmospheric CO\ :sub:`2` partial pressure (Pa) calculated
from CO\ :sub:`2` concentration (ppmv), :math:`e_{i}`  is the
saturation vapor pressure (Pa) evaluated at the leaf temperature
:math:`T_{v}` , and :math:`e_{a}`  is the vapor pressure of air (Pa).
The vapor pressure of air in the plant canopy :math:`e_{a}`  (Pa) is
determined from

.. math::
   :label: 9.33) 

   e_{a} =\frac{P_{atm} q_{s} }{0.622}

where :math:`q_{s}`  is the specific humidity of canopy air (kg
kg\ :sup:`-1`) (section 5.3). Equations and are solved for
:math:`c_{s}`  and :math:`e_{s}` 

.. math::
   :label: 9.34) 

   c_{s} =c_{a} -1.4r_{b} P_{atm} A_{n}

.. math::
   :label: ZEqnNum123611 

   e_{s} =\frac{e_{a} r_{s} +e_{i} r_{b} }{r_{b} +r_{s} }

Substitution of equation into equation gives an expression for stomatal
resistance (:math:`r_{s}` ) as a function of photosynthesis
(:math:`A_{n}` ), given here in terms of conductance with
:math:`g_{s} =1/r_{s}`  and :math:`g_{b} =1/r_{b}` 

.. math::
   :label: 9.36) 

   c_{s} g_{s}^{2} +\left[c_{s} \left(g_{b} -b\right)-m{\it A}_{n} P_{atm} \right]g_{s} -g_{b} \left[c_{s} b+mA_{n} P_{atm} {e_{a} \mathord{\left/ {\vphantom {e_{a}  e_{\*} \left(T_{v} \right)}} \right. \kern-\nulldelimiterspace} e_{\*} \left(T_{v} \right)} \right]=0.

Stomatal conductance is the larger of the two roots that satisfy the
quadratic equation. Values for :math:`c_{i}`  are given by

.. math::
   :label: 9.37) 

   c_{i} =c_{a} -\left(1.4r_{b} +1.6r_{s} \right)P_{atm} A{}_{n}

The equations for :math:`c_{i}` , :math:`c_{s}` , :math:`r_{s}` , and
:math:`A_{n}`  are solved iteratively until :math:`c_{i}`  converges.
Sun et al. (2012) pointed out that the CLM4 numerical approach does not
always converge. Therefore, the model uses a hybrid algorithm that
combines the secant method and Brent’s method to solve for
:math:`c_{i}` . The equation set is solved separately for sunlit
(:math:`A_{n}^{sun}` , :math:`r_{s}^{sun}` ) and shaded
(:math:`A_{n}^{sha}` , :math:`r_{s}^{sha}` ) leaves.

The model has an optional (though not supported) multi-layer canopy, as
described by Bonan et al. (2012). The multi-layer model is only intended
to address the non-linearity of light profiles, photosynthesis, and
stomatal conductance in the plant canopy. In the multi-layer canopy,
sunlit (:math:`A_{n}^{sun}` , :math:`r_{s}^{sun}` ) and shaded
(:math:`A_{n}^{sha}` , :math:`r_{s}^{sha}` ) leaves are explicitly
resolved at depths in the canopy using a light profile (Chapter 4). In
this case, :math:`V_{c\max 25}`  is not integrated over the canopy, but
is instead given explicitly for each canopy layer using equation . This
also uses the Lloyd et al. (2010) relationship whereby
K\ :sub:`n` scales with V\ :sub:`cmax` as

.. math::
   :label: 9.38) 

   K_{n} =\exp \left(0.00963V_{c\max } -2.43\right)

such that higher values of V\ :sub:`cmax` imply steeper declines
in photosynthetic capacity through the canopy with respect to cumulative
leaf area.
