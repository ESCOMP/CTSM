.. _rst_Plant Respiration:

Plant Respiration
=================

Introduction
-----------------

CLM includes a fully-prognostic treatment of the terrestrial carbon and
nitrogen cycles including interactions between these cycles as mediated
by biological mechanisms of plants and soil heterotrophs. The model is
fully prognostic with respect to all carbon and nitrogen state variables
in the vegetation, litter, and soil organic matter. The seasonal timing
of new vegetation growth and litterfall is also prognostic, responding
to soil and air temperature, soil water availability, and daylength, in
varying degrees depending on a specified phenology type for each PFT
(see Chapter 14). The prognostic LAI, SAI, and vegetation heights are
utilized by the biophysical model.

Separate state variables for C and N are tracked for leaf, live stem,
dead stem, live coarse root, dead coarse root, and fine root pools
(Figure 13.1 and 13.2). Each of these pools has two corresponding
storage pools representing, respectively, short-term and long-term
storage of non-structural carbohydrates and labile nitrogen. There are
two additional carbon pools, one for the storage of growth respiration
reserves, and another used to meet excess demand for maintenance
respiration during periods with low photosynthesis. One additional
nitrogen pool tracks retranslocated nitrogen, mobilized from leaf tissue
prior to abscission and litterfall. Altogether there are 20 state
variables for vegetation carbon, and 19 for vegetation nitrogen.

Figure 13.1. Vegetation fluxes and pools.

.. image:: image1.png
    :width: 800px
    :height: 800px

Figure 13.2: Carbon and nitrogen pools.

.. image:: image2.png
    :width: 400px
    :height: 400px

Carbon Allocation for Maintenance Respiration Costs
--------------------------------------------------------

The carbon and nitrogen allocation routines in CLM determine the fate of
newly assimilated carbon, coming from the calculation of photosynthesis,
and available mineral nitrogen, coming from plant uptake of mineral
nitrogen in the soil or being drawn out of plant reserves. Allocation
fluxes are determined in three steps: first :math:`CF_{GPPpot}` is
used to evaluate the potential allocation of carbon and nitrogen
assuming an unlimited nitrogen supply, then the actual nitrogen supply
is compared against the demand, and finally allocation of carbon and
nitrogen are reduced, if necessary, to match nitrogen supply and demand.

Allocation of available carbon on each time step is prioritized, with
first priority given to the demand for carbon to support maintenance
respiration of live tissues (section 13.7). Second priority is to
replenish the internal plant carbon pool that supports maintenance
respiration during times when maintenance respiration exceeds
photosynthesis (e.g. at night, during winter for perennial vegetation,
or during periods of drought stress) (Sprugel et al., 1995). Third
priority is to support growth of new tissues, including allocation to
storage pools from which new growth will be displayed in subsequent time
steps.

The total maintenance respiration demand (:math:`CF_{mr}`, gC
m\ :sup:`-2` s\ :sup:`-1`) is calculated as a function of
tissue mass and nitrogen concentration, and temperature (section 13.7).
The carbon supply to support this demand is composed of fluxes allocated
from carbon assimilated in the current timestep
(:math:`CF_{GPP,mr}`, gC m\ :sup:`-2` s\ :sup:`-1`)
and from a storage pool that is drawn down when total demand exceeds
photosynthesis ( :math:`CF_{xs,mr}`, gC m\ :sup:`-2`
s\ :sup:`-1`):

.. math::
   :label: 17.1) 

   CF_{mr} =CF_{GPP,mr} +CF_{xs,mr}

.. math::
   :label: 17.2) 

   CF_{GPP,mr} =\_ \left\{\begin{array}{l} {CF_{mr} \qquad \qquad {\rm for\; }CF_{mr} \le CF_{GPPpot} } \\ {CF_{GPPpot} \qquad {\rm for\; }CF_{mr} >CF_{GPPpot} } \end{array}\right.

.. math::
   :label: 17.3) 

   CF_{xs,mr} =\_ \left\{\begin{array}{l} {0\qquad \qquad \qquad {\rm for\; }CF_{mr} \le CF_{GPPpot} } \\ {CF_{mr} -CF_{GPPpot} \qquad {\rm for\; }CF_{mr} >CF_{GPPpot} } \end{array}\right.

The storage pool that supplies carbon for maintenance respiration in
excess of current  :math:`CF_{GPPpot}` ( :math:`CS_{xs}`, gC
m\ :sup:`-2`) is permitted to run a deficit (negative state), and
the magnitude of this deficit determines an allocation demand which
gradually replenishes  :math:`CS_{xs}`. The logic for allowing a
negative state for this pool is to eliminate the need to know in advance
what the total maintenance respiration demand will be for a particular
combination of climate and plant type. Using the deficit approach, the
allocation to alleviate the deficit increases as the deficit increases,
until the supply of carbon into the pool balances the demand for carbon
leaving the pool in a quasi-steady state, with variability driven by the
seasonal cycle, climate variation, disturbance, and internal dynamics of
the plant-litter-soil system. In cases where the combination of climate
and plant type are not suitable to sustained growth, the deficit in this
pool increases until the available carbon is being allocated mostly to
alleviate the deficit, and new growth approaches zero. The allocation
flux to  :math:`CS_{xs}` (:math:`CF_{GPP,xs}`, gC
m\ :sup:`-2` s\ :sup:`-1`) is given as

.. math::
   :label: 17.4) 

   CF_{GPP,xs,pot} =\left\{\begin{array}{l} {0\qquad \qquad \qquad {\rm for\; }CS_{xs} \ge 0} \\ {-CS_{xs} /(86400\tau _{xs} )\qquad {\rm for\; }CS_{xs} <0} \end{array}\right.

.. math::
   :label: 17.5) 

   CF_{GPP,xs} =\left\{\begin{array}{l} {CF_{GPP,xs,pot} \qquad \qquad \qquad {\rm for\; }CF_{GPP,xs,pot} \le CF_{GPPpot} -CF_{GPP,mr} } \\ {\max (CF_{GPPpot} -CF_{GPP,mr} ,0)\qquad {\rm for\; }CF_{GPP,xs,pot} >CF_{GPPpot} -CF_{GPP,mr} } \end{array}\right.

where :math:`\tau_{xs}` is the time constant (currently
set to 30 days) controlling the rate of replenishment of :math:`CS_{xs}`.

Note that these two top-priority carbon allocation fluxes
(:math:`CF_{GPP,mr}` and :math:`CF_{GPP,xs}`) are not
stoichiometrically associated with any nitrogen fluxes, and so this
initial allocation step can proceed without reference to (or limitation
from) the available mineral nitrogen supply.

Carbon and Nitrogen Stoichiometry of New Growth
----------------------------------------------------

After accounting for the carbon cost of maintenance respiration, the
remaining carbon flux from photosynthesis which can be allocated to new
growth (:math:`CF_{avail}`, gC m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 17.6) 

   CF_{avail\_ alloc} =CF_{GPPpot} -CF_{GPP,mr} -CF_{GPP,xs} .

Potential allocation to new growth is calculated for all of the plant
carbon and nitrogen state variables based on specified C:N ratios for
each tissue type and allometric parameters that relate allocation
between various tissue types. The allometric parameters are defined as
follows:

.. math::
   :label: ZEqnNum650137 

   \begin{array}{l} {a_{1} ={\rm \; ratio\; of\; new\; fine\; root\; :\; new\; leaf\; carbon\; allocation}} \\ {a_{2} ={\rm \; ratio\; of\; new\; coarse\; root\; :\; new\; stem\; carbon\; allocation}} \\ {a_{3} ={\rm \; ratio\; of\; new\; stem\; :\; new\; leaf\; carbon\; allocation}} \\ {a_{4} ={\rm \; ratio\; new\; live\; wood\; :\; new\; total\; wood\; allocation}} \\ {g_{1} ={\rm ratio\; of\; growth\; respiration\; carbon\; :\; new\; growth\; carbon.\; }} \end{array}

Parameters :math:`a_{1}`, :math:`a_{2}`, and :math:`a_{4}` are defined as constants for a given PFT (Table
13.1), while  :math:`g_{l }` = 0.3 (unitless) is prescribed as a
constant for all PFTs, based on construction costs for a range of woody
and non-woody tissues (Larcher, 1995).

The model includes a dynamic allocation scheme for woody vegetation
(parameter :math:`a_{3}` = -1, Table 13.1), in which case the
ratio for carbon allocation between new stem and new leaf increases with
increasing net primary production (NPP), as

.. math::
   :label: 17.8) 

   a_{3} =\frac{2.7}{1+e^{-0.004NPP_{ann} -300} } -0.4

where :math:`NPP_{ann}` is the annual sum of NPP from the previous
year. This mechanism has the effect of increasing woody allocation in
favorable growth environments (Allen et al., 2005; Vanninen and Makela,
2005) and during the phase of stand growth prior to canopy closure
(Axelsson and Axelsson, 1986).

Table 13.1. Allocation and carbon:nitrogen ratio parameters

+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Plant functional type            | :math:`a_{1}`         | :math:`a_{2}`         | :math:`a_{3}`         | :math:`a_{4}`         |  :math:`CN_{leaf}`        |  :math:`CN_{fr}`        | :math:`CN_{lw}`         | :math:`CN_{dw}`         |
+==================================+=======================+=======================+=======================+=======================+===========================+=========================+=========================+=========================+
| NET Temperate                    | 1                     | 0.3                   | -1                    | 0.1                   | 35                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| NET Boreal                       | 1                     | 0.3                   | -1                    | 0.1                   | 40                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| NDT Boreal                       | 1                     | 0.3                   | -1                    | 0.1                   | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BET Tropical                     | 1                     | 0.3                   | -1                    | 0.1                   | 30                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BET temperate                    | 1                     | 0.3                   | -1                    | 0.1                   | 30                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BDT tropical                     | 1                     | 0.3                   | -1                    | 0.1                   | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BDT temperate                    | 1                     | 0.3                   | -1                    | 0.1                   | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BDT boreal                       | 1                     | 0.3                   | -1                    | 0.1                   | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BES temperate                    | 1                     | 0.3                   | 0.2                   | 0.5                   | 30                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BDS temperate                    | 1                     | 0.3                   | 0.2                   | 0.5                   | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| BDS boreal                       | 1                     | 0.3                   | 0.2                   | 0.1                   | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| C\ :sub:`3` arctic grass         | 1                     | 0                     | 0                     | 0                     | 25                        | 42                      | 0                       | 0                       |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| C\ :sub:`3` grass                | 2                     | 0                     | 0                     | 0                     | 25                        | 42                      | 0                       | 0                       |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| C\ :sub:`4` grass                | 2                     | 0                     | 0                     | 0                     | 25                        | 42                      | 0                       | 0                       |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Crop R                           | 2                     | 0                     | 0                     | 0                     | 25                        | 42                      | 0                       | 0                       |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Crop I                           | 2                     | 0                     | 0                     | 0                     | 25                        | 42                      | 0                       | 0                       |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Corn R                           | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Corn I                           | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Temp Cereal R                    | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Temp Cereal I                    | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Winter Cereal R                  | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Winter Cereal I                  | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Soybean R                        | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
| Soybean I                        | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
+----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+

Carbon to nitrogen ratios are defined for different tissue types as
follows:

.. math::
   :label: ZEqnNum413927 

   \begin{array}{l} {CN_{leaf} =\_ {\rm \; C:N\; for\; leaf}} \\ {CN_{fr} =\_ {\rm \; C:N\; for\; fine\; root}} \\ {CN_{lw} =\_ {\rm \; C:N\; for\; live\; wood\; (in\; stem\; and\; coarse\; root)}} \\ {CN_{dw} =\_ {\rm \; C:N\; for\; dead\; wood\; (in\; stem\; and\; coarse\; root)}} \end{array}

where all C:N parameters are defined as constants for a given PFT (Table
13.1).

Given values for the parameters in and , total carbon and nitrogen
allocation to new growth ( :math:`CF_{alloc}`, gC
m\ :sup:`-2` s\ :sup:`-1`, and :math:`NF_{alloc}`, gN
m\ :sup:`-2` s\ :sup:`-1`, respectively) can be expressed as
functions of new leaf carbon allocation (:math:`CF_{GPP,leaf}`, gC
m\ :sup:`-2` s\ :sup:`-1`):

.. math::
   :label: ZEqnNum555154 

   \begin{array}{l} {CF_{alloc} =CF_{GPP,leaf} {\kern 1pt} C_{allom} } \\ {NF_{alloc} =CF_{GPP,leaf} {\kern 1pt} N_{allom} } \end{array}

where

.. math::
   :label: 17.11) 

   \begin{array}{l} {C_{allom} =\left\{\begin{array}{l} {\left(1+g_{1} \right)\left(1+a_{1} +a_{3} \left(1+a_{2} \right)\right)\qquad {\rm for\; woody\; PFT}} \\ {1+g_{1} +a_{1} \left(1+g_{1} \right)\qquad \qquad {\rm for\; non-woody\; PFT}} \end{array}\right. } \\ {} \end{array}

.. math::
   :label: 17.12) 

   N_{allom} =\left\{\begin{array}{l} {\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } +\frac{a_{3} a_{4} \left(1+a_{2} \right)}{CN_{lw} } +} \\ {\qquad \frac{a_{3} \left(1-a_{4} \right)\left(1+a_{2} \right)}{CN_{dw} } \qquad {\rm for\; woody\; PFT}} \\ {\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } \qquad \qquad \qquad {\rm for\; non-woody\; PFT.}} \end{array}\right.

Since the C:N stoichiometry for new growth allocation is defined, from
Eq. , as :math:`C_{allom}`/ :math:`N_{allom}`, the total carbon available for new growth allocation
(:math:`CF_{avail\_alloc}`) can be used to calculate the total
plant nitrogen demand for new growth ( :math:`NF_{plant\_demand}`,
gN m\ :sup:`-2` s\ :sup:`-1`) as:

.. math::
   :label: 17.13) 

   NF_{plant\_ demand} =CF_{avail\_ alloc} \frac{N_{allom} }{C_{allom} } .

Deployment of retranslocated nitrogen
------------------------------------------

In many plants, some portion of the nitrogen used to construct new
tissues is mobilized from senescing tissues, especially leaves, and
retained within the plant when the tissues are lost as litter. This
store of retranslocated nitrogen is used to supply part of the nitrogen
demand for subsequent growth (Magill et al., 1997; Oikawa et al., 2005;
Son and Gower, 1991). CLM includes one pool of retranslocated nitrogen
(:math:`NS_{retrans}`, gN m\ :sup:`-2`), and the
availability of nitrogen from this pool to support new growth
(:math:`NF_{avail\_retrans}`, gN m\ :sup:`-2`
s\ :sup:`-1`) is proportional to the plant nitrogen demand, as:

.. math::
   :label: 17.14) 

   NF_{avail\_ retrans} =\min \left(\frac{NF_{retrans\_ ann} \frac{NF_{plant\_ demand} }{NF_{plant\_ demand\_ ann} } }{\Delta t} ,\; \frac{NS_{retrans} }{\Delta t} \right)

where :math:`NF_{retrans\_ann}` (gN m\ :sup:`-2` y\ :sup:`-1`) is the previous year’s annual sum of retranslocated
nitrogen extracted from senescing tissues,
:math:`NF_{plant\_demand\_ann}` (gN m\ :sup:`-2` y\ :sup:`-1`) is the previous year’s annual sum of
:math:`NF_{plant\_demand}`, and :math:`\Delta`\ *t* (s) is the
model’s biogeochemistry time step. This formulation produces an annual
cycle in the extraction of nitrogen from :math:`NS_{retrans}`
which corresponds to the annual cycle of plant nitrogen demand, and
which is scaled to give :math:`NS_{retrans}` approximately a
one-year turnover time. The minimum function prevents extraction of more
than the remaining pool of retranslocated nitrogen, which can be an
important constraint under conditions where high rates of mortality are
modifying the size of the pool. During the first year of an initial
simulation, before :math:`NF_{plant\_demand\_ann}` and
:math:`NF_{retrans\_ann}` have valid values,
:math:`NF_{avail\_retrans}` is set to 0.0.

The actual flux of nitrogen from the retranslocated N pool into
allocation of new growth (:math:`NF_{retrans,alloc}`, gN
m\ :sup:`-2` s\ :sup:`-1`) is never greater than the plant
demand for new nitrogen:

.. math::
   :label: 17.15) 

   NF_{retrans,alloc} =\min \left(NF_{plant\_ demand} ,NF_{avail\_ retrans} \right)

Plant nitrogen uptake from soil mineral nitrogen pool
----------------------------------------------------------

The total plant nitrogen demand is reduced by the nitrogen flux from
:math:`NS_{retrans}` to give the plant demand for mineral nitrogen
from the soil (:math:`NF_{plant\_demand\_soil}`, gN
m\ :sup:`-2` s\ :sup:`-1`):

.. math::
   :label: ZEqnNum491412 

   NF_{plant\_ demand\_ soil} =NF_{plant\_ demand} -NF_{retrans,alloc} .

The combined demand from all PFTs sharing space on a soil column and the
demand from the heterotrophic community in the soil (nitrogen
immobilization demand) compete for the available soil mineral nitrogen
pool. The result of this competition is passed back to the allocation
algorithm as :math:`f_{plant\_demand}`, the fraction (from 0 to 1)
of the plant nitrogen demand which can be met given the current soil
mineral nitrogen supply and competition with heterotrophs. Plant uptake
from the soil mineral nitrogen pool is then given as:

.. math::
   :label: 17.17) 

   NF_{sminn,alloc} =NF_{plant\_ demand\_ soil} f_{plant\_ demand}

Final carbon and nitrogen allocation
-----------------------------------------

The total flux of allocated nitrogen is given as:

.. math::
   :label: 17.18) 

   NF_{alloc} =NF_{retrans,alloc} +NF_{sminn,alloc}

From the stoichiometric relationship in Eq. , the associated carbon
allocation flux is:

.. math::
   :label: 17.19) 

   CF_{alloc} =NF_{alloc} \frac{C_{allom} }{N_{allom} } .

The downregulation of photosynthesis can then be calculated as:

.. math::
   :label: 17.20) 

   f_{dreg} =\frac{CF_{alloc} -CF_{avail\_ alloc} }{CF_{GPPpot} } .

Total allocation to new leaf carbon
(:math:`CF_{alloc,leaf\_tot}`, gC m\ :sup:`-2` s\ :sup:`-1`) is calculated as:

.. math::
   :label: 17.21) 

   CF_{alloc,leaf\_ tot} =\frac{CF_{alloc} }{C_{allom} } .

There are two carbon pools associated with each plant tissue – one which
represents the currently displayed tissue, and another which represents
carbon stored for display in a subsequent growth period. The nitrogen
pools follow this same organization. The model keeps track of stored
carbon according to which tissue type it will eventually be displayed
as, and the separation between display in the current timestep and
storage for later display depends on the parameter :math:`f_{cur}`
(values 0 to 1). Given :math:`CF_{alloc,leaf}` and :math:`f_{cur}`, the allocation fluxes of carbon to display and
storage pools (where storage is indicated with *\_stor*) for the various
tissue types are given as:

.. math::
   :label: 17.22) 

   CF_{alloc,leaf} \_ =CF_{alloc,leaf\_ tot} f_{cur}

.. math::
   :label: 17.23) 

   CF_{alloc,leaf\_ stor} \_ =CF_{alloc,leaf\_ tot} \left(1-f_{cur} \right)

.. math::
   :label: 17.24) 

   CF_{alloc,froot} \_ =CF_{alloc,leaf\_ tot} a_{1} f_{cur}

.. math::
   :label: 17.25) 

   CF_{alloc,froot\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{1} \left(1-f_{cur} \right)

.. math::
   :label: 17.26) 

   CF_{alloc,livestem} \_ =CF_{alloc,leaf\_ tot} a_{3} a_{4} f_{cur}

.. math::
   :label: 17.27) 

   CF_{alloc,livestem\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{3} a_{4} \left(1-f_{cur} \right)

.. math::
   :label: 17.28) 

   CF_{alloc,deadstem} \_ =CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)f_{cur}

.. math::
   :label: 17.29) 

   CF_{alloc,deadstem\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)\left(1-f_{cur} \right)

.. math::
   :label: 17.30) 

   CF_{alloc,livecroot} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} f_{cur}

.. math::
   :label: 17.31) 

   CF_{alloc,livecroot\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} \left(1-f_{cur} \right)

.. math::
   :label: 17.32) 

   CF_{alloc,deadcroot} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} \left(1-a_{4} \right)f_{cur}

.. math::
   :label: 17.33) 

   CF_{alloc,deadcroot\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} \left(1-a_{4} \right)\left(1-f_{cur} \right).

The corresponding nitrogen allocation fluxes are given as:

.. math::
   :label: 17.34) 

   NF_{alloc,leaf} \_ =\frac{CF_{alloc,leaf\_ tot} }{CN_{leaf} } f_{cur}

.. math::
   :label: 17.35) 

   NF_{alloc,leaf\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} }{CN_{leaf} } \left(1-f_{cur} \right)

.. math::
   :label: 17.36) 

   NF_{alloc,froot} \_ =\frac{CF_{alloc,leaf\_ tot} a_{1} }{CN_{fr} } f_{cur}

.. math::
   :label: 17.37) 

   NF_{alloc,froot\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{1} }{CN_{fr} } \left(1-f_{cur} \right)

.. math::
   :label: 17.38) 

   NF_{alloc,livestem} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} a_{4} }{CN_{lw} } f_{cur}

.. math::
   :label: 17.39) 

   NF_{alloc,livestem\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} a_{4} }{CN_{lw} } \left(1-f_{cur} \right)

.. math::
   :label: 17.40) 

   NF_{alloc,deadstem} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)}{CN_{dw} } f_{cur}

.. math::
   :label: 17.41) 

   NF_{alloc,deadstem\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)}{CN_{dw} } \left(1-f_{cur} \right)

.. math::
   :label: 17.42) 

   NF_{alloc,livecroot} \_ =\frac{CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} }{CN_{lw} } f_{cur}

.. math::
   :label: 17.43) 

   NF_{alloc,livecroot\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} }{CN_{lw} } \left(1-f_{cur} \right)

.. math::
   :label: 17.44) 

   NF_{alloc,deadcroot} \_ =\frac{CF_{alloc,leaf\_ tot} a_{2} a_{3} \left(1-a_{4} \right)}{CN_{dw} } f_{cur}

.. math::
   :label: 17.45) 

   NF_{alloc,deadcroot\_ stor} \_ =\frac{CF_{alloc,leaf} a_{2} a_{3} \left(1-a_{4} \right)}{CN_{dw} } \left(1-f_{cur} \right).

Autotrophic Respiration
----------------------------

The model treats maintenance and growth respiration fluxes separately,
even though it is difficult to measure them as separate fluxes (Lavigne
and Ryan, 1997; Sprugel et al., 1995). Maintenance respiration is
defined as the carbon cost to support the metabolic activity of existing
live tissue, while growth respiration is defined as the additional
carbon cost for the synthesis of new growth.

13.7.1 Maintenance Respiration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Under the assumption that tissue nitrogen content is a suitable index of
cellular metabolic rate, maintenance respiration costs for live tissues
(leaf, live stem, live coarse root, and fine root) are calculated as
functions tissue nitrogen content and the relevant temperature,
following the empirical relationship reported by Ryan (1991):

.. math::
   :label: 17.46) 

   CF_{mr\_ leaf} \_ =NS_{leaf} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

.. math::
   :label: 17.47) 

   CF_{mr\_ livestem} \_ =NS_{livestem} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

.. math::
   :label: 17.48) 

   CF_{mr\_ livecroot} \_ =NS_{livecroot} MR_{base} MR_{Q10} ^{(T_{2m} -20)/10}

.. math::
   :label: 17.49) 

   CF_{mr\_ froot} \_ =\sum _{j=1}^{nlevsoi}NS_{froot} rootfr_{j} MR_{base} MR_{Q10} ^{(Ts_{j} -20)/10}

where :math:`MR_{base}` (= 2.525e\ :sup:`-6` gC gN\ :sup:`-1` s\ :sup:`-1`) is the base rate of maintenance
respiration per unit nitrogen content, :math:`MR_{q10}` (= 2.0) is
the temperature sensitivity for maintenance respiration,
:math:`T_{2m}` (:sup:`o`\ C) is the air temperature at 2m
height, :math:`Ts_{j}`* (:sup:`o`\ C) is the soil
temperature at level *j*, and :math:`rootfr_{j}` is the fraction
of fine roots distributed in soil level *j*.

Note that, for woody vegetation, maintenance respiration costs are not
calculated for the dead stem and dead coarse root components. These
components are assumed to consist of dead xylem cells, with no metabolic
function. By separating the small live component of the woody tissue
(ray parenchyma, phloem, and sheathing lateral meristem cells) from the
larger fraction of dead woody tissue, it is reasonable to assume a
common base maintenance respiration rate for all live tissue types.

The total maintenance respiration cost is then given as:

.. math::
   :label: 17.50) 

   CF_{mr} =CF_{mr\_ leaf} +CF_{mr\_ froot} +CF_{mr\_ livestem} +CF_{mr\_ livecroot} .

13.7.2 Growth Respiration
^^^^^^^^^^^^^^^^^^^^^^^^^

Growth respiration is calculated as a factor of 0.3 times the total
carbon in new growth on a given timestep, based on construction costs
for a range of woody and non-woody tissues (Larcher, 1995). For new
carbon and nitrogen allocation that enters storage pools for subsequent
display, it is not clear what fraction of the associated growth
respiration should occur at the time of initial allocation, and what
fraction should occur later, at the time of display of new growth from
storage. Eddy covariance estimates of carbon fluxes in forest ecosystems
suggest that the growth respiration associated with transfer of
allocated carbon and nitrogen from storage into displayed tissue is not
significant (Churkina et al., 2003), and so it is assumed in CLM that
all of the growth respiration cost is incurred at the time of initial
allocation, regardless of the fraction of allocation that is displayed
immediately (i.e. regardless of the value of :math:`f_{cur}`,
section 13.5). This behavior is parameterized in such a way that if
future research suggests that some fraction of the growth respiration
cost should be incurred at the time of display from storage, a simple
parameter modification will effect the change. [1]_

.. [1]
   Parameter :math:`\text{grpnow}`  in routines CNGResp and  CNAllocation, currently set to 1.0, could be changed to a smaller
   value to transfer some portion (1 - :math:`\text{grpnow}` ) of the growth respiration forward in time to occur at the time of growth
   display from storage.

