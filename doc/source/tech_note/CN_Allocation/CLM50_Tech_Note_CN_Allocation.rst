.. _rst_CN Allocation:

Carbon and Nitrogen Allocation
==============================

Introduction
-----------------

The carbon and nitrogen allocation routines in CLM determine the fate of newly assimilated carbon, coming from the calculation of photosynthesis, and available mineral nitrogen, coming from plant uptake of mineral nitrogen in the soil or being drawn out of plant reserves. A significant change to CLM5 relative to prior versions is that allocation of carbon and nitrogen proceed independently rather than in a sequential manner.

Carbon Allocation for Maintenance Respiration Costs
--------------------------------------------------------

Allocation of available carbon on each time step is prioritized, with first priority given to the demand for carbon to support maintenance respiration of live tissues (section 13.7). Second priority is to replenish the internal plant carbon pool that supports maintenance respiration during times when maintenance respiration exceeds photosynthesis (e.g. at night, during winter for perennial vegetation, or during periods of drought stress) (Sprugel et al., 1995). Third priority is to support growth of new tissues, including allocation to storage pools from which new growth will be displayed in subsequent time steps.

The total maintenance respiration demand (:math:`CF_{mr}`, gC m\ :sup:`-2` s\ :sup:`-1`) is calculated as a function of tissue mass and nitrogen concentration, and temperature (section 13.7) The carbon supply to support this demand is composed of fluxes allocated from carbon assimilated in the current timestep (:math:`CF_{GPP,mr}`, gC m\ :sup:`-2` s\ :sup:`-1` and from a storage pool that is drawn down when total demand exceeds photosynthesis ( :math:`CF_{xs,mr}`, gC m\ :sup:`-2` s\ :sup:`-1`):

.. math::
   :label: 19.1

   CF_{mr} =CF_{GPP,mr} +CF_{xs,mr}

.. math::
   :label: 19.2

   CF_{GPP,mr} =\_ \left\{\begin{array}{l} {CF_{mr} \qquad \qquad {\rm for\; }CF_{mr} \le CF_{GPP} } \\ {CF_{GPP} \qquad {\rm for\; }CF_{mr} >CF_{GPP} } \end{array}\right.

.. math::
   :label: 19.3

   CF_{xs,mr} =\_ \left\{\begin{array}{l} {0\qquad \qquad \qquad {\rm for\; }CF_{mr} \le CF_{GPP} } \\ {CF_{mr} -CF_{GPP} \qquad {\rm for\; }CF_{mr} >CF_{GPP} } \end{array}\right.

The storage pool that supplies carbon for maintenance respiration in excess of current :math:`CF_{GPP}` ( :math:`CS_{xs}`, gC m\ :sup:`-2`) is permitted to run a deficit (negative state), and the magnitude of this deficit determines an allocation demand which gradually replenishes :math:`CS_{xs}`. The logic for allowing a negative state for this pool is to eliminate the need to know in advance what the total maintenance respiration demand will be for a particular combination of climate and plant type. Using the deficit approach, the allocation to alleviate the deficit increases as the deficit increases, until the supply of carbon into the pool balances the demand for carbon leaving the pool in a quasi-steady state, with variability driven by the seasonal cycle, climate variation, disturbance, and internal dynamics of the plant-litter-soil system. In cases where the combination of climate and plant type are not suitable to sustained growth, the deficit in this pool increases until the available carbon is being allocated mostly to alleviate the deficit, and new growth approaches zero. The allocation flux to :math:`CS_{xs}` (:math:`CF_{GPP,xs}`, gC m\ :sup:`-2` s\ :sup:`-1`) is given as

.. math::
   :label: 19.4

   CF_{GPP,xs,pot} =\left\{\begin{array}{l} {0\qquad \qquad \qquad {\rm for\; }CS_{xs} \ge 0} \\ {-CS_{xs} /(86400\tau _{xs} )\qquad {\rm for\; }CS_{xs} <0} \end{array}\right.

.. math::
   :label: 19.5

   CF_{GPP,xs} =\left\{\begin{array}{l} {CF_{GPP,xs,pot} \qquad \qquad \qquad {\rm for\; }CF_{GPP,xs,pot} \le CF_{GPP} -CF_{GPP,mr} } \\ {\max (CF_{GPP} -CF_{GPP,mr} ,0)\qquad {\rm for\; }CF_{GPP,xs,pot} >CF_{GPP} -CF_{GPP,mr} } \end{array}\right.

where :math:`\tau_{xs}` is the time constant (currently set to 30 days) controlling the rate of replenishment of :math:`CS_{xs}`.

Note that these two top-priority carbon allocation fluxes (:math:`CF_{GPP,mr}` and :math:`CF_{GPP,xs}`) are not stoichiometrically associated with any nitrogen fluxes.

Carbon and Nitrogen Stoichiometry of New Growth
----------------------------------------------------

After accounting for the carbon cost of maintenance respiration, the remaining carbon flux from photosynthesis which can be allocated to new growth (:math:`CF_{avail}`, gC m\ :sup:`-2` s\ :sup:`-1`) is

.. math::
   :label: 19.6

   CF_{avail\_ alloc} =CF_{GPP} -CF_{GPP,mr} -CF_{GPP,xs} .

Potential allocation to new growth is calculated for all of the plant carbon and nitrogen state variables based on specified C:N ratios for each tissue type and allometric parameters that relate allocation between various tissue types. The allometric parameters are defined as follows:

.. math::
   :label: 19.7

   \begin{array}{l} {a_{1} ={\rm \; ratio\; of\; new\; fine\; root\; :\; new\; leaf\; carbon\; allocation}} \\ {a_{2} ={\rm \; ratio\; of\; new\; coarse\; root\; :\; new\; stem\; carbon\; allocation}} \\ {a_{3} ={\rm \; ratio\; of\; new\; stem\; :\; new\; leaf\; carbon\; allocation}} \\ {a_{4} ={\rm \; ratio\; new\; live\; wood\; :\; new\; total\; wood\; allocation}} \\ {g_{1} ={\rm ratio\; of\; growth\; respiration\; carbon\; :\; new\; growth\; carbon.\; }} \end{array}

Parameters :math:`a_{1}`, :math:`a_{2}`, and :math:`a_{4}` are defined as constants for a given PFT (Table 13.1), while :math:`g_{l }` = 0.3 (unitless) is prescribed as a constant for all PFTs, based on construction costs for a range of woody and non-woody tissues (Larcher, 1995).

The model includes a dynamic allocation scheme for woody vegetation (parameter :math:`a_{3}` = -1, :numref:`Table Allocation and CN ratio parameters`), in which case the ratio for carbon allocation between new stem and new leaf increases with increasing net primary production (NPP), as

.. math::
   :label: 19.8

   a_{3} =\frac{2.7}{1+e^{-0.004NPP_{ann} -300} } -0.4

where :math:`NPP_{ann}` is the annual sum of NPP from the previous year. This mechanism has the effect of increasing woody allocation in favorable growth environments (Allen et al., 2005; Vanninen and Makela, 2005) and during the phase of stand growth prior to canopy closure (Axelsson and Axelsson, 1986).

.. _Table Allocation and CN ratio parameters:

.. table:: Allocation and target carbon\:nitrogen ratio parameters

 +----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
 | Plant functional type            | :math:`a_{1}`         | :math:`a_{2}`         | :math:`a_{3}`         | :math:`a_{4}`         |  :math:`Target CN_{leaf}` |  :math:`Target CN_{fr}` | :math:`Target CN_{lw}`  | :math:`Target CN_{dw}`  |
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
 | Miscanthus R                     | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
 +----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
 | Miscanthus I                     | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
 +----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
 | Switchgrass R                    | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
 +----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+
 | Switchgrass I                    | 2                     | 0                     | 0                     | 1                     | 25                        | 42                      | 50                      | 500                     |
 +----------------------------------+-----------------------+-----------------------+-----------------------+-----------------------+---------------------------+-------------------------+-------------------------+-------------------------+

Carbon to nitrogen ratios are defined for different tissue types as follows:

.. math::
   :label: 19.9

   \begin{array}{l} {CN_{leaf} =\_ {\rm \; C:N\; for\; leaf}} \\ {CN_{fr} =\_ {\rm \; C:N\; for\; fine\; root}} \\ {CN_{lw} =\_ {\rm \; C:N\; for\; live\; wood\; (in\; stem\; and\; coarse\; root)}} \\ {CN_{dw} =\_ {\rm \; C:N\; for\; dead\; wood\; (in\; stem\; and\; coarse\; root)}} \end{array}

where all C:N parameters are defined as constants for a given PFT (:numref:`Table Allocation and CN ratio parameters`).

Given values for the parameters in and, total carbon and nitrogen allocation to new growth ( :math:`CF_{alloc}`, gC m\ :sup:`-2` s\ :sup:`-1`, and :math:`NF_{alloc}`, gN m\ :sup:`-2` s\ :sup:`-1`, respectively) can be expressed as functions of new leaf carbon allocation (:math:`CF_{GPP,leaf}`, gC m\ :sup:`-2` s\ :sup:`-1`):

.. math::
   :label: 19.10

   \begin{array}{l} {CF_{alloc} =CF_{GPP,leaf} {\kern 1pt} C_{allom} } \\ {NF_{alloc} =CF_{GPP,leaf} {\kern 1pt} N_{allom} } \end{array}

where

.. math::
   :label: 19.11

   \begin{array}{l} {C_{allom} =\left\{\begin{array}{l} {\left(1+g_{1} \right)\left(1+a_{1} +a_{3} \left(1+a_{2} \right)\right)\qquad {\rm for\; woody\; PFT}} \\ {1+g_{1} +a_{1} \left(1+g_{1} \right)\qquad \qquad {\rm for\; non-woody\; PFT}} \end{array}\right. } \\ {} \end{array}

.. math::
   :label: 19.12

   N_{allom} =\left\{\begin{array}{l} {\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } +\frac{a_{3} a_{4} \left(1+a_{2} \right)}{CN_{lw} } +} \\ {\qquad \frac{a_{3} \left(1-a_{4} \right)\left(1+a_{2} \right)}{CN_{dw} } \qquad {\rm for\; woody\; PFT}} \\ {\frac{1}{CN_{leaf} } +\frac{a_{1} }{CN_{fr} } \qquad \qquad \qquad {\rm for\; non-woody\; PFT.}} \end{array}\right.

Since the C:N stoichiometry for new growth allocation is defined, from Eq., as :math:`C_{allom}`/ :math:`N_{allom}`, the total carbon available for new growth allocation (:math:`CF_{avail\_alloc}`) can be used to calculate the total plant nitrogen demand for new growth ( :math:`NF_{plant\_demand}`, gN m\ :sup:`-2` s\ :sup:`-1`) as:

.. math::
   :label: 19.13

   NF_{plant\_ demand} =CF_{avail\_ alloc} \frac{N_{allom} }{C_{allom} } .

.. _Carbon Allocation to New Growth:

Carbon Allocation to New Growth
-----------------------------------------

There are two carbon pools associated with each plant tissue – one which represents the currently displayed tissue, and another which represents carbon stored for display in a subsequent growth period. The nitrogen pools follow this same organization. The model keeps track of stored carbon according to which tissue type it will eventually be displayed as, and the separation between display in the current timestep and storage for later display depends on the parameter :math:`f_{cur}` (values 0 to 1). Given :math:`CF_{alloc,leaf}` and :math:`f_{cur}`, the allocation fluxes of carbon to display and storage pools (where storage is indicated with *\_stor*) for the various tissue types are given as:

.. math::
   :label: 19.14

   CF_{alloc,leaf} \_ =CF_{alloc,leaf\_ tot} f_{cur}

.. math::
   :label: 19.15

   CF_{alloc,leaf\_ stor} \_ =CF_{alloc,leaf\_ tot} \left(1-f_{cur} \right)

.. math::
   :label: 19.16

   CF_{alloc,froot} \_ =CF_{alloc,leaf\_ tot} a_{1} f_{cur}

.. math::
   :label: 19.17

   CF_{alloc,froot\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{1} \left(1-f_{cur} \right)

.. math::
   :label: 19.18

   CF_{alloc,livestem} \_ =CF_{alloc,leaf\_ tot} a_{3} a_{4} f_{cur}

.. math::
   :label: 19.19

   CF_{alloc,livestem\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{3} a_{4} \left(1-f_{cur} \right)

.. math::
   :label: 19.20

   CF_{alloc,deadstem} \_ =CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)f_{cur}

.. math::
   :label: 19.21

   CF_{alloc,deadstem\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)\left(1-f_{cur} \right)

.. math::
   :label: 19.22

   CF_{alloc,livecroot} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} f_{cur}

.. math::
   :label: 19.23

   CF_{alloc,livecroot\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} \left(1-f_{cur} \right)

.. math::
   :label: 19.24

   CF_{alloc,deadcroot} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} \left(1-a_{4} \right)f_{cur}

.. math::
   :label: 19.25

   CF_{alloc,deadcroot\_ stor} \_ =CF_{alloc,leaf\_ tot} a_{2} a_{3} \left(1-a_{4} \right)\left(1-f_{cur} \right).

Nitrogen allocation
-----------------------------------------

The total flux of nitrogen to be allocated is given by the FUN model (Chapter :numref:`rst_FUN`). This gives a total N to be allocated within a given timestep, :math:`N_{supply}`. The total N allocated for a given tissue :math:`i` is the minimum between the supply and the demand:

.. math::
   :label: 19.26

   NF_{alloc,i} = min \left( NF_{demand, i}, NF_{supply, i} \right)

The demand for each tissue, calculated for the tissue to remain on stoichiometry during growth, is:

.. math::
   :label: 19.27

   NF_{demand,leaf} \_ =\frac{CF_{alloc,leaf\_ tot} }{CN_{leaf} } f_{cur}

.. math::
   :label: 19.28

   NF_{demand,leaf\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} }{CN_{leaf} } \left(1-f_{cur} \right)

.. math::
   :label: 19.29

   NF_{demand,froot} \_ =\frac{CF_{alloc,leaf\_ tot} a_{1} }{CN_{fr} } f_{cur}

.. math::
   :label: 19.30

   NF_{demand,froot\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{1} }{CN_{fr} } \left(1-f_{cur} \right)

.. math::
   :label: 19.31

   NF_{demand,livestem} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} a_{4} }{CN_{lw} } f_{cur}

.. math::
   :label: 19.32

   NF_{demand,livestem\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} a_{4} }{CN_{lw} } \left(1-f_{cur} \right)

.. math::
   :label: 19.33

   NF_{demand,deadstem} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)}{CN_{dw} } f_{cur}

.. math::
   :label: 19.34

   NF_{demand,deadstem\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{3} \left(1-a_{4} \right)}{CN_{dw} } \left(1-f_{cur} \right)

.. math::
   :label: 19.35

   NF_{demand,livecroot} \_ =\frac{CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} }{CN_{lw} } f_{cur}

.. math::
   :label: 19.36

   NF_{demand,livecroot\_ stor} \_ =\frac{CF_{alloc,leaf\_ tot} a_{2} a_{3} a_{4} }{CN_{lw} } \left(1-f_{cur} \right)

.. math::
   :label: 19.37

   NF_{demand,deadcroot} \_ =\frac{CF_{alloc,leaf\_ tot} a_{2} a_{3} \left(1-a_{4} \right)}{CN_{dw} } f_{cur}

.. math::
   :label: 19.38

   NF_{demand,deadcroot\_ stor} \_ =\frac{CF_{alloc,leaf} a_{2} a_{3} \left(1-a_{4} \right)}{CN_{dw} } \left(1-f_{cur} \right).

After each pool's demand is calculated, the total plant N demand is then the sum of each individual pool :math:`i` corresponding to each tissue:

.. math::
   :label: 19.39

   NF_{demand,tot} = \sum _{i=tissues} NF_{demand,i}

and the total supply for each tissue :math:`i` is the product of the fractional demand and the total available N, calculated as the term :math:`N_{uptake}` equal to the sum of the eight N uptake streams described in the FUN model (Chapter :numref:`rst_FUN`).

.. math::
   :label: 19.40

   NF_{alloc,i} = N_{uptake} NF_{demand,i} / NF_{demand,tot}
