.. _rst_Crops and Irrigation:

Crops and Irrigation
========================

.. _Summary of CLM4.5 updates relative to the CLM4.0:

Summary of CLM4.5 updates relative to the CLM4.0
-----------------------------------------------------

We describe here the complete crop and irrigation parameterizations that
appear in CLM4.5. Corresponding information for CLM4.0 appeared on the
CLM4.0 web site in a pdf document independent of the CLM4.0 Technical
Note (Oleson et al. 2010a). The CLM4.0 crop model description also
appeared in Levis et al. (2012).

CLM4.5 includes the following updates to the CROP option, where CROP
refers to the interactive crop management model:

- Interactive irrigation

- Interactive fertilization

- Biological nitrogen fixation for soybeans

- Modified C:N ratios for crops

- Nitrogen retranslocation for crops

- Separate reproductive pool

These updates appear in detail in section 20.4. Most also appear in
Drewniak et al. (2013).

.. _The crop model:

The crop model
-------------------

Introduction
^^^^^^^^^^^^^^^^^^^

Groups developing Earth System Models generally account for the human
footprint on the landscape in simulations of historical and future
climates. Traditionally we have represented this footprint with natural
vegetation types and particularly grasses because they resemble many
common crops. Most modeling efforts have not incorporated more explicit
representations of land management such as crop type, planting,
harvesting, tillage, fertilization, and irrigation, because global scale
datasets of these factors have lagged behind vegetation mapping. As this
begins to change, we increasingly find models that will simulate the
biogeophysical and biogeochemical effects not only of natural but also
human-managed land cover.

AgroIBIS is a state-of-the-art land surface model with options to
simulate dynamic vegetation (Kucharik et al. 2000) and interactive
crop management (Kucharik and Brye 2003). The interactive crop
management parameterizations from AgroIBIS (March 2003 version) were
coupled as a proof-of-concept to the Community Land Model version 3
[CLM3.0, Oleson et al. (2004)] (not published), then coupled to the
CLM3.5 (Levis et al. 2009) and later released to the community with
CLM4CN (Levis et al. 2012).

With interactive crop management and, therefore, a more accurate
representation of agricultural landscapes, we hope to improve the CLM’s
simulated biogeophysics and biogeochemistry. These advances may improve
fully coupled simulations with the Community Earth System Model (CESM),
while helping human societies answer questions about changing food,
energy, and water resources in response to climate, environmental, land
use, and land management change (e.g., Kucharik and Brye 2003; Lobell et al. 2006).

.. _Crop plant functional types:

Crop plant functional types
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CLM’s default list of plant functional types (pfts) includes an
unmanaged crop (Table 2.1) treated as a second C3 grass. The unmanaged
crop has grid cell coverage assigned from satellite data, as do all
natural pfts when CLM’s dynamic vegetation model (CNDV; Castillo et al. 2012) is not active.

The new crop pfts used in the CLM get grid cell coverage from the
present-day crop dataset of Portmann et al. (2010). We assign these
managed crops in the proportions given by Portmann et al. (2010) without
exceeding the area previously assigned to the unmanaged crop. The
unmanaged crop continues to occupy any of its original area that remains
and continues to be handled just by the carbon/nitrogen cycling part of
the CLM (i.e., CN). The managed crop types (corn, soybean, and temperate
cereals) were chosen based on the availability of corresponding
algorithms in AgroIBIS. Temperate cereals include wheat, barley, and rye
here. We treat all temperate cereals as summer crops (like spring wheat,
for example) at this time. We may introduce winter cereals (such as
winter wheat) in a future version of the model.

To allow crops to coexist with natural vegetation in a grid cell and be
treated by separate models (i.e., CLM4.5CNcrop versus CLM4.5CNDV), we
separate the vegetated land unit into a naturally vegetated land unit
and a human managed land unit. Plant functional types in the naturally
vegetated land unit share one soil column and compete for water (default
CLM setting). Managed crop PFTs in the human managed land unit do not
share soil columns and thus permit for differences in land management
between crops.

.. _Phenology:

Phenology
^^^^^^^^^^^^^^^^

CLM4.5CN includes evergreen, seasonally deciduous (responding to changes
in day length), and stress deciduous (responding to changes in
temperature and/or soil moisture) phenology algorithms (Chapter 14). In
CLM4.5CNcrop we have added the AgroIBIS crop phenology algorithm,
consisting of three distinct phases.

Phase 1 starts at planting and ends with leaf emergence, phase 2
continues from leaf emergence to the beginning of grain fill, and phase
3 starts from the beginning of grain fill and ends with physiological
maturity and harvest.

.. _Planting:

Planting
'''''''''''''''''

Corn and temperate cereals must meet the following requirements between
April 1\ :sup:`st` and June 14\ :sup:`th` for planting in the northern hemisphere (NH):

.. math::
   :label: ZEqnNum568682 

   \begin{array}{l} 
   {T_{10d} >T_{p} } \\ 
   {T_{10d}^{\min } >T_{p}^{\min } }  \\ 
   {GDD_{8} \ge GDD_{\min } } 
   \end{array}

where :math:`{T}_{10d}` is the 10-day running mean of :math:`{T}_{2m}`, (the simulated 2-m air
temperature at every model time step) and :math:`T_{10d}^{\min}`  is
the 10-day running mean of :math:`T_{2m}^{\min }`  (the daily minimum of
:math:`{T}_{2m}`. :math:`{T}_{p}` and :math:`T_{p}^{\min }`  are crop-specific coldest planting temperatures
(Table 20.1), :math:`{GDD}_{8}` is the 20-year running mean growing
degree-days (units are degree-days or :sup:`o` days) tracked
from April through September (NH) base 8\ :sup:`o` C with
maximum daily increments of 30\ :sup:`o` days (see Eq.XXX ), and
:math:`{GDD}_{min }`\ is the minimum growing degree day requirement
(Table 20.1). Soy must meet the same requirements but between May
1\ :sup:`st` and June 14\ :sup:`th` for planting. If the
requirements in Eq. are not met by June 14\ :sup:`th`, then corn,
soybean, and temperate cereals are still planted on June
15\ :sup:`th` as long as  :math:`{GDD}_{8} > 0`. In
the southern hemisphere (SH) the NH requirements apply 6 months later.

:math:`{GDD}_{8}` does not change as quickly as :math:`{T}_{10d}` and :math:`T_{10d}^{\min }`, so
it determines whether the crop can be planted in a grid cell, while the
two faster-changing variables determine when the crop may be planted.

At planting, each crop is assigned 1 g leaf C m\ :sup:`-2` pft
column area to be transferred to the leaves upon leaf emergence. An
equivalent amount of seed leaf N is assigned given the pft’s C to N
ratio for leaves (:math:`{CN}_{leaf}`). (This differs from AgroIBIS,
which uses a seed leaf area index instead of seed C.)

At planting, the model updates the average growing degree-days necessary
for the crop to reach vegetative and physiological maturity,
:math:`{GDD}_{mat}`, according to the following AgroIBIS rules:

.. math::
   :label: 25.2) 

   \begin{array}{l} {GDD_{{\rm mat}}^{{\rm corn}} =0.85GDD_{{\rm 8}} {\rm \; \; \; and\; \; \; 950}<GDD_{{\rm mat}}^{{\rm corn}} <1850{}^\circ {\rm days}} \\ {GDD_{{\rm mat}}^{{\rm temp.\; cereals}} =GDD_{{\rm 0}} {\rm \; \; \; and\; \; \; }GDD_{{\rm mat}}^{{\rm temp.\; cereals}} <1700{}^\circ {\rm days}} \\ {GDD_{{\rm mat}}^{{\rm soy}} =GDD_{{\rm 10}} {\rm \; \; \; and\; \; \; }GDD_{{\rm mat}}^{{\rm soy}} <1700{}^\circ {\rm days}} \end{array}

where :math:`{GDD}_{10}` is the 20-year running mean growing
degree-days tracked from April through September (NH) base
10\ :math:`{}^\circ`\ C with maximum daily increments of
30\ :math:`{}^\circ`\ days. Eq. shows how we calculate
:math:`{GDD}_{0}`, :math:`{GDD}_{8}`, and :math:`{GDD}_{10}`:

.. math::
   :label: ZEqnNum977351 

   \begin{array}{l} {GDD_{{\rm 0}} =GDD_{0} +T_{2{\rm m}} -T_{f} {\rm \; \; \; where\; \; \; 0}\le T_{2{\rm m}} -T_{f} \le 26{}^\circ {\rm days}} \\ {GDD_{{\rm 8}} =GDD_{8} +T_{2{\rm m}} -T_{f} -8{\rm \; \; \; where\; \; \; 0}\le T_{2{\rm m}} -T_{f} -8\le 30{}^\circ {\rm days}} \\ {GDD_{{\rm 10}} =GDD_{10} +T_{2{\rm m}} -T_{f} -10{\rm \; \; \; where\; \; \; 0}\le T_{2{\rm m}} -T_{f} -10\le 30{}^\circ {\rm days}} \end{array}

where, if :math:`{T}_{2m}` -  :math:`{T}_{f}` takes on values
outside the above ranges, then it equals the minimum or maximum value in
the range. Also  :math:`{T}_{f}` equals 273.15 K,
:math:`{T}_{2m}` has units of K, and *GDD* has units of :sup:`o`\ days.

.. _Leaf emergence:

Leaf emergence
'''''''''''''''''''''''

According to AgroIBIS, leaves may emerge when the growing degree-days of
soil temperature to 0.05 m depth tracked since planting
(:math:`GDD_{T_{soi} }` ) reaches 3 to 5% of :math:`{GDD}_{mat}`
(Table 20.1). :math:`GDD_{T_{soi} }` is base 8, 0, and
10\ :math:`{}^\circ`\ C for corn, soybean, and temperate cereals. 
Leaf onset, as defined in the CN part of the model, occurs in the first
time step of phase 2, at which moment all seed C is transferred to leaf
C. Subsequently, the leaf area index generally increases and reaches
a maximum value during phase 2.

.. _Grain fill:

Grain fill
'''''''''''''''''''

Phase 3 begins in a similar way to phase 2. A variable tracked since
planting like :math:`GDD_{T_{soi} }`  but for 2-m air temperature,
:math:`GDD_{T_{{\rm 2m}} }`, must reach a heat unit threshold, *h*,
of 40 to 70% of  :math:`{GDD}_{mat}` (Table 20.1). For corn the
percentage itself is an empirical function of :math:`{GDD}_{mat}`
(not shown). In phase 3, the leaf area index begins to decline in
response to a background litterfall rate calculated as the inverse of
leaf longevity for the pft as done in the CN part of the model.

.. _Harvest:

Harvest
''''''''''''''''

Harvest is assumed to occur as soon as the crop reaches maturity. When
:math:`GDD_{T_{{\rm 2m}} }`  reaches 100% of :math:`{GDD}_{mat}` or
the number of days past planting reaches a crop-specific maximum (Table
20.1), then the crop is harvested. Harvest occurs in one time step using
CN’s leaf offset algorithm. New variables track the flow of grain C and
N to food and of live stem C and N to litter. Currently, food C and N
are routed directly to litter using the labile, cellulose, and lignin
fractions for leaves. The same fractions for leaves are used for the
flow of live stem C and N to litter for corn, soybean, and temperate
cereals. This is in contrast to the approach for unmanaged PFTs which
puts live stem C and N to dead stems first, rather than to litter.

.. _Allocation:

Allocation
^^^^^^^^^^^^^^^^^

Allocation responds to the same phases as phenology (section 20.2.3).
Simulated C assimilation begins every year upon leaf emergence in phase
2 and ends with harvest at the end of phase 3; therefore, so does the
allocation of such C to the crop’s leaf, live stem, fine root, and
reproductive pools.

.. _Leaf emergence to grain fill:

Leaf emergence to grain fill
'''''''''''''''''''''''''''''''''''''

During phase 2, the allocation coefficients (fraction of available C) to
each C pool are defined as:

.. math::
   :label: 25.4) 

   \begin{array}{l} {a_{repr} =0} \\ {a_{froot} =a_{froot}^{i} -(a_{froot}^{i} -a_{froot}^{f} )\frac{GDD_{T_{{\rm 2m}} } }{GDD_{{\rm mat}} } {\rm \; \; \; where\; \; \; }\frac{GDD_{T_{{\rm 2m}} } }{GDD_{{\rm mat}} } \le 1} \\ {a_{leaf} =(1-a_{froot} )\cdot \frac{a_{leaf}^{i} (e^{-b} -e^{-b\frac{GDD_{T_{{\rm 2m}} } }{h} } )}{e^{-b} -1} {\rm \; \; \; where\; \; \; }b=0.1} \\ {a_{livestem} =1-a_{repr} -a_{froot} -a_{leaf} } \end{array}

where :math:`a_{leaf}^{i}` , :math:`a_{froot}^{i}` , and
:math:`a_{froot}^{f}`  are initial and final values of these
coefficients (Table 20.2), and *h* is a heat unit threshold defined in
section 20.2.3. At a crop-specific maximum leaf area index,
:math:`{L}_{max}` (Table 20.2), carbon allocation is directed
exclusively to the fine roots.

.. _Grain fill to harvest:

Grain fill to harvest
''''''''''''''''''''''''''''''

The calculation of :math:`a_{froot}`  remains the same from phase 2 to
phase 3. Other allocation coefficients change to:

.. math::
   :label: ZEqnNum833921 

   \begin{array}{l} 
   {a_{leaf} =a_{leaf}^{i,3} {\rm \; \; \; when\; \; \; }a_{leaf}^{i,3} \le a_{leaf}^{f} {\rm \; \; \; else...}} \\ 
   {a_{leaf} =a_{leaf} \left(1-\frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \right)^{d_{alloc}^{leaf} } \ge a_{leaf}^{f} {\rm \; \; \; where\; \; \; }\frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \le 1} \\ 
   {} \\ 
   {a_{livestem} =a_{livestem}^{i,3} {\rm \; \; \; when\; \; \; }a_{livestem}^{i,3} \le a_{livestem}^{f} {\rm \; \; \; else...}} \\ 
   {a_{livestem} =a_{livestem} \left(1-\frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \right)^{d_{alloc}^{stem} } \ge a_{livestem}^{f} {\rm \; \; \; where\; \; \; }\frac{GDD_{T_{{\rm 2m}} } -h}{GDD_{{\rm mat}} d_{L} -h} \le 1} \\ 
   {} \\ 
   {a_{repr} =1-a_{froot} -a_{livestem} -a_{leaf} } 
   \end{array}

where :math:`a_{leaf}^{i,3}`  and :math:`a_{livestem}^{i,3}`  (initial
values) equal the last :math:`a_{leaf}`  and :math:`a_{livestem}` 
calculated in phase 2, :math:`d_{L}` , :math:`d_{alloc}^{leaf}`  and
:math:`d_{alloc}^{stem}`  are leaf area index and leaf and stem
allocation decline factors, and :math:`a_{leaf}^{f}`  and
:math:`a_{livestem}^{f}`  are final values of these allocation
coefficients (Table 20.2).

.. _General comments:

General comments
^^^^^^^^^^^^^^^^^^^^^^^

C and N accounting now includes new pools and fluxes pertaining to live
stems and reproductive tissues. For example, the calculations of growth
respiration, above ground net primary production, litter fall, and
displayed vegetation all now account for reproductive C.

We track allocation to reproductive C separately from CN’s allocation to
other C pools but within the CN framework. CN uses
:math:`{\textstyle\frac{a_{root} }{a_{leaf} }}`  and :math:`{\textstyle\frac{a_{livestem} }{a_{leaf} }}`  to calculate C and
N allometry and plant N demand.

Stem area index (*S*) is equal to 0.1\ *L* for corn and 0.2\ *L* for
other crops, as in AgroIBIS, where *L* is the leaf area index. All live
C and N pools go to 0 after crop harvest, but the *S* is kept at 0.25 to
simulate a post-harvest “stubble” on the ground.

Crop heights at the top and bottom of the canopy, :math:`{z}_{top}`
and :math:`{z}_{bot}` (m), come from the AgroIBIS formulation:

.. math::
   :label: 25.6) 

   \begin{array}{l} 
   {z_{top} =z_{top}^{\max } \left(\frac{L}{L_{\max } -1} \right)^{2} \ge 0.05{\rm \; where\; }\frac{L}{L_{\max } -1} \le 1} \\ 
   {z_{bot} =0.02{\rm m}} 
   \end{array}

The CN part of the model keeps track of a term representing excess
maintenance respiration that for perennial pfts or pfts with C storage
may be extracted from later gross primary production. Later extraction
cannot continue to happen after harvest for annual crops, so at harvest
we turn the excess respiration pool into a flux that extracts
CO\ :sub:`2` directly from the atmosphere. This way we eliminate
any excess maintenance respiration remaining at harvest as if such
respiration had not taken place.

In the list of plant physiological and other parameters used by the CLM,
we started the managed crops with the existing values assigned to the
unmanaged C3 crop. Then we changed the following parameters to
distinguish corn, soybean, and temperate cereals from the unmanaged C3
crop and from each other:

#. Growth respiration coefficient from 0.30 to the AgroIBIS value of
   0.25.

#. Fraction of leaf N in the Rubisco enzyme from 0.1 to 0.2 g N Rubisco
   g\ :sup:`-1` N leaf for temperate cereals to increase
   productivity (not chosen based on AgroIBIS).

#. Fraction of current photosynthesis displayed as growth changed from
   0.5 to 1 (not chosen based on AgroIBIS).

#. CLM4.5CN curve for the effect of temperature on photosynthesis
   instead of crop-specific curves from AgroIBIS.

#. Quantum efficiency at 25\ :sup:`o`\ C,
   :math:`\alpha` , from 0.06 to 0.04 *µ*\ mol CO\ :sub:`2`  *µ*\ mol\ :sup:`-1` photon for C4 crops (corn and unmanaged C4
   crop), using CLM4.5CN’s C4 grass value.

#. Slope, *m*, of conductance-to-photosynthesis relationship from 9 to 4 for C4 crops as in AgroIBIS.

#. Specific leaf areas, *SLA*, to the AgroIBIS values (Table 20.1).

#. Leaf orientation, :math:`\chi _{L}`, to the AgroIBIS values (Table 20.1).

#. Soil moisture photosynthesis limitation factor,
   :math:`\beta _{t}`, for soybeans multiplied as in AgroIBIS by 1.25
   for increased drought tolerance.

Table 20.1. Crop plant functional types (pfts) in CLM4.5CNcrop and their
parameters relating to phenology and morphology. Numbers in the first
column correspond to the list of pfts in Table 2.1.

+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
|       | Phenological                                 | :math:`T_{p}`        | :math:`T_{p}^{\min }`        | :math:`{GDD}_{min}`       | :math:`{GDD}_{mat}`       | Phase 2                      | Phase 3                      | Harvest: days   | :math:`z_{top}^{\max }`      | *SLA*                                           | :math:`\chi _{L}`      |    |
|       | Type                                         | K                    | K                            | ºdays                     | ºdays                     | %\ :math:`{GDD}_{mat}`       | %\ :math:`{GDD}_{mat}`       | past planting   | m                            | m\ :sup:`2`\ leaf g\ :sup:`-1`\ C               | index                  |    |
+=======+==============================================+======================+==============================+===========================+===========================+==============================+==============================+=================+==============================+=================================================+========================+====+
| 15.   | C\ :sub:`3` unmanaged rainfed crop           |                      |                              |                           |                           |                              |                              | 0.03            | -0.30                        |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 16.   | C\ :sub:`3` unmanaged irrigated crop         |                      |                              |                           |                           |                              |                              | 0.03            | -0.30                        |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 17.   | Rainfed Corn (also known as Maize)           |                      | 950-1850                     | 3                         | 55-65                     | :math:`\mathrm{\le}`\ 165    | 2.50                         | 0.05            | -0.50                        |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 18.   | Irrigated Corn (also known as Maize)         |                      | 950-1850                     | 3                         | 55-65                     | :math:`\mathrm{\le}`\ 165    | 2.50                         | 0.05            | -0.50                        |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 19.   | Rainfed Temperate Cereals                    |                      | :math:`\mathrm{\le}`\ 1700   | 5                         | 60                        | :math:`\mathrm{\le}`\ 150    | 1.20                         | 0.07            | 0.65                         |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 20.   | Irrigated Temperate Cereals                  |                      | :math:`\mathrm{\le}`\ 1700   | 5                         | 60                        | :math:`\mathrm{\le}`\ 150    | 1.20                         | 0.07            | 0.65                         |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 21.   | Rainfed Winter cereals (place holder)        |                      | 1900                         | 5                         | 40                        | :math:`\mathrm{\le}`\ 265    | 1.20                         | 0.07            | 0.65                         |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 22.   | Irrigated Winter cereals (place holder)      |                      | 1900                         | 5                         | 40                        | :math:`\mathrm{\le}`\ 265    | 1.20                         | 0.07            | 0.65                         |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 23.   | Rainfed Soybean                              |                      | :math:`\mathrm{\le}`\ 1700   | 3                         | 70                        | :math:`\mathrm{\le}`\ 150    | 0.75                         | 0.07            | -0.50                        |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+
| 24.   | Irrigated Soybean                            |                      | :math:`\mathrm{\le}`\ 1700   | 3                         | 70                        | :math:`\mathrm{\le}`\ 150    | 0.75                         | 0.07            | -0.50                        |                                                 |                        |    |
+-------+----------------------------------------------+----------------------+------------------------------+---------------------------+---------------------------+------------------------------+------------------------------+-----------------+------------------------------+-------------------------------------------------+------------------------+----+

Notes: :math:`T_{p}` and :math:`T_{p}^{\min }` are coldest
planting temperatures but for winter cereals :math:`T_{p}^{\min }`
is a warmest planting temperature. :math:`{GDD}_{min}` is the lowest
(for planting) 20-year running mean growing degree-days base 0ºC (winter
cereals) or 8 (other crops) tracked from April to September (NH).
:math:`{GDD}_{mat}` is a crop’s 20-year running mean growing
degree-days needed for vegetative and physiological maturity. Harvest
occurs at 100%\ :math:`{GDD}_{mat}` or when the days past planting
reach the number in the 10\ :sup:`th` column. Crop growth phases
are described in the text. :math:`z_{top}^{\max }`  is the maximum
top-of-canopy height of a crop, *SLA* is specific leaf area, and leaf
orientation index, :math:`\chi _{L}` , equals -1 for vertical, 0 for
random, and 1 for horizontal leaf orientation.

Table 20.2. Crop pfts in CLM4.5CNcrop and their parameters relating to
allocation. Numbers in the first column correspond to the list of pfts in Table 2.1.

+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
|       | :math:`a_{leaf}^{i}`                         | :math:`{L}_{max}`                        | :math:`a_{froot}^{i}`    | :math:`a_{froot}^{f}`    | :math:`a_{leaf}^{f}`    | :math:`a_{livestem}^{f}`    | :math:`d_{L}`        | :math:`d_{alloc}^{stem}`        | :math:`d_{alloc}^{leaf}`    |     |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
|       | fraction                                     | m\ :sup:`2`  m\ :sup:`-2`                |                          |                          |                         |                             |                      |                                 |                             |     |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 17.   | Rainfed Corn (also referred to as Maize)     | 0.800                                    | 5                        | 0.400                    | 0.050                   | 0.000                       | 0.000                | 1.05                            | 2                           | 5   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 18.   | Irrigated Corn (also referred to as Maize)   | 0.800                                    | 5                        | 0.400                    | 0.050                   | 0.000                       | 0.000                | 1.05                            | 2                           | 5   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 19.   | Rainfed Temperate Cereals                    | 0.750                                    | 7                        | 0.300                    | 0.000                   | 0.000                       | 0.050                | 1.05                            | 1                           | 3   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 20.   | Irrigated Temperate Cereals                  | 0.750                                    | 7                        | 0.300                    | 0.000                   | 0.000                       | 0.050                | 1.05                            | 1                           | 3   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 21.   | Rainfed Winter cereals (place holder)        | 0.425                                    | 7                        | 0.300                    | 0.000                   | 0.000                       | 0.050                | 1.05                            | 1                           | 3   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 22.   | Irrigated Winter cereals (place holder)      | 0.425                                    | 7                        | 0.300                    | 0.000                   | 0.000                       | 0.050                | 1.05                            | 1                           | 3   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 23.   | Rainfed Soybean                              | 0.850                                    | 6                        | 0.500                    | 0.200                   | 0.000                       | 0.300                | 1.05                            | 5                           | 2   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+
| 24.   | Irrigated Soybean                            | 0.850                                    | 6                        | 0.500                    | 0.200                   | 0.000                       | 0.300                | 1.05                            | 5                           | 2   |
+-------+----------------------------------------------+------------------------------------------+--------------------------+--------------------------+-------------------------+-----------------------------+----------------------+---------------------------------+-----------------------------+-----+

Notes: Crop growth phases and corresponding variables are described in
the text

.. _The irrigation model:

The irrigation model
-------------------------

The CLM includes the option to irrigate cropland areas that are equipped
for irrigation. The application of irrigation responds dynamically to
the soil moisture conditions simulated by the CLM. This irrigation
algorithm is based loosely on the implementation of Ozdogan et al.
(2010).

When irrigation is enabled, the crop areas of each grid cell are divided
into irrigated and rainfed fractions according to a dataset of areas
equipped for irrigation (Portmann et al. 2010). Irrigated and rainfed
crops are placed on separate soil columns, so that irrigation is only
applied to the soil beneath irrigated crops.

In irrigated croplands, a check is made once per day to determine
whether irrigation is required on that day. This check is made in the
first time step after 6 AM local time. Irrigation is required if crop
leaf area :math:`>` 0, and :math:`\beta_{t} < 1`, i.e., water is limiting for photosynthesis (see section 8.4).

If irrigation is required, the model computes the deficit between the
current soil moisture content and a target soil moisture content; this
deficit is the amount of water that will be added through irrigation.
The target soil moisture content in each soil layer *i*
(:math:`{s}_{target,i}`, kg m\ :sup:`-2`) is a weighted
average of the minimum soil moisture content that results in no water
stress in that layer (:math:`{s}_{o,i}`, kg m\ :sup:`-2`) and
the soil moisture content at saturation in that layer (:math:`{w}_{sat,i}`, kg m\ :sup:`-2`):

.. math::
   :label: 25.7) 

   w_{target,i} =(1-0.7)\cdot w_{o,i} +0.7\cdot w_{sat,i}

:math:`{w}_{o,i}` is determined by inverting equation 8.19 in Oleson
et al. (2010a) to solve for the value of :math:`{s}_{i}` (soil
wetness) that makes :math:`\Psi_{i} = \Psi_{o}` (where :math:`\Psi_{i}` is
the soil water matric potential and :math:`\Psi_{o}` is
the soil water potential when stomata are fully open), and then
converting this value to units of kg m\ :sup:`-2`.
:math:`{w}_{sat,i}` is calculated simply by converting effective
porosity (section 7.4) to units of kg m\ :sup:`-2`. The value 0.7
was determined empirically, in order to give global, annual irrigation
amounts that approximately match observed gross irrigation water use
around the year 2000 (i.e., total water withdrawals for irrigation:
:math:`\sim` 2500 – 3000 km\ :sup:`3` year\ :sup:`-1`
(Shiklomanov 2000)). The total water deficit (:math:`{w}_{deficit}`,
kg m\ :sup:`-2`) of the column is then determined by:

.. math::
   :label: 25.8) 

   w_{deficit} =\sum _{i}\max \left(w_{target,i} -w_{liq,i} ,0\right)

where  :math:`{w}_{liq,i}` (kg m\ :sup:`-2`) is the current
soil water content of layer *i* (Chapter 7). The max function means that
a surplus in one layer cannot make up for a deficit in another layer.
The sum is taken only over soil layers that contain roots. In addition,
if the temperature of any soil layer is below freezing, then the sum
only includes layers above the top-most frozen soil layer.

The amount of water added to this column through irrigation is then
equal to :math:`{w}_{deficit}`. This irrigation is applied at a
constant rate over the following four hours. Irrigation water is applied
directly to the ground surface, bypassing canopy interception (i.e.,
added to  :math:`{q}_{grnd,liq}`: section 7.1). Added irrigation is
removed from total liquid runoff ( :math:`{R}_{liq}`: Chapter 11),
simulating removal from nearby rivers.

.. _The details about what is new in CLM4.5:

The details about what is new in CLM4.5
--------------------------------------------

.. _Interactive irrigation for corn, temperate cereals, and soybean:

Interactive irrigation for corn, temperate cereals, and soybean
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CLM4.0 included interactive irrigation only for the generic C3 crops,
i.e. plant functional types (pfts) 15 (rainfed) and 16 (irrigated) in
the CLM list of pfts and not for the additional crops of the interactive
crop management model (CROP). Irrigation and CROP were mutually
exclusive in CLM4.0.

In CLM4.5 we have reversed this situation. Now the irrigation model can
be used only while running with CROP. To accomplish this we downloaded
data of percent irrigated and percent rainfed corn, soybean, and
temperate cereals (wheat, barley, and rye) (Portmann et al. 2010),
available online from

*ftp://ftp.rz.uni-frankfurt.de/pub/uni-frankfurt/physische\_geographie/hydrologie/public/data/MIRCA2000/harvested\_area\_grids.*

We embedded this data in CLM’s high-resolution pft data for use with the
tool mksurfdat to generate surface datasets at any desired resolution.
Now this data includes percent cover for 24 pfts:

1-16 as in the standard list of pfts, plus six more:

17 corn

18 irrigated\_corn

19 spring\_temperate\_cereal

20 irrigated\_spring\_temperate\_cereal

21 winter\_temperate\_cereal

22 irrigated\_winter\_temperate\_cereal

23 soybean

24 irrigated\_soybean

We intend surface datasets with 24 pfts only for CROP simulations with
or without irrigation. In simulations without irrigation, the rainfed
and irrigated crops merge into just rainfed crops at run time. Surface
datasets with 16 pfts can be used for all other CLM simulations.

.. _Interactive fertilization:

Interactive fertilization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

CLM adds nitrogen directly to the soil mineral nitrogen pool to meet
crop nitrogen demands. CLM’s separate crop land unit ensures that
natural vegetation will not access the fertilizer applied to crops.
Fertilizer amounts are obtained from the Agro-IBIS model (Kucharik and
Brye 2003), but can be modified in CLM’s pft-physiology input dataset.
Fertilizer is reported in g N/m\ :sup:`2` by plant functional
type. Total nitrogen fertilizer amounts are 150 g N/m\ :sup:`2`
for maize, 80 g N/m\ :sup:`2` for temperate cereals, and 25 g
N/m\ :sup:`2` for soybean, representative of central U.S. annual
fertilizer application amounts. Since CLM’s denitrification rate is high
and results in a 50% loss of the unused available nitrogen each day,
fertilizer is applied slowly to minimize the loss and maximize plant
uptake. Fertilizer application begins during the emergence phase of crop
development and continues for 20 days, which helps reduce large losses
of nitrogen from leaching and denitrification during the early stage of
crop development. The 20-day period is chosen as an optimization to
limit fertilizer application to the emergence stage. A fertilizer
counter in seconds, *f*, is set as soon as the onset growth for crops
initiates:

*f* = *n* \* 86400 [20.9)]

where *n* is set to 20 fertilizer application days. When the crop enters
phase 2 (leaf emergence to the beginning of grain fill) of its growth
cycle, fertilizer application begins by initializing fertilizer amount
to the total fertilizer divided by the initialized *f*. Fertilizer is
applied and *f* is decremented each time step until a zero balance on
the counter is reached.

The crop fertilization scheme was developed in versions of the CLM prior
to CLM4.5. In CLM4.5, crops with fertilization may be simulated over
productive.

.. _Biological nitrogen fixation for soybeans:

Biological nitrogen fixation for soybeans
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nitrogen fixation by soybeans is similar to that in the SWAT model
(Neitsch et al. 2005) and depends on soil moisture, nitrogen
availability, and growth stage. Soybean fixation is calculated only for
unmet nitrogen demand; if soil nitrogen meets soybean demand, there will
be no fixation during the time step. Soybean fixation is determined by

.. math::
   :label: 25.10) 

   N_{fix} \; =\; N_{plant\_ ndemand} \; *\; min\; \left(\; 1,\; fxw,\; fxn\; \right)*\; fxg

where :math:`{N}_{plant\_demand}` is the balance of nitrogen needed
to reach potential growth that cannot be supplied from the soil mineral
nitrogen pool, *fxw* is the soil water factor, *fxn* is the soil
nitrogen factor, and *fxg* is the growth stage factor calculated by

.. math::
   :label: 25.11) 

   fxw=\frac{wf}{0.85}

.. math::
   :label: 25.12) 

   fxn=\; \left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad {\rm for\; }sminn\le 10} \\ {1.5-0.005\left(sminn\times 10\right)\qquad {\rm for\; 10\; <\; }sminn{\rm \; }\ge 30} \\ {1\qquad \qquad \qquad \qquad {\rm for\; }sminn>30} \end{array}\right\}

.. math::
   :label: 25.13) 

   fxg=\left\{\begin{array}{l} {0\qquad \qquad \qquad \qquad \qquad {\rm for\; }GDD_{T_{2m} } \le 0.15} \\ {6.67\times GDD_{T_{2m} } -1\qquad \qquad \qquad {\rm for\; }0.15<GDD_{T_{2m} } \ge 0.30} \\ {1\qquad \qquad \qquad \qquad \qquad {\rm for\; }0.30<GDD_{T_{2m} } \ge 0.55} \\ {3.75-5\times GDD_{T_{2m} } \qquad \qquad \qquad {\rm for\; }0.55<GDD_{T_{2m} } \ge 0.75} \\ {0\qquad \qquad \qquad \qquad \qquad {\rm for\; }GDD_{T_{2m} } >0.75} \end{array}\right\}

where *wf* is the soil water content as a fraction of the water holding
capacity for the top 0.05 m, *sminn* is the total nitrogen in the soil
pool (g/m:sup:`2`), and :math:`{GDD}_{T_{2m}}` is the fraction of
growing degree-days accumulated during the growing season.
:math:`N\mathrm{fix}` is added directly to the soil mineral nitrogen
pool for use that time step. Nitrogen fixation occurs after the plant
has accumulated 15%\ :math:`{GDD}_{mat}` and before
75%\  :math:`{GDD}_{mat}`, so before grain fill begins.

.. _Modified C\:N ratios for crops:

Modified C:N ratios for crops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Typically, C:N ratios in plant tissue vary throughout the growing season
and tend to be lower during early growth stages and higher in later
growth stages. In order to account for this change, two sets of C:N
ratios are established in CLM for the leaf, stem, and fine root of
crops. This modified C:N ratio approach accounts for the nitrogen
retranslocation that occurs during phase 3 of crop growth. Leaf and stem
(and root for temperate cereals) C:N ratios for phases 1 and 2 are lower
than measurements (Table 20.3) to allow excess nitrogen storage in plant
tissue. During grain fill (phase 3) of the crop growth cycle, the
nitrogen in the plant tissues is moved to a storage pool to fulfill
nitrogen demands of organ (reproductive pool) development, such that the
resulting C:N ratio of the plant tissue is reflective of measurements at
harvest. All C:N ratios were determined by calibration process, through
comparisons of model output versus observations of plant carbon
throughout the growth season.

.. _Nitrogen retranslocation for crops:

Nitrogen retranslocation for crops
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Nitrogen retranslocation in crops occurs when nitrogen that was used for
tissue growth of leaves, stems, and fine roots during the early growth
season is remobilized and used for grain development (Pollmer et al.
1979; Crawford et al. 1982; Simpson et al. 1983; Ta and Weiland 1992;
Barbottin et al. 2005; Gallais et al. 2006, 2007). Nitrogen allocation
for crops follows that of natural vegetation, is supplied in CLM by the
soil mineral nitrogen pool, and depends on C:N ratios for leaves, stems,
roots, and organs. Nitrogen demand during organ development is fulfilled
through retranslocation from leaves, stems, and roots. Nitrogen
retranslocation is initiated at the beginning of the grain fill stage
for corn and temperate cereals, but not until after LAI decline in
soybean. Nitrogen stored in the leaf and stem is moved into a storage
retranslocation pool. For temperate cereals, nitrogen in roots is also
released into the retranslocation storage pool. The quantity of nitrogen
mobilized depends on the C:N ratio of the plant tissue, and is
calculated as

.. math::
   :label: 25.14) 

   leaf\_ to\_ retransn=\frac{c_{leaf} }{CN_{leaf} } -\frac{c_{leaf} }{CN_{leaf}^{f} }

.. math::
   :label: 25.15) 

   stemn\_ to\_ retransn=\frac{c_{stem} }{CN_{stem} } -\frac{c_{stem} }{CN_{stem}^{f} }

.. math::
   :label: 25.16) 

   frootn\_ to\_ retransn=\frac{c_{froot} }{CN_{froot} } -\frac{c_{froot} }{CN_{froot}^{f} }

where :math:`{C}_{leaf}`, :math:`{C}_{stem}`, and :math:`{C}_{froot}` is the carbon in the plant leaf, stem, and fine
root, respectively, :math:`{CN}_{leaf}`, :math:`{CN}_{stem}`, and :math:`{CN}_{froot}` is the pre-grain fill C:N ratio of the
leaf, stem, and fine root respectively, and :math:`CN^f_{leaf}`,
:math:`CN^f_{stem}`, and :math:`CN^f_{froot}` is the post-grain fill C:N
ratio of the leaf, stem, and fine root respectively (Table 20.3). Since
C:N measurements are taken from mature crops, pre-grain development C:N
ratios for leaves, stems, and roots are optimized to allow maximum
nitrogen accumulation for later use during organ development. Post-grain
fill C:N ratios are assigned the same as crop residue. Once excess
nitrogen is moved into the retranslocated pool, during the remainder of
the growing season the retranslocated pool is used first to meet plant
nitrogen demand by assigning the available nitrogen from the
retranslocated pool equal to the plant nitrogen demand. Once the
retranslocation pool is depleted, soil mineral nitrogen pool is used to
fulfill plant nitrogen demands.

Table 20.3. Pre- and post-grain fill C:N ratios for crop leaf, stem,
fine root, and reproductive pools.

+----------------------------+--------+---------------------+-----------+
| Pre-grain fill stage       | Corn   | Temperate Cereals   | Soybean   |
+============================+========+=====================+===========+
| :math:`{CN}_{leaf}`        | 10     | 15                  | 25        |
+----------------------------+--------+---------------------+-----------+
| :math:`{CN}_{stem}`        | 50     | 50                  | 50        |
+----------------------------+--------+---------------------+-----------+
| :math:`{CN}_{froot}`       | 42     | 30                  | 42        |
+----------------------------+--------+---------------------+-----------+
| Post-grain fill stage      |        |                     |           |
+----------------------------+--------+---------------------+-----------+
| :math:`CN_{leaf}^{f}`      | 65     | 65                  | 65        |
+----------------------------+--------+---------------------+-----------+
| :math:`CN_{stem}^{f}`      | 120    | 100                 | 130       |
+----------------------------+--------+---------------------+-----------+
| :math:`CN_{froot}^{f}`     | 42     | 40                  | 42        |
+----------------------------+--------+---------------------+-----------+
| :math:`CN_{repr}^{f}`      | 50     | 40                  | 60        |
+----------------------------+--------+---------------------+-----------+

.. _Separate reproductive pool:

Separate reproductive pool
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One notable difference between natural vegetation and crops is the
presence of a reproductive carbon pool (and nitrogen pool). Accounting
for the reproductive pool helps determine whether crops are performing
reasonably, through yield calculations, seasonal GPP and NEE changes,
etc. The reproductive pool is maintained similarly to the leaf, stem,
and fine root pools, but allocation of carbon and nitrogen does not
begin until the grain fill stage of crop development. Eq. shows the
carbon and nitrogen allocation coefficients to the reproductive pool. In
the CLM4.0, allocation of carbon to the reproductive pool was calculated
but merged with the stem pool. In the model, as allocation declines
during the grain fill stage of growth, increasing amounts of carbon and
nitrogen are available for grain development.
